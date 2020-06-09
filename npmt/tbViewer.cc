/* uses TBaconRun Class */
//////////////////////////////////////////////////////////
//  M.Gold May 2020 
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex>//includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"

#include "TBaconRun.hxx"

/*
 0 PPM 10000_10100
 1 PPM 20000_20021
 2 PPM 20022_20041
 5 PPM 20042_20062
10 PPM 20063_29999
*/



static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double binwidth = par[2];
  double tau = par[0];
  double f = par[1]/tau*binwidth*TMath::Exp(-x/tau);
  return f; 
}


class tbViewer {
  public :
    tbViewer();
    virtual ~tbViewer(){;}
    TFile *inFile;
    TBaconRun *bRun;
    TH1D* getLife(int run);
    TH1D* getMpv(int run);
    void histNorm(int norm, TH1D *hist);
    TH1D* hRes;
    TH1D* hResMean;
    TH1D* hResSig;
    void calcResidual(int irun, TH1D* h, double start, double end);
};

void tbViewer::histNorm(int norm, TH1D *hist){
  for(int ibin =0; ibin< hist->GetNbinsX() ; ++ibin) {
    double dnorm = double(norm);
    double cbin = hist->GetBinContent(ibin);
    double c,e;
    if(cbin==0) { 
      c=1E-9;
      e=0;
    } else {
      double ebin = hist->GetBinError(ibin);
      c = cbin/dnorm;
      e = c*sqrt( pow(ebin/cbin,2.)+1./dnorm);
      if(c-e<0) e=c;
    }
    hist->SetBinContent(ibin,c);
    hist->SetBinError(ibin,e);
  }
}

TH1D* tbViewer::getLife(int run) {
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  TH1D* hlife = NULL;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1D *h = (TH1D*)key->ReadObj();
    TString hname(h->GetName());
    if(hname.Contains("LifeRun")) {
      TString tnum = hname(7,hname.Length());
      int irun = atoi(tnum.Data());
      if(irun==run ) hlife=h;
    }
  }
  return hlife;
}

TH1D* tbViewer::getMpv(int run) {
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  TH1D* hmvd = NULL;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1D *h = (TH1D*)key->ReadObj();
    TString hname(h->GetName());
    if(hname.Contains("ChargeSumRun")){ 
      TString tnum = hname(12,hname.Length());
      int irun = atoi(tnum.Data());
      if(irun==run) hmvd=h;
      if(hmvd) break;
    }
  }
  return hmvd;
}



void tbViewer::calcResidual(int irun, TH1D* h, double start, double end)
{
  Double_t binwidth = h->GetBinWidth(1);
  TF1* f1 = new TF1(Form("fexp-%i",irun),fexp,start,end,3);
  cout << "  xxxx " << h->GetName() << "  " << f1->GetName() << endl;
  f1->SetNpx(1000); // numb points for function
  f1->SetParNames("lifetime","integral","binwidth");
  f1->SetParameters(1,1,binwidth);
  f1->FixParameter(2,binwidth);
  h->Fit(f1,"R0");
  TF1 *func = h->GetFunction(f1->GetName());
  double mean=0;
  double meanErr=0;
  if(func) {
    mean = func->GetParameter(0);
    meanErr = func->GetParError(0);
  }

 
  TH1D *hResRun = (TH1D*) hRes->Clone(Form("hResRun-%i-from-%i-to-%i",irun,int(start),int(end)));

  int ibinStart = h->FindBin(start);
  int ibinEnd = h->FindBin(end);
  for(int ibin=ibinStart; ibin< ibinEnd; ++ibin ) {
    double t = h->GetBinLowEdge(ibin);
    double r= (h->GetBinContent(ibin) - func->Eval(t))/ h->GetBinError(ibin);
    hResRun->Fill(r);
    //printf(" irun %i ibin %i res %f\n",irun,ibin,r);
  }
  gStyle->SetOptFit();
  hResRun->Fit("gaus","0");
  TF1 *fgaus = hResRun->GetFunction("gaus");

  double gmean=0;
  double gmeanErr=0;
  double gsig=0;
  double gsigErr=0;

  if(fgaus) {
    gmean = fgaus->GetParameter(1);
    gmeanErr = fgaus->GetParError(1);
    gsig = fgaus->GetParameter(2);
    gsigErr = fgaus->GetParError(2);

  }
  hResMean->Fill(fgaus->GetParameter(1));
  hResSig->Fill(fgaus->GetParameter(2));
  printf(" Exp fit run %i range (%.0f,%.0f) exp fit mean %f +/- %f  residual %f +/- %f  sig %f +/- %f  \n",irun, start,end, mean, meanErr,gmean,gmeanErr,gsig,gsigErr);


  TCanvas *can = new TCanvas(Form("ResFitRun-%i-from-%i-to-%i",irun,int(start),int(end)),"");
  hResRun->Draw();
  can->Print(".pdf");

}

tbViewer::tbViewer() {
  /* data sets */
  double maxLife=10.0;
  int lifeBins = int(1.0E4/8.0); // ns bins

  enum {MAXSETS=6};
  
  TString runTag[MAXSETS];
  runTag[0]=TString("00PPM-MU");
  runTag[1]=TString("01PPM-Ran");
  runTag[2]=TString("02PPM-Ran");
  runTag[3]=TString("05PPM-Ran");
  runTag[4]=TString("10PPM-Ran");
  runTag[5]=TString("10PPM-MU");

  TString runRange[MAXSETS];
  runRange[0]=TString("20000-20000");
  runRange[1]=TString("20001-20020");
  runRange[2]=TString("20022-20040");
  runRange[3]=TString("20040-20060");
  runRange[4]=TString("20060-20080");
  runRange[5]=TString("20200-20215");

   // prototypes 
  hRes = new TH1D("hres","",40,-10,10);

  inFile = new TFile("tbReader-3000-29999-1250-bins.root");
  if(!inFile) return;
  TTree *treeRun = NULL;
  inFile->GetObject("tRun",treeRun);
  if(!treeRun) return;

  bRun = new TBaconRun();
  treeRun->SetBranchAddress("brun",&bRun);

  printf(" treeRun has %lld entries \n",treeRun->GetEntries());

  // prototypes 
  hRes = new TH1D("hres","",40,-10,10);
  TH1D* hMpv = new TH1D("hmpv","",500,200,1200);
  TH1D* hChargeSum = new TH1D("ChargeSum"," good summed event charge > 10 pulses  ",1000,0,10000);
  hChargeSum->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");
  TH1D *hLife = new TH1D("LifeCut"," lifetime PMT >10 hits ",lifeBins,0,maxLife);
  hLife->GetXaxis()->SetTitle(" micro-seconds ");
  hLife->SetMarkerColor(kBlack);
  hLife->SetMarkerStyle(22);
  hLife->SetMarkerSize(.2);


  hResMean = new TH1D("ResMean","",100,-1.,1.);
  hResSig = new TH1D("ResSig","",100,-10,10.);

  enum {MAXHIST=100};
  TH1D* hLifeSum[MAXHIST];
  TH1D* hMpvSum[MAXHIST];
  int runNum[MAXHIST];
  int nhist=0;
  int nhist2=0;


  // scale errors and look at residuals
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1D *h = (TH1D*)key->ReadObj();
    TString hname(h->GetName());
    TString tnum = hname(7,hname.Length());
    TString tnum2 = hname(12,hname.Length());
    int irun = atoi(tnum.Data());
    int irun2 = atoi(tnum2.Data());
    if(hname.Contains("LifeRun")&&irun<4000) {
      hLifeSum[nhist]=h;
      runNum[nhist++]=irun;
    }
    if(hname.Contains("ChargeCutRun")&&irun2<4000)  {
      //cout << irun2 << "  " << hname << endl;
      hMpvSum[nhist2++]=h;
    }
  }

  cout << " got " << nhist << " " << nhist2 << endl;
  if(nhist!=nhist2) return;

  TCanvas *cjunk = new TCanvas("junk","junk");

  

  bool doRes=false;
  // scale errors
  for(int iset=1; iset<nhist; ++iset) {
    for(int ibin =0; ibin<  hLifeSum[iset]->GetNbinsX() ; ++ibin) {
      double e =  hLifeSum[iset]->GetBinError(ibin);
      hLifeSum[iset]->SetBinError(ibin, e*3.4);
    }
    if(doRes) calcResidual(runNum[iset],hLifeSum[iset],7,8);
    //if(doRes) calcResidual(irun,h,8,9);
  }

  TFile *fout = new TFile("tbViewerOut.root","RECREATE");

  TH1D* hLifeAll = (TH1D*) hLifeSum[0]->Clone();
  hLifeAll->Reset();

  for(int iset=0; iset<nhist; ++iset) {
    int norm = int( hMpvSum[iset]->GetEntries());
    double preIntegral = hLifeSum[iset]->Integral();
    histNorm(norm,hLifeSum[iset]);
    printf(" %i %s %s norm  %i pre  integral %f norm integral %f \n",runNum[iset],hLifeSum[iset]->GetName(), hMpvSum[iset]->GetName(),norm,preIntegral,hLifeSum[iset]->Integral() );
    fout->Append(hLifeSum[iset]);
    hLifeAll->Add(hLifeSum[iset]);
  }

  histNorm(nhist,hLifeAll);

  TCanvas *canLifeAll=new TCanvas("Life-ALL-set1","Life-ALL-set1");
  canLifeAll->SetLogy();
  hLifeAll->SetMarkerStyle(24);
  hLifeAll->SetMarkerColor(kRed);
  hLifeAll->SetLineColor(kRed);
  hLifeAll->SetMarkerSize(0.7);

  hLifeAll->Draw();
  hLifeAll->SetTitle("SPE versus Time  ");
  for(int iset=0; iset<nhist; ++iset) {
    hLifeSum[iset]->SetTitle("SPE versus Time  ");
    hLifeSum[iset]->Draw("sames");
  }
  canLifeAll->Print(".pdf");



 
  fout->ls();
  fout->Write();


}

