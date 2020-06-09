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


class tbAna {
  public :
    tbAna();
    virtual ~tbAna(){;}
    TFile *inFile;
    TBaconRun *bRun;
    void calcResidual(int irun, TH1D* h, double start, double end);
    TH1D* hRes;
    TH1D* getLife(int run);
    TH1D* getMpv(int run);
    TH1D* hResMean;
    TH1D* hResSig;
    void histNorm(int norm, TH1D *hist);
};

void tbAna::histNorm(int norm, TH1D *hist){
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

TH1D* tbAna::getLife(int run) {
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

TH1D* tbAna::getMpv(int run) {
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



void tbAna::calcResidual(int irun, TH1D* h, double start, double end)
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

tbAna::tbAna() {
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


  TFile *fout = new TFile("tbAnaOut.root","RECREATE");

  TH1D* hMpvSet[MAXSETS];
  TH1D* hMpvSum[MAXSETS];
  TH1D* hLifeSum[MAXSETS];
  hResMean = new TH1D("ResMean","",100,-1.,1.);
  hResSig = new TH1D("ResSig","",100,-10,10.);


  for(int j=0; j<MAXSETS; ++j ) {
    hMpvSet[j]= (TH1D*) hMpv->Clone(Form("MpvSet-%s-%s",runTag[j].Data(),runRange[j].Data()) );
    hMpvSet[j]->SetTitle(Form("MpvSet %i %s %s",j+1,runTag[j].Data(),runRange[j].Data()) );
    hMpvSum[j]= (TH1D*) hChargeSum->Clone(Form("MpvSum-%s-%s",runTag[j].Data(),runRange[j].Data()) );
    hMpvSum[j]->SetTitle(Form("Mpv set %i Sum %s %s",j,runTag[j].Data(),runRange[j].Data()) );
    hLifeSum[j]= (TH1D*) hLife->Clone(Form("LifeSum-%s-%s",runTag[j].Data(),runRange[j].Data()) );
    hLifeSum[j]->SetTitle(Form("Lifetime set %i Sum %s %s",j+1,runTag[j].Data(),runRange[j].Data()) );
  }

  // scale errors and look at residuals
  bool doRes=false;
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1D *h = (TH1D*)key->ReadObj();
    TString hname(h->GetName());
    TString tnum = hname(7,hname.Length());
    int irun = atoi(tnum.Data());
    if(hname.Contains("LifeRun")) {
      // scale errors
      for(int ibin =0; ibin< h->GetNbinsX() ; ++ibin) {
        double e = h->GetBinError(ibin);
        h->SetBinError(ibin, e*3.4);
      }
      if(doRes) calcResidual(irun,h,7,8);
      if(doRes) calcResidual(irun,h,8,9);
    }
  }

  vector<double> vrun;
  vector<double> vrunErr;
  vector<double> vmpv;
  vector<double> vmpvErr;
  vector<double> vspe;
  vector<double> vspeErr;
  vrun.clear();
  vrunErr.clear();
  vmpv.clear();
  vmpvErr.clear();
  vspe.clear();
  vspeErr.clear();

  // fill vectors and collect sets
  int setNorm[MAXSETS]={0,0,0,0,0};
  TH1D* hadd=NULL;
  for(Long64_t entry =0; entry< treeRun->GetEntries(); ++entry) {
    treeRun->GetEntry(entry);
    int irun=bRun->run;
    int iset=-1;
    if(bRun->run>=3000&&bRun->run<=3020) {
      irun = irun - 3000 + 20000 - 20;
      iset=0;
    }
    if(bRun->run>=20005&&bRun->run<20020) iset=1;
    if(bRun->run>=20025&&bRun->run<=20040) iset=2;
    if(bRun->run>=20045&&bRun->run<=20060) iset=3;
    if(bRun->run>=20065&&bRun->run<=20080) iset=4;
    if(bRun->run>=20205&&bRun->run<=20215) iset=5;

    if(iset<0) {
      printf(" ..... skipping run  %i  \n",bRun->run);
      continue;
    }

    /*
    if(iset<MAXSETS-1&&bRun->run>runStart[iset]) { 
      ++iset;
      printf(" at run %i starting %i iset %i \n",bRun->run,runStart[iset],iset);
    }
    */

    printf(" set %i run %5i good events %5i MPV %5.1f +/- %5.1f ave SPE %5.1f +/- %5.1f\n",iset,
        bRun->run, bRun->nevc , bRun->mpvc, bRun->mpvcErr, bRun->aveSpe, bRun->aveSpeErr);
    vrun.push_back(double(irun));
    vrunErr.push_back(0.0);
    vmpv.push_back(bRun->mpvc);
    vmpvErr.push_back(bRun->mpvcErr);
    vspe.push_back(bRun->aveSpe);
    vspeErr.push_back(bRun->aveSpeErr);
    setNorm[iset]+= bRun->nevc;
    hMpvSet[iset]->Fill(bRun->mpvc);
    hadd = getMpv(bRun->run); 
    if(hadd) hMpvSum[iset]->Add(hadd);
    else printf("\t Mpv not found %i \n",bRun->run);
    hadd = getLife(bRun->run); 
    if(hadd) hLifeSum[iset]->Add(hadd);
    else printf("\t Life not found %i \n",bRun->run);
  }



  TGraphErrors *gmpv = new TGraphErrors(vrun.size(),&vrun[0],&vmpv[0],&vrunErr[0],&vmpvErr[0]);
  TCanvas *cmpv = new TCanvas("mpv","mpv");
  cmpv->SetGridx(); cmpv->SetGridy();
  gmpv->SetTitle("Landau MPV");
  gmpv->SetMarkerColor(kBlue);
  gmpv->SetMarkerStyle(22);
  gmpv->SetMarkerSize(1);
  gmpv->GetXaxis()->SetTitle(" run ");
  gmpv->GetYaxis()->SetTitle(" Landau MPV ");
  gmpv->Draw("ap");
  cmpv->Print(".png");


  TGraphErrors *gspe = new TGraphErrors(vrun.size(),&vrun[0],&vspe[0],&vrunErr[0],&vspeErr[0]);
  TCanvas *cspe = new TCanvas("spe","spe");
  cspe->SetGridx(); cspe->SetGridy();
  gspe->SetTitle("Ave SPE");
  gspe->SetMarkerColor(kBlue);
  gspe->SetMarkerStyle(22);
  gspe->SetMarkerSize(1);
  gspe->GetXaxis()->SetTitle(" run ");
  gspe->GetYaxis()->SetTitle(" ave SPE ");
  gspe->Draw("ap");
  cspe->Print(".png");


  TGraphErrors *gcomp = new TGraphErrors(vmpv.size(),&vmpv[0],&vspe[0],&vmpvErr[0],&vspeErr[0]);
  TCanvas *ccomp = new TCanvas("comp","comp");
  ccomp->SetGridx(); ccomp->SetGridy();
  gcomp->SetTitle("Ave SPE");
  gcomp->SetMarkerColor(kBlue);
  gcomp->SetMarkerStyle(22);
  gcomp->SetMarkerSize(1);
  gcomp->GetXaxis()->SetTitle("  MPV ");
  gcomp->GetYaxis()->SetTitle(" ave SPE ");
  gcomp->Draw("ap");
  ccomp->Print(".png");


  // normalize
  printf(" normalize \n");
  for(int iset=0; iset<MAXSETS; ++iset) {
    double mpvTot = hMpvSum[iset]->Integral();
    histNorm(setNorm[iset],hMpvSum[iset]);
    histNorm(setNorm[iset],hLifeSum[iset]);
    printf("\t  set %i %s %s %i mpv %f (%f) \n",iset, runTag[iset].Data(),runRange[iset].Data(),setNorm[iset],mpvTot,hMpvSum[iset]->Integral());
  }

  /*
  printf("look for neg bins \n");
  for(int iset=0; iset<MAXSETS; ++iset) {
   for(int ibin =0; ibin< hLifeSum[iset]->GetNbinsX() ; ++ibin) 
      if(hLifeSum[iset]->GetBinContent(ibin)<=0) printf(" %s %i %E \n",hLifeSum[iset]->GetName(),ibin,hLifeSum[iset]->GetBinContent(ibin));
  }
  */


  fout->ls();
  fout->Write();


  //deriv

  /*
     vector<double> vderiv;
  vector<double> vderivErr;

  for(unsigned i=0; i<vmpv.size()-1; ++i) {
    vderiv.push_back(vmpv[i]-vmpv[i+1]);
    vderivErr.push_back( sqrt( pow(vmpvErr[i+1],2) + pow(vmpvErr[i],2) ));
  }


  
  TGraphErrors *gderiv = new TGraphErrors(vderiv.size(),&vrun[0],&vderiv[0],&vrunErr[0],&vderivErr[0]);
  TCanvas *cderiv = new TCanvas("deriv","deriv");
  cderiv->SetGridx(); cderiv->SetGridy();
  gderiv->SetTitle("MPV derivative ");
  gderiv->SetMarkerColor(kBlue);
  gderiv->SetMarkerStyle(22);
  gderiv->SetMarkerSize(1);
  gderiv->GetXaxis()->SetTitle(" run ");
  gderiv->GetYaxis()->SetTitle(" MPV derivative  ");
  gderiv->Draw("ap");
  cderiv->Print(".png");
  */
  
}

