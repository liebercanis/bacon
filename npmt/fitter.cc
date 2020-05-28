/* uses TBaconEvent class */
/*#include "fitter.hh"*/
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

//
#include "TBaconEvent.hxx"
using namespace TMath;

static double xstart =1.1;
static double xstop =10;
static double tiny =.000021; 
/*I fix all of the lifetimes: microseconds*/
static double pmodtS = 7E-3;
static double pmodtT = 1.96;
static double pmodtG = 3.2;
static double pmodtX = 20E-3;

// 9 parameter model !
static double pmodel(double *xx, double *par) 
{
  int type = int(par[9]);
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double S0 = binwidth*par[1];
  double T0 = binwidth*par[2];
  double G0 = binwidth*par[3];
  double ts = par[4];
  double tt = par[5];
  double k1 = par[6];
  double tg = par[7];
  double tx = par[8];
  double S = S0*Exp(-t*(1.0/ts+k1));
  double T = T0*Exp(-t*(1.0/tt+k1));
  double G = G0*Exp(-t*(1.0/tg));
  //
  double M = k1*(S0*ts+T0*tt)*Exp(-t*k1) -k1*(S*ts+T*tt);
  //
  double x1 = tx-tt+k1*tt*tx;
  double x2 = tx-tt+k1*ts*tx;
  double x3 = k1*tx-1.0;
  //
  double y1= (S0*ts*x1+T0*tt*x2)*k1*k1*tx*tx/x1/x2/x3;
  double y2=k1*k1*ts*ts*tx/x2;
  double y3=k1*k1*tt*tt*tx/x1;;
  double y4=k1*k1*tx*(S0*ts+T0*tt)/x3;
  double X = y1*Exp(-t/tx) + y2*S + y3*T - y4*Exp(-t*k1);
 //
  double f;
  if(type ==0) f=T+S;
  else if(type ==1) f=M;
  else if(type ==2) f=X;
  else if(type ==3) f=G;
  else f =  S+T+G+X;  
  return f;
}


static double kvalFit(double x) {
  double k0 = 3.53;  // offset for undoped argon
  return Max(1E-3,0.196*(x+k0)+0.689);
}

static void ratioE(double x, double xe, double y, double ye, double& r,double& re){
  r = x/y;
  re = r * sqrt( pow(xe/x,2.)+pow(ye/y,2.) );
}

// fmodelFit[ifit]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
static double fmodel(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  int type = int(par[7]);
  if(t<0) return tiny;
  double binwidth = par[0];
  double A0 = binwidth*par[1];
  double ta = par[2];
  double k = par[3];
  double tx = par[4];
  double Am = binwidth*par[5];
  double tm = par[6];
  double tap = 1./( 1./ta + k);
  double t1=1./( 1./tx - k);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta*k;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t*k) - Exp(-t/tap) );
  double X = C*k*(t1*Exp(-t*k) - t2*Exp(-t/tap)) + C*k*(t2-t1)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f;
  if(type ==0) f=A;
  else if(type ==1) f=D;
  else if(type ==2) f=X;
  else if(type ==3) f=M;
  else f =  A+X+M; 
  return f;
}



static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double binwidth = par[2];
  double tau = par[0];
  double f = par[1]/tau*binwidth*Exp(-x/tau);
  return f; 
}


static double ftriple(double *xx, double *par)
{
  double x= xx[0];
  double bw = par[0];
  double A1 = par[1];
  double A2 = par[2];
  double A3 = par[3];
  double ts = par[4];
  double td = par[5];
  double tl = par[6];
  double f =  bw*( A1/ts*Exp(-x/ts) - A2/td*Exp(-x/td) + A3/tl*Exp(-x/tl) );
  return f; 
}


class fitter {
  public :
    fitter(Long64_t nloop=0);
    virtual ~fitter(){;}
    enum {MAXFILE=6};
    enum {PFITS=6};
    TF1 *fp[MAXFILE];
    TH1D* hLife[MAXFILE];
    TH1D* hQSum[MAXFILE];
    TH1D* hQSum2[MAXFILE];
};


fitter::fitter(Long64_t nloop)
{
  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    hLife[ifile]=NULL;
    fp[ifile]=NULL;
  }

  TString tag("DS_2");

  TString runTag[MAXFILE+1];
  runTag[0]=TString("00_PPM");
  runTag[1]=TString("01_PPM");
  runTag[2]=TString("02_PPM");
  runTag[3]=TString("05_PPM");
  runTag[4]=TString("10_PPM");
  runTag[5]=TString("10_PPMN2");
  runTag[6]=TString("10_PPMN10");

  TString runRange[MAXFILE+1];
  runRange[0]=TString("10000_20000");
  runRange[1]=TString("20001_20021");
  runRange[2]=TString("20022_20041");
  runRange[3]=TString("20042_20062");
  runRange[4]=TString("20063_29999");
  runRange[5]=TString("30000_40001");
  runRange[6]=TString("30100_40001");

  double runLow[MAXFILE+1];
  double runHigh[MAXFILE+1];

  runLow[0]=10000;
  runLow[1]=20001;//1
  runLow[2]=20022;//2
  runLow[3]=20042;//5
  runLow[4]=20061;//10
  runLow[5]=20081;//11
  runLow[6]=20219;//11

  runHigh[0]=20000;
  runHigh[1]=20021;
  runHigh[2]=20041;
  runHigh[3]=20060;
  runHigh[4]=20080;
  runHigh[5]=20100;
  runHigh[6]=40001;

  
  TFile *fin[MAXFILE];

  TString fname[MAXFILE];
  TString trigTag("All");
  for(int ifile=0; ifile<MAXFILE+1; ++ifile) {
    //fname[ifile]=TString(Form("life-%s_DS_2-%s-bins-1250.root",runRange[ifile].Data(),trigTag.Data()));
    fname[ifile]=TString(Form("life-%s_DS_2-%s-bins-400.root",runRange[ifile].Data(),trigTag.Data()));
    cout << ifile << ")  " << "  " <<fname[ifile] << endl;
  }

  TH1D* hist;
  TH1D* hist2;
  TH1D* hist3;
  TNtuple *nt;
  TList* ntList = new TList;
  gStyle->SetOptStat(0);
  for(int ifile =0; ifile<MAXFILE; ++ifile ) {
    fin[ifile] = new TFile(fname[ifile],"READONLY");
    fin[ifile]->GetObject("LifeAll",hist);
    fin[ifile]->GetObject("EventSumCut",hist2);
    fin[ifile]->GetObject("EventSumNeil",hist3);
    fin[ifile]->GetObject("ntRun",nt);
    ntList->Add(nt);
    if(!hist) continue;
    cout << "file " << fin[ifile]->GetName() << "  " << hist->GetName() << " " << hist2->GetName() << " " << hist3->GetName() <<endl;
    hLife[ifile]= (TH1D*) hist->Clone(Form("LifeAll%i",ifile));
    hQSum[ifile]= (TH1D*) hist2->Clone(Form("EventQ%i",ifile));
    hQSum2[ifile]= (TH1D*) hist3->Clone(Form("EventQNeil%i",ifile));

    TObject * fobj= hLife[ifile]->GetListOfFunctions()->FindObject("modelFit0");
    hLife[ifile]->GetListOfFunctions()->ls();
    if(fobj) cout << " remove  fobj " << fobj->GetName() << endl;
    hLife[ifile]->GetListOfFunctions()->RecursiveRemove(fobj);
    hLife[ifile]->SetTitle(Form("life-%s-pass-%s",runTag[ifile].Data(),trigTag.Data()));
    hLife[ifile]->SetStats(0);

    double total =  hLife[ifile]->Integral();
    printf(" file %i total %f integral %f bw %f \n",ifile,total,hLife[ifile]->Integral(), hLife[ifile]->GetBinWidth(1));
  }

  for(int ifile=0; ifile<MAXFILE; ++ifile) if(!hLife[ifile]) return;

  TTree * tRunSum = TTree::MergeTrees(ntList);
  tRunSum->GetListOfBranches()->ls();
  cout<< " tRunSum entries " << tRunSum->GetEntries() << endl;
  Float_t frun,fnev,fmean,fmeanErr,fsig,fsigErr,fmpv,fmpvErr,fwid,fwidErr,fSpe;
  tRunSum->SetBranchAddress("run",&frun);
  tRunSum->SetBranchAddress("nev",&fnev);
  tRunSum->SetBranchAddress("mean",&fmean);
  tRunSum->SetBranchAddress("meanErr",&fmeanErr);
  tRunSum->SetBranchAddress("sig",&fsig);
  tRunSum->SetBranchAddress("sigErr",&fsigErr);
  tRunSum->SetBranchAddress("mpv",&fmpv);
  tRunSum->SetBranchAddress("mpvErr",&fmpvErr);
  tRunSum->SetBranchAddress("wid",&fwid);
  tRunSum->SetBranchAddress("widErr",&fwidErr);
  tRunSum->SetBranchAddress("sumSPE",&fSpe);
  TFile *fout = new TFile(Form("fitter-%s-%s.root",tag.Data(),trigTag.Data()),"RECREATE");

  TTree *treeOut = (TTree*) tRunSum->Clone("runSummary");
  fout->Append(treeOut);


  TH1D* hSumMean[MAXFILE+1];
  TH1D* hSumMpv[MAXFILE+1];
  for(unsigned irun=0; irun<MAXFILE+1; ++irun) { 
    hSumMean[irun] = new TH1D(Form("SumMean%i",irun),Form("sunMean %f to %f ",runLow[irun],runHigh[irun]),200,200,1200);
    hSumMpv[irun] = new TH1D(Form("SumMpv%i",irun),Form("sunMpv %f to %f ",runLow[irun],runHigh[irun]),200,200,1200);
  }




  vector<float> vrun;
  vector<float> vrunErr;
  vector<float> vmean;
  vector<float> vmeanErr;
  vector<float> vmpv;
  vector<float> vmpvErr;
  vector<float> vSpe;
  
  vector<float> vrunSum;
  vector<float> vmeanSum;
  vector<float> vmpvSum;
  vector<float> vnumSum;
  vrunSum.push_back(0);
  vrunSum.push_back(1);
  vrunSum.push_back(2);
  vrunSum.push_back(5);
  vrunSum.push_back(10);
  vrunSum.push_back(11);
  vrunSum.push_back(12);


  float meanSum=0;
  float mpvSum=0;
  bool isNew = false;
  unsigned numSum=0;
  unsigned bcount=0;
  for(Long64_t entry=0; entry< tRunSum->GetEntries() ; ++entry) {
    tRunSum->GetEntry(entry);
    if(frun<20000) vrun.push_back(frun+9990);
    else if(frun>=30000) vrun.push_back(frun-30000+20100);
    else  vrun.push_back(frun);
    printf("..... %lu %f %f count %u %f %f \n",vrun.size(),frun,vrun[vrun.size()-1],bcount,runLow[bcount],runHigh[bcount]);
    vrunErr.push_back(0);
    vmean.push_back(fmean);
    vmeanErr.push_back(fmeanErr);
    vmpv.push_back(fmpv);
    vmpvErr.push_back(fmpvErr);
    meanSum+= fmean;
    mpvSum+= fmpv;
    ++numSum;
    if(fnev>1) vSpe.push_back(fSpe/fnev);
    else vSpe.push_back(0);
    
    hSumMean[bcount]->Fill(fmean);
    hSumMpv[bcount]->Fill(fmpv);

    // break for runs 
    if(frun>=runLow[bcount+1]) {
      printf("xxxxxxxxxxxxxxxxx %u break at %f (%f,%f) %u %f %f  %f \n",bcount,frun,runLow[bcount],runHigh[bcount], numSum, meanSum, mpvSum,
          hSumMean[bcount]->GetEntries());
      vmeanSum.push_back(meanSum);
      vmpvSum.push_back(mpvSum);
      vnumSum.push_back(float(numSum));
      meanSum=0;
      mpvSum=0;
      numSum=0;
      ++bcount;
    }
  }
  printf("LLLxxxxxxxxxxxxxx %u break at %f (%f,%f) %u %f %f %f \n",bcount,frun,runLow[bcount],runHigh[bcount], numSum, meanSum, mpvSum,
        hSumMean[bcount]->GetEntries());
  vmeanSum.push_back(meanSum);
  vmpvSum.push_back(mpvSum);


  
  return;

  TGraphErrors *gLanMean = new TGraphErrors(vrun.size(), &vrun[0],&vmean[0],&vrunErr[0],&vmeanErr[0]);
  TGraphErrors *gLanMpv = new TGraphErrors(vrun.size(), &vrun[0],&vmpv[0],&vrunErr[0],&vmpvErr[0]);
  TGraph *gSpe = new TGraph(vrun.size(), &vrun[0],&vSpe[0]);
  TGraph *gComp = new TGraph(vmpv.size(), &vmpv[0],&vmean[0]);

  TCanvas *cLanSig = new TCanvas("mpvByRun","mpvByRun");
  cLanSig->SetGridx(); cLanSig->SetGridy();
  gLanMean->GetHistogram()->GetYaxis()->SetRangeUser(100,1200);
  gLanMean->SetName("Landau-sigByRun");
  gLanMean->SetTitle("Landau-sigByRun");
  gLanMean->SetMarkerColor(kRed);
  gLanMean->SetMarkerStyle(22);
  gLanMean->SetMarkerSize(1);
  gLanMean->GetXaxis()->SetTitle(" Xe ppm ");

  gLanMpv->SetTitle("Landau-mpvByRun");
  gLanMpv->SetName("Landau-mpvByRun");
  gLanMpv->SetMarkerColor(kBlue);
  gLanMpv->SetMarkerStyle(21);
  gLanMpv->SetMarkerSize(1);
  gLanMpv->GetXaxis()->SetTitle(" Xe ppm ");

  gLanMean->Draw("ap");
  gLanMpv->Draw("psame");
  cLanSig->Print(".png");

  TCanvas *cSpe = new TCanvas("speByRun","speByRun");
  cSpe->SetGridx(); cSpe->SetGridy();
  gSpe->SetName("SPE-ByRun");
  gSpe->SetTitle("mean SPE By Run");
  gSpe->SetMarkerColor(kRed);
  gSpe->SetMarkerStyle(21);
  gSpe->SetMarkerSize(1);

  gSpe->Draw("ap");
  cSpe->Print(".png");

  // comp plot
  
  TCanvas *cComp = new TCanvas("compByRun","compByRun");
  cComp->SetGridx(); cLanSig->SetGridy();
  gComp->GetHistogram()->GetYaxis()->SetRangeUser(500,1200);
  gComp->SetName("LandauCompByRun");
  gComp->SetTitle("Landau Comp By Run");
  gComp->SetMarkerColor(kRed);
  gComp->SetMarkerStyle(22);
  gComp->SetMarkerSize(1);
  gComp->GetXaxis()->SetTitle(" MPV ");
  gComp->GetYaxis()->SetTitle(" mean ");
  gComp->Draw("ap");
  cComp->Print(".png");


  fout->Append(gLanMpv);
  fout->Append(gLanMean);
  fout->Append(gSpe);

  printf("\t\t summary \n");
  for(unsigned ir=0; ir<vrun.size(); ++ir ) printf(" run %.f mean %f mpv %f spe %f \n",  vrun[ir],vmean[ir], vmpv[ir], vSpe[ir]);

  return;
  

  Double_t binwidth = hLife[0]->GetBinWidth(1);
  Double_t integral = hLife[0]->Integral();

  // pmodel 
  TF1* pmodelFit[MAXFILE][PFITS];
  int pcolor[PFITS] = {2,3,4,45,1,30};
  int pstyle[PFITS] = {21,22,23,24,25,26};
  double runPPM[MAXFILE] = {0,1,2,5,10,11};
  double Single0=0.3;
  double Triplet0=0.7*integral;
  double Gas0=.001*integral;
  //TCanvas *canPmodel[MAXFILE];

  // fmodel 
  double tauX=pmodtX;
  double tauA = pmodtT;
  double tauM = pmodtG; 
  double Amystery = Gas0;
  TF1* fmodelFit[MAXFILE][PFITS];
  TCanvas *canFmodel[MAXFILE];


  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    double k1= kvalFit(runPPM[ifile]);
    for(int ifit=0; ifit<PFITS; ++ifit) {
      pmodelFit[ifile][ifit] = new TF1(Form("pmodelFit-%.f-%i",runPPM[ifile],ifit),pmodel,xstart,xstop,10);
      pmodelFit[ifile][ifit]->SetNpx(1000); // numb points for function
      pmodelFit[ifile][ifit]->SetParNames("binwidth","S0","T0","G0","tauS","tauT","k1","tauG","tauX","type");
      pmodelFit[ifile][ifit]->SetParameters(binwidth,Single0,Triplet0,Gas0,pmodtS,pmodtT,k1,pmodtG,pmodtX,ifit);
      pmodelFit[ifile][ifit]->FixParameter(9,ifit);
      pmodelFit[ifile][ifit]->SetLineColor(pcolor[ifit]);
      pmodelFit[ifile][ifit]->SetMarkerStyle(pstyle[ifit]);
      pmodelFit[ifile][ifit]->Print();
      fmodelFit[ifile][ifit] = new TF1(Form("fmodelFit-%.f-%i",runPPM[ifile],ifit),fmodel,xstart,xstop,8);
      fmodelFit[ifile][ifit]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
      fmodelFit[ifile][ifit]->SetParameters(binwidth,integral, pmodtT, k1 ,pmodtX, Amystery ,pmodtG, ifit );
      fmodelFit[ifile][ifit]->SetNpx(1000); // numb points for function
      fmodelFit[ifile][ifit]->FixParameter(0,binwidth);
      fmodelFit[ifile][ifit]->FixParameter(7,ifit);
      fmodelFit[ifile][ifit]->SetLineColor(pcolor[ifit]);
      fmodelFit[ifile][ifit]->SetLineStyle(10);
      fmodelFit[ifile][ifit]->SetMarkerStyle(pstyle[ifit]);
      fmodelFit[ifile][ifit]->Print();
    }

    canFmodel[ifile] = new TCanvas(Form("fmodel-%.f-PPM",runPPM[ifile]),Form("fmodel-%.f-PPM",runPPM[ifile]));
    canFmodel[ifile]->SetLogy();
    fmodelFit[ifile][4]->Draw("");
    for(int ifit=0; ifit<PFITS; ++ifit) {
      fmodelFit[ifile][ifit]->Draw("same");
      pmodelFit[ifile][ifit]->Draw("same");
    }
  }

  int icolor[MAXFILE]={1,2,3,6,7,8};
  int imark[MAXFILE]={20,21,22,23,24,25};
  TCanvas* canOne[MAXFILE];
  TCanvas* canSumOne[MAXFILE];
  TCanvas* canSumNeil[MAXFILE];

  vector<double> lanFitVal;
  vector<double> lanFitWid;
  vector<double> lanFitValErr;
  vector<double> lanFitWidErr;
  vector<double> lanFitVal2;
  vector<double> lanFitWid2;
  vector<double> lanFitValErr2;
  vector<double> lanFitWidErr2;
  vector<double> vppm;
  vector<double> vppmErr;

  
  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    canOne[ifile] = new TCanvas(Form("input-%.0f-%s",runPPM[ifile],trigTag.Data()),Form("input-%.0f-%s",runPPM[ifile],trigTag.Data()));
    vppm.push_back( runPPM[ifile]);
    vppmErr.push_back(0);
    canOne[ifile]->SetLogy();
    hLife[ifile]->SetMarkerSize(0.5);
    hLife[ifile]->SetMarkerStyle(imark[ifile]);
    hLife[ifile]->SetMarkerColor(icolor[ifile]);
    hLife[ifile]->SetLineColor(icolor[ifile]);
    hLife[ifile]->Draw("E");
    canOne[ifile]->Print(".pdf");

    // landau
    hQSum[ifile]->SetMarkerSize(0.5);
    hQSum[ifile]->SetMarkerStyle(imark[ifile]);
    hQSum[ifile]->SetMarkerColor(icolor[ifile]);
    hQSum[ifile]->SetLineColor(icolor[ifile]);

    canSumOne[ifile] = new TCanvas(Form("EventQSum-%.0f-%s",runPPM[ifile],trigTag.Data()),Form("EventQSum-%.0f-%s",runPPM[ifile],trigTag.Data()));
    hQSum[ifile]->Fit("landau");
    TF1 *lfit =  hQSum[ifile]->GetFunction("landau");
    lanFitVal.push_back(lfit->GetParameter(1));
    lanFitValErr.push_back(lfit->GetParError(1));
    lanFitWid.push_back(lfit->GetParameter(2));
    lanFitWidErr.push_back(lfit->GetParError(2));
    hQSum[ifile]->Draw("E");
    canSumOne[ifile]->Print(".pdf");
  }


  for(int ifile=0; ifile<MAXFILE; ++ifile) {  
    hQSum2[ifile]->SetMarkerSize(0.5);
    hQSum2[ifile]->SetMarkerStyle(imark[ifile]);
    hQSum2[ifile]->SetMarkerColor(icolor[ifile]);
    hQSum2[ifile]->SetLineColor(icolor[ifile]);
    canSumNeil[ifile] = new TCanvas(Form("EventQNeil-%.0f-%s",runPPM[ifile],trigTag.Data()),Form("EventQNeil-%.0f-%s",runPPM[ifile],trigTag.Data()));
    hQSum2[ifile]->Fit("landau");
    TF1 *lfit2 =  hQSum2[ifile]->GetFunction("landau");
    lanFitVal2.push_back(lfit2->GetParameter(1));
    lanFitValErr2.push_back(lfit2->GetParError(1));
    lanFitWid2.push_back(lfit2->GetParameter(2));
    lanFitWidErr2.push_back(lfit2->GetParError(2));
    hQSum2[ifile]->Draw("E");
    canSumNeil[ifile]->Print(".pdf");
  }

  printf(" ........... PPM  ............. \n");
  for(int ifile=0; ifile<MAXFILE; ++ifile)  {
    double r1,re1;
    double r2,re2;
    ratioE(lanFitVal[ifile],lanFitValErr[ifile],lanFitVal[0],lanFitValErr[0],r1,re1);
    ratioE(lanFitVal2[ifile],lanFitValErr2[ifile],lanFitVal2[0],lanFitValErr2[0],r2,re2);

    printf(" PPM %2.0f mean %.2f +/- %.2f  wid %.2f  +/- %.2f mean2 %.2f +/- %.2f  wid %.2f  +/- %.2f  blue %.3f +/- %.3f red %.3f +/- %.3f \n",
        runPPM[ifile],lanFitVal[ifile],lanFitValErr[ifile],lanFitWid[ifile],lanFitWidErr[ifile],
        lanFitVal2[ifile],lanFitValErr2[ifile],lanFitWid2[ifile],lanFitWidErr2[ifile],r1,re1,r2,re2);
  }

  TGraphErrors *glandau = new TGraphErrors(vppm.size(), &vppm[0],&lanFitVal[0],&vppmErr[0],&lanFitValErr[0]);
  TGraphErrors *glandau2 = new TGraphErrors(vppm.size(), &vppm[0],&lanFitVal2[0],&vppmErr[0],&lanFitValErr2[0]);
  TCanvas *clandau = new TCanvas(Form("landau-%s",trigTag.Data()),Form("landau-%s",trigTag.Data()));
  clandau->SetGridx(); clandau->SetGridy();
  glandau->SetTitle(Form("landau-%s",trigTag.Data()));
  glandau->SetMarkerColor(kRed);
  glandau->SetLineColor(kRed);
  glandau->SetMarkerStyle(22);
  glandau->SetMarkerSize(1);
  glandau->GetXaxis()->SetTitle(" Xe ppm ");
  glandau->GetHistogram()->GetYaxis()->SetRangeUser(100,2000);
  glandau2->SetTitle("landauNeil");
  glandau2->SetMarkerColor(kBlue);
  glandau2->SetLineColor(kBlue);
  glandau2->SetMarkerStyle(21);
  glandau2->SetMarkerSize(1);
  glandau2->GetXaxis()->SetTitle(" Xe ppm ");
  glandau->Draw("apl");
  glandau2->Draw("plsame");
  clandau->Print(".png");


  TCanvas* canAll = new TCanvas(Form("InputAll-%s",trigTag.Data()),Form("InputAll-%s",trigTag.Data()));
  canAll->SetLogy();
  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    if(ifile==0) hLife[ifile]->Draw("E");
    else hLife[ifile]->Draw("Esame");
  }
  canAll->Print(".png");

  
  TCanvas* canSumAll = new TCanvas(Form("EventQSum-%s",trigTag.Data()),Form("EventQSum-%s",trigTag.Data()));
  canSumAll->SetLogy();
  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    if(ifile==0) hQSum[ifile]->Draw("E");
    else hQSum[ifile]->Draw("Esame");
  }
  canSumAll->Print(".png");



  ///fin->GetObject("LifePass",hLife[1]);
  //fin->GetObject("LifeFail",hLife[2]);

  for(int ih=0; ih<MAXFILE; ++ih) {
    hLife[ih]->ls(); fout->Append(hLife[ih]); fout->Append(hQSum[ih]);
  }
 
  TCanvas *can[MAXFILE];
  //gStyle->SetOptFit();


  TH1D* hLifeResiduals = new TH1D("LifeAllResiduals"," residuals ",2000,-200,200);
  TNtuple *nResidual = new TNtuple("nResidual","","bin:x:res:fit:val:error");

  
  double maxLife=10.0;
  double qCutLow=5;
  double qCut=100;

  
    /*
  double A0 =  fmodelFit[4][4]->GetParameter(2);
  double ta = fmodelFit[4][4]->GetParameter(3);
  double kone = fmodelFit[4][4]->GetParameter(4);
  double tx = fmodelFit[4][4]->GetParameter(5);
  double tap = 1./( 1./ta + kone);
  double t1inv = ( -kone + 1./tx);
  double t2inv = ( -1./tap + 1./tx);


  printf(" tap %f tapinv %f tdpinv %f txinv %f 1/t1=taudpinv-tauxinv %f 1/t2=taudpinv-tauxinv %f \n",tap,1./tap,kone,1./tx,t1inv,t2inv);
  */

  /*
  for(unsigned j=0; j<100; ++j ) {
    double x = xstart+double(j)/10.;
    double xx = double(j)/10.;
    double DD =  A0*ta/tdp*(Exp(-xx/tdp) - Exp(-xx/tap)) ;
    printf(" x=%f A %.3E D %.3E (%.3E) X %.3E T %.3E \n",x, modelFitA->Eval(x),modelFitD->Eval(x),DD,modelFitX->Eval(x),modelFit->Eval(x));
    }
    */

  std::vector<double> vprompt;
  std::vector<double> vtot;
  double lowVal=0.9;
  double highVal=1.10;
  Int_t plowBin = hLife[0]->FindBin(lowVal);
  Int_t phighBin = hLife[0]->FindBin(highVal);
  double afterLow=1.45;
  double afterHigh=1.6;
  Int_t afterLowBin = hLife[0]->FindBin(afterLow);
  Int_t afterHighBin = hLife[0]->FindBin(afterHigh);


  for(int ih=0; ih<MAXFILE; ++ih) {
    if(hLife[ih]==NULL) continue;
    hLife[ih]->SetStats(1);
    //if(ih>0) continue;
    Double_t tguess  = (hLife[ih]->GetBinWidth(1))*Double_t(hLife[ih]->GetNbinsX())/10;
    Int_t lowBin = hLife[ih]->FindBin(xstart);
    Int_t highBin = hLife[ih]->FindBin(xstop);

    integral = hLife[ih]->Integral();
    double prompt = 0; //hLife[ih]->Integral(plowBin,phighBin);
    double tot =0;

    for(int jbin=plowBin; jbin< hLife[ih]->GetNbinsX(); ++jbin) {
      double jval = hLife[ih]->GetBinContent(jbin);
      if(jbin<phighBin) prompt += jval;
      if(jbin<afterLowBin||jbin>afterHighBin) tot += jval;
    }
    vprompt.push_back(prompt);
    vtot.push_back(tot);

    double kval = kvalFit(runPPM[ih]);
    printf(" ffffffff  fit %i xppm %f k1 %f  range xstart %f bin %i xstop %f bin %i %f \n",ih,runPPM[ih],kval,xstart,lowBin,xstop,highBin,integral);
    bool PFIT=true;
    if(!PFIT) {
      //TF1 fit =createModelFit(ih,binwidth,runPPM[ih],integral,tauA,tauD,tauX,Amystery,tauM);
      int ifit = ih+100;
      fp[ih] = new TF1(Form("fmodel%i",ifit),fmodel,xstart,xstop,8);
      fp[ih]->SetNpx(1000); // numb points for function
      fp[ih]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
      fp[ih]->SetParameters(binwidth,integral, tauA, kval , tauX, Amystery , tauM, ifit );
      fp[ih]->FixParameter(0,binwidth);
      fp[ih]->FixParameter(7,ifit);
      fp[ih]->FixParameter(2,tauA);
      //fp[ih]->FixParameter(3,kval);
      fp[ih]->FixParameter(4,tauX);
      fp[ih]->FixParameter(6,tauM);
      //fp[ih]->FixParameter(7,tauM);
      fp[ih]->SetParLimits(5,0,0.1*integral);

    } else {
      int ifit = ih+100;
      fp[ih] = new TF1(Form("pmodelFit-%i",ifit),pmodel,xstart,xstop,10);
      fp[ih]->SetNpx(1000); // numb points for function
      fp[ih]->SetParNames("binwidth","S0","T0","G0","tauS","tauT","k1","tauG","tauX","type");
      fp[ih]->SetParameters(binwidth,Single0,Triplet0,Gas0,pmodtS,pmodtT,kval,pmodtG,pmodtX,ifit);
      fp[ih]->FixParameter(0,binwidth);
      fp[ih]->FixParameter(4,pmodtS);
      fp[ih]->FixParameter(5,pmodtT);
      //fp[ih]->FixParameter(6,kval);
      fp[ih]->FixParameter(7,pmodtG);
      fp[ih]->FixParameter(8,pmodtX);
      fp[ih]->FixParameter(9,ifit);
      fp[ih]->SetLineColor(pcolor[ih]);
      fp[ih]->SetMarkerStyle(pstyle[ih]);
      fp[ih]->SetParLimits(1,0,2*integral);
      fp[ih]->SetParLimits(2,0,2*integral);
      fp[ih]->SetParLimits(3,0,0.1*integral);
      fp[ih]->SetParLimits(6,0,10);
    }

    printf(" \t number %i  %s nbins %i lowbin %i xstart %E xstop %E integral %E width %E \n",
        ih,hLife[ih]->GetName(),hLife[ih]->GetNbinsX(),lowBin,xstart,xstop,integral,binwidth); //,hsum1->GetNa
    hLife[ih]->Fit(fp[ih],"RLE0+");
    printf(" fit %i tauL = %f +- %f micro-sec \n",ih,fp[ih]->GetParameter(3),fp[ih]->GetParError(3));
    fp[ih]->Print();
    //for(int ii=0; ii<8 ; ++ii) printf("\t  param %i %s %.3f \n",ii,fp[ih]->GetParName(ii),fp[ih]->GetParameter(ii));
    int ii;
    printf(" xxxx  %.f PPM \n",runPPM[ih]);
    ii=1; printf("\t xxxx  param %i %s %.3f +/- %.3f \n",ii,fp[ih]->GetParName(ii),fp[ih]->GetParameter(ii)/integral,
        fp[ih]->GetParError(ii)/integral);
    ii=2; printf("\t xxxx  param %i %s %.3f +/- %.3f \n",ii,fp[ih]->GetParName(ii),fp[ih]->GetParameter(ii)/integral,
        fp[ih]->GetParError(ii)/integral);
    ii=3; printf("\t xxxx  param %i %s %.3f +/- %.3f \n",ii,fp[ih]->GetParName(ii),fp[ih]->GetParameter(ii)/integral,
        fp[ih]->GetParError(ii)/integral);
    ii=6; printf("\t xxxx  param %i %s %.3f +/- %.3f \n",ii,fp[ih]->GetParName(ii),fp[ih]->GetParameter(ii),
        fp[ih]->GetParError(ii));
  }

  /*
  fp[0]->Print();
  for(int ii=0; ii<8 ; ++ii) {
    printf(" param %i %s %.3f \n",ii,fp[0]->GetParName(ii),modelFit->GetParameter(ii));
  }
  */

  double promptNorm=vprompt[0];
  double totNorm=vtot[0];

  for(unsigned ii=0 ; ii<vppm.size() ; ++ii) {
    double vt = vtot[ii];
    vprompt[ii]=vprompt[ii]/promptNorm;
    vtot[ii]=vtot[ii]/totNorm;
    printf(" %i %f %f %F \n",ii,vprompt[ii],vtot[ii],vt);
  }
 
  TGraph *gprmpt = new TGraph(vppm.size(),&vppm[0],&vprompt[0]);
  TCanvas *cprompt = new TCanvas(Form("promptQ-%.1f-%.1f",lowVal,highVal),Form("promptQ-%.1f-%.1f",lowVal,highVal));
  cprompt->SetGridx(); cprompt->SetGridy();
  gprmpt->SetTitle(Form("promptQ-%.1f-%.1f",lowVal,highVal));
  gprmpt->SetMarkerColor(kRed);
  gprmpt->SetMarkerStyle(22);
  gprmpt->SetMarkerSize(1);
  gprmpt->GetXaxis()->SetTitle(" Xe ppm ");
  gprmpt->Draw("ap");
  cprompt->Print(".png");


   
  TGraph *gtot = new TGraph(vppm.size(),&vppm[0],&vtot[0]);
  TCanvas *ctot = new TCanvas("totQ","totQ");
  ctot->SetGridx(); ctot->SetGridy();
  gtot->SetTitle("totQ");
  gtot->SetMarkerColor(kBlue);
  gtot->SetMarkerStyle(22);
  gtot->SetMarkerSize(1);
  gtot->GetXaxis()->SetTitle(" Xe ppm ");
  gtot->Draw("ap");
  ctot->Print(".png");




  for(int ibin =1; ibin< hLife[0]->GetNbinsX() -1 ; ++ibin) {
    double x =  hLife[0]->GetBinLowEdge(ibin);
    if(x<xstart) continue;
    double fitval = fp[0]->Eval(x);
    double residual = 0;
    if(hLife[0]->GetBinError(ibin)>0) residual=(hLife[0]->GetBinContent(ibin)-fitval)/hLife[0]->GetBinError(ibin);
    //if(ibin<130) printf(" %i %f %f %f %f \n",ibin,x,fitval,hLife[0]->GetBinContent(ibin),residual);
    nResidual->Fill(float(ibin),x,residual,fitval,hLife[0]->GetBinContent(ibin), hLife[0]->GetBinError(ibin));
    hLifeResiduals->Fill(residual);
  }

  gStyle->SetOptFit();

  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    can[ifile] = new TCanvas(Form("life-%s-pass-%i-hits",runTag[ifile].Data(),10),Form("life-%s-passs-%i-hits",runTag[ifile].Data(),10));
    hLife[ifile]->SetTitle( Form("life-%s-pass-%i-hits",runTag[ifile].Data(),10) );
    can[ifile]->SetLogy();
    //hLife[ifile]->SetStats(1);
    hLife[ifile]->Draw("E");
    if(fp[ifile]) fp[ifile]->Draw("same");
    can[ifile]->Print(".png");
    if(fp[ifile]) fout->Append(fp[ifile]);
  }


  
  fout->ls();
  fout->Write();

}
