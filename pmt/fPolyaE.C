// test polya distributions
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TKey.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"

using namespace TMath;
TString htag;
std::vector<TH1D*> hlist;
static double nsum;

static double background(double *x, double *par)
{
   return par[0]*exp(-x[0]/par[1]);
}
// polya + exp 
static double fpolyaE(double *x, double *par)
{
  double e = par[0];
  double g = Gamma(par[0]);
  double qn = par[1];
  double y = e*x[0]/qn;
  double b = par[2];
  double xmin = par[5];
  double xmax = par[6];
  double norm = Gamma(e,xmax/qn*e)-Gamma(e,xmin/qn*e);
  double norme = Exp(-xmin/par[2]) -Exp(-xmax/par[2]);
  double C=Exp(-y);
  double sig  = par[3]*e*pow(y,e-1)*C/g/par[1]/norm;
  double bg = (1-par[3])*exp(-x[0]/par[2])/par[2]/norme;
  //printf(" %f %f \n ",sig,bg);
  double f= par[4]*(sig+bg); 
  return TMath::Max(0.0,f);
}


/*****************************************************
   read all histograms 
*******************************************************/

void reading(TDirectory *fdir) 
{
  TObject *obj=NULL;
  TKey* key;
  TIter nextkey(fdir->GetListOfKeys());
  htag = TString("Raw");
  //htag = TString("QhitOff");
  while ( (key = (TKey*)nextkey()) )  {
    fdir->cd();
    obj = key->ReadObj();
    if(obj->IsA()->InheritsFrom("TH1D")) { //case of TH1 or TProfile
      TH1D* h1 = dynamic_cast<TH1D*>(obj);
      if( TString(h1->GetName()).Contains(htag) ) hlist.push_back(h1);
    } 
  }
}

void fPolyaE(TString tag="07-31-1555_0")
{
  enum {NPMT=21};
  TString inputFileName = TString("pmtAna_")+tag+TString(".root");
  printf(" opening file %s \n",inputFileName.Data()); 
  TFile *infile = new TFile(inputFileName);
  //infile->ls();
  reading(infile);
  for(UInt_t ih=0; ih<hlist.size(); ++ih) {
    printf(" %i %s \n",ih,hlist[ih]->GetName());
    hlist[ih]->GetListOfFunctions()->Clear();
    hlist[ih]->Sumw2();
  }


  double en=2;
  double  qn=40;
  int nbins = hlist[0]->GetNbinsX();
  double xmin = hlist[0]->GetBinLowEdge(3);
  double xmax = hlist[0]->GetBinLowEdge(nbins-1);
  double binWidth = hlist[0]->GetBinWidth(0); // all same
  double enorm = Exp(-xmin/qn*en)- Exp(-xmax/qn*en);

  printf(" fitting from %f to %f \n",xmin,xmax);

  TF1 *fp[NPMT+1];
  for(int i=0; i<NPMT+1; ++i ) {
    if(i==7||i==16||i==5||i==17) xmin = hlist[0]->GetBinLowEdge(5);
    else xmin = hlist[0]->GetBinLowEdge(3);
    fp[i] = new TF1(Form("poyaE-%i",i),fpolyaE,xmin,xmax,7);
    fp[i]->SetNpx(1000); // numb points for function
    fp[i]->SetParName(0,"primary charge en");
    fp[i]->SetParName(1,"gain Qn");
    fp[i]->SetParName(2,"norm");
    fp[i]->SetParName(2,"exp const");
    fp[i]->SetParName(3,"sig frac");
    fp[i]->SetParName(4,"norm");
    fp[i]->SetParName(5,"xmin");
    fp[i]->SetParName(6,"xmax");
    fp[i]->SetParameter(0,en);
    fp[i]->SetParameter(1,qn);
    fp[i]->SetParameter(2,1);
    fp[i]->SetParLimits(2,0.5,10); // limit range of exp factor
    fp[i]->SetParameter(3,0.7);
    fp[i]->SetParLimits(3,0,1); // limit range signal 
    fp[i]->SetParameter(4,1.0);
    fp[i]->FixParameter(5,xmin);
    fp[i]->FixParameter(6,xmax);
    fp[i]->SetLineColor(kBlue);
  }
  
  
  gStyle->SetOptStat();
  gStyle->SetOptFit();

  TCanvas *cpolyaE = new TCanvas("polyaE","polyaE");
  fp[NPMT]->SetTitle("poyaE-start");
  fp[NPMT]->Draw();
  fp[NPMT]->Print("V");

 
  TString canTitle;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  TCanvas* can[NPMT];
  for(UInt_t ih=0; ih<NPMT; ++ih) { 
    printf(" \n *********** fitting %s to pmt %i ************ \n",htag.Data(),ih);
    canTitle.Form("%s-pmt%i",htag.Data(),ih);
    can[ih] = new TCanvas(canTitle,canTitle);
    gPad->SetLogy();
    hlist[ih]->Fit(fp[ih],"R+");
    hlist[ih]->SetLineColor(kBlack);
    hlist[ih]->Draw();
    can[ih]->Update();
    can[ih]->Print(".pdf");
  }

  cout << " fits " << htag << endl;
  for(int ih=0; ih<NPMT; ++ih ) {
    double width = hlist[ih]->GetBinWidth(0); // all same
    double qmax = hlist[ih]->GetBinLowEdge(hlist[ih]->GetMaximumBin()) + 0.5*width;
    double qfit=fp[ih]->GetParameter(1); 
    double qfite=fp[ih]->GetParError(1); 
    printf(" %i max bin %.3f qfit %.3f +/- %.3f \n",ih,qmax,qfit,qfite);
  }

  cout << " // fits " << htag << endl;
  if(htag.Contains("NoBeam")) {
      for(int ih=0; ih<NPMT; ++ih ) printf(" fitNoBeam[%i]=  %.3f ;\n",ih,fp[ih]->GetParameter(1));
  } else {
      for(int ih=0; ih<NPMT; ++ih ) printf(" fitRaw[%i]=  %.3f ; fitRawE[%i]= %.3f ;\n",ih,fp[ih]->GetParameter(1), ih, fp[ih]->GetParError(1));
  }
      
}

