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
// single polyastatic 
double xmax=10;
static double xmin=0;

static double fpolya(double *x, double *par)
{
  double y = x[0]/par[1];
  double e = par[0];
  double norm = Gamma(e,xmax/par[1]*e)-Gamma(e,xmin/par[1]*e);   
  double f = par[2]*e*pow(e*y,e-1)*Exp(-e*y)/Gamma(e)/par[1]/norm;
  return f; 
}

/*****************************************************
   read all histograms 
*******************************************************/

void reading(TDirectory *fdir) 
{
  TObject *obj=NULL;
  TKey* key;
  TIter nextkey(fdir->GetListOfKeys());
  htag = TString("QhitNoBeam");
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

void fitP(TString tag="simEvents_20190218_10000_Ev_0_derivative")
{
  enum {NPMT=2};
  TString fileName ; fileName.Form("%s.root",tag.Data());

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(!_file0) {
    cout << " not finding file " << fileName << endl;
    return;
  }

   bool isSimulation=false;
   if(fileName.Contains("sim")) isSimulation=true;

   if(isSimulation) printf(" \t\t this is simulation \n");
   else printf(" \t\t this is real data \n");

  TNtuple* ntHit=NULL;
  _file0->GetObject("ntHit",ntHit);
  if(ntHit==NULL) {
    printf(" did not get ntuple ntHit\n");
    return;
  }

  double qnominal = 0.05;  //mV/(single electron)

  TH1D* hinput = new TH1D("charge","charge",400,-2,8);
  ntHit->Draw("q/0.05>>charge","time>2");
  TH1D* hcharge2 = new TH1D("charge2","charge",400,-2,8);
  TH1D* hcharge0 = new TH1D("charge0","charge",400,-2,8);
  if(isSimulation)  {
    ntHit->Draw("q/0.05>>charge2","time>2&&good==2");
    ntHit->Draw("q/0.05>>charge0","time>2&&good==0");
    hcharge0->SetFillColor(kRed);
    hcharge2->SetFillColor(kGreen);
  }

  TString chName; chName.Form("charge-%s",tag.Data());
  TH1D* hq = (TH1D*) hinput->Clone(chName);

  TCanvas *charge = new TCanvas(chName,chName);
  hq->Draw();
  if(isSimulation)  {
    hcharge0->Draw("sames");
    hcharge2->Draw("sames");
  }

 

  double en=5;
  double  qn=2;

  TF1 *fp[NPMT];
  for(int i=0; i<NPMT; ++i ) {
    fp[i] = new TF1(Form("poya-%i",i),fpolya,xmin,xmax,3);
    fp[i]->SetNpx(1000); // numb points for function
    fp[i]->SetParName(0,"primary charge en");
    fp[i]->SetParName(1,"gain Qn");
    fp[i]->SetParName(2,"norm");
    fp[i]->SetParameter(0,en);
    fp[i]->SetParameter(1,qn);
    fp[i]->SetParameter(2,1);
    fp[i]->SetLineColor(kRed);
  }
  TCanvas *cpolyaE = new TCanvas("polyaE","polyaE");
  fp[0]->Draw();
  fp[0]->Print();
 
  
  gStyle->SetOptStat();
  gStyle->SetOptFit();

  TString chFit; chFit.Form("fit-%s",tag.Data());
  TCanvas *cfit = new TCanvas(chFit,chFit);
  hq->Fit(fp[0],"R+");
  hq->Draw();
  cfit->Print(".pdf");

 
  //TCanvas *cpolya = new TCanvas("polya","polya");
  double epar=fp[0]->GetParameter(0); 
  double qpar=fp[0]->GetParameter(1); 
  double pnorm = Gamma(epar,xmax/qpar*epar)-Gamma(epar,xmin/qpar*epar); // from lower limit of xmin
  double econst=fp[0]->GetParameter(2); 
  double norme = Exp(-xmin/econst) - Exp(-xmax/econst);
  printf(" range (%f,%f) exp norm is %f and polya norm is %f \n",xmin,xmax,norme,pnorm );
  float spole=fp[0]->GetHistogram()->Integral("width");
  printf(" bin width %f pol integrates in the range (%f,%f) to  %e\n",fp[0]->GetHistogram()->GetBinWidth(0),xmin,xmax,spole);
  return;
 
  TString canTitle;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  for(UInt_t ih=0; ih<hlist.size(); ++ih) { 
    double hintegral = hlist[ih]->Integral();
    fp[ih]->SetParameter(4,hintegral);
    hlist[ih]->SetNormFactor(hintegral);
    printf(" \n *********** fitting %s to pmt %i ************ \n",htag.Data(),ih);
    hlist[ih]->Fit(fp[ih],"R+");
    canTitle.Form("%s-pmt%i",htag.Data(),ih);
    TCanvas *c0 = new TCanvas(canTitle,canTitle);
    gPad->SetLogy();
    hlist[ih]->SetLineColor(kBlack);
    hlist[ih]->Draw();
    c0->Update();
    c0->Print(".pdf");
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
      for(int ih=0; ih<NPMT; ++ih ) printf(" fitOff[%i]=  %.3f ;\n",ih,fp[ih]->GetParameter(1));
  }
      
}

