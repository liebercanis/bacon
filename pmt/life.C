using namespace TMath;
enum {NPMT=1};
static  double qmax = 0.5;
static double qmin = 0.0;
static Int_t nbins = 500;



static double fexp(double *xx, double *par)
{
  double binwidth = par[2];
  double  x= xx[0];
  double tau = par[0];
  double f = par[1]/tau*binwidth*Exp(-x/tau);
  return f; 
}


//void life(TString tag = "250kohms_0_5Sigma")
void life(TString tag = "simEvents_20181218_10000_Ev_100_derivative") 
{
  TString fileName ; fileName.Form("%s.root",tag.Data());

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(!_file0) {
    cout << " not finding file " << fileName << endl;
    return;
  }

  TH1D* hlife0;
  _file0->GetObject("Life0",hlife0);

  hlife0->Rebin(10);
  // fit
   // fit
  Double_t binwidth = hlife0->GetBinWidth(1);
  Int_t lowBin = hlife0->FindBin(0);
  Double_t xlow=hlife0->GetBinLowEdge(lowBin);
  Double_t xhigh=hlife0->GetBinLowEdge(hlife0->GetNbinsX());
  Double_t integral = hlife0->Integral();
  printf(" %s nbins %i lowbin %i xlow %E xhigh %E integral %E width %E \n",hlife0->GetName(),hlife0->GetNbinsX(),lowBin,xlow,xhigh,integral,binwidth); //,hsum1->GetNa
 
  TF1 *fp[NPMT];  // on for each pmt
  for(int i=0; i<NPMT; ++i ) {
    fp[i] = new TF1(Form("fexp-%i",i),fexp,xlow,xhigh,3);
    fp[i]->SetNpx(1000); // numb points for function
    fp[i]->SetParName(0,"lifetime");
    fp[i]->SetParName(1,"integral");
    fp[i]->SetParName(2,"binwidth");
    fp[i]->SetParameter(0,1E-6);
    fp[i]->SetParameter(1,integral);
    fp[i]->FixParameter(2,binwidth);
  }

  gStyle->SetOptFit();
  TString all0Name;
  all0Name.Form("LifePmt0-%s",tag.Data());
  TCanvas *call0= new TCanvas(all0Name,all0Name);
  call0->SetLogy();
  hlife0->Fit(fp[0],"RL");
  hlife0->Draw();
  call0->Print(".pdf");
  printf(" pmt %i tau = %f +- %f micro-sec \n",0,fp[0]->GetParameter(0)*1E6,fp[0]->GetParError(0)*1E6);

  return;

}
