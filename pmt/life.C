using namespace TMath;
enum {NPMT=2};
static  double qmax = 0.5;
static double qmin = 0.0;
static Int_t nbins = 500;



static double fexp(double *xx, double *par)
{
  double binwidth = 4.0E-9;
  double  x= xx[0];
  double tau = par[0];
  double f = par[1]/tau*binwidth*Exp(-x/tau);
  return f; 
}


//void life(TString tag = "250kohms_0_5Sigma")
void life(TString tag = "LAr_1_0_5Sigma") 
{
  TString fileName ; fileName.Form("baconRunAna_%s.root",tag.Data());

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(!_file0) {
    cout << " not finding file " << fileName << endl;
    return;
  }

  TNtuple *ntEvent;
  _file0->GetObject("ntEvent",ntEvent);

  TH1D *hfsum0;
  _file0->GetObject("SumPmt0",hfsum0);

  TH1D *hfsum1;
  //_file0->GetObject("SumPmt1",hfsum1);


  TH1D* hsum0 = (TH1D*) hfsum0->Clone("sum0-sub");
  //TH1D* hsum1 = (TH1D*) hfsum1->Clone("sum1-sub");


  Double_t base[NPMT];
  for(int i=0; i<NPMT; ++i ) base[i]=0;

  Int_t istart =  hsum0->FindBin(0);
  Int_t iend =  hsum0->FindBin(.2E-6);
  for(int ibin = istart; ibin < iend; ++ibin) base[0] += hfsum0->GetBinContent(ibin);
  //for(int ibin = istart; ibin < iend; ++ibin) base[1] += hfsum1->GetBinContent(ibin);
  
  base[0] /= Double_t(  iend - istart);
  base[1] /= Double_t(  iend - istart);


  printf("  base 0 %f base 1 %f \n",base[0],base[1]);


  for(int ibin = 0; ibin < hsum0->GetNbinsX(); ++ibin) hsum0->SetBinContent( ibin,  hfsum0->GetBinContent(ibin) - base[0]);
  //for(int ibin = 0; ibin < hsum1->GetNbinsX(); ++ibin) hsum1->SetBinContent( ibin,  hfsum1->GetBinContent(ibin) - base[1]);

  printf(" total number of events is %lli \n",ntEvent->GetEntries());
  if(hsum0==NULL ) return;

  printf(" %s %s \n",hsum0->GetName()); //,hsum1->GetName());


  // fit
  Double_t xlow=1.00E-6;
  Double_t xhigh=2.0E-6;

  TF1 *fp[NPMT];  // on for each pmt
  for(int i=0; i<NPMT; ++i ) {
    fp[i] = new TF1(Form("fexp-%i",i),fexp,xlow,xhigh,2);
    fp[i]->SetNpx(1000); // numb points for function
    fp[i]->SetParName(0,"lifetime");
    fp[i]->SetParName(1,"integral");
    fp[i]->SetParameter(0,1E-6);
  }

  gStyle->SetOptFit();
  TString all0Name;
  all0Name.Form("qsumPmt0-%s",tag.Data());
  TCanvas *call0= new TCanvas(all0Name,all0Name);
  //call0->SetLogy();
  hsum0->Fit(fp[0],"R+");
  hsum0->Draw();
  call0->Print(".pdf");
  printf(" pmt %i tau = %f +- %f micro-sec \n",0,fp[0]->GetParameter(0)*1E6,fp[0]->GetParError(0)*1E6);

  return;

  TString all1Name;
  all1Name.Form("qsumPmt1-%s",tag.Data());
  TCanvas *call1= new TCanvas(all1Name,all1Name);
  //call1->SetLogy();
  //hsum1->Fit(fp[1],"R+");
  //hsum1->Draw();
  //call1->Print(".pdf");
  //printf(" pmt %i tau = %f +- %f micro-sec \n",1,fp[1]->GetParameter(0)*1E6,fp[1]->GetParError(0)*1E6);


  return;

}
