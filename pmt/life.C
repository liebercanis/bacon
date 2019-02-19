using namespace TMath;
enum {NPMT=2};
enum {NHIST=6};
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


void life(TString tag = "baconRun_run_4000_0_Ev_0_derivative")
//void life(TString tag = "simEvents_20190211_10000_Ev_0_derivative")
{
  TString fileName ; fileName.Form("%s.root",tag.Data());

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(!_file0) {
    cout << " not finding file " << fileName << endl;
    return;
  }

  TString htag[NHIST];
  TH1D* hlife[NHIST];

  htag[0]=TString("Life0");
  htag[1]=TString("LifeCut0");
  htag[2]=TString("LifeCount0");
  htag[3]=TString("LifeSim0");
  htag[4]=TString("LifeNoise0");
  htag[5]=TString("LifeTrue0");

 
  Double_t tguess;
  for(int ih=0; ih<NHIST; ++ih) {
    _file0->GetObject(htag[ih],hlife[ih]);
    if(hlife[ih]==NULL) { printf(" no histo %s \n",htag[ih].Data()); continue; } 
    tguess = (hlife[ih]->GetBinWidth(1))*Double_t(hlife[ih]->GetNbinsX())/10;
    printf(" name %s bin width %f nbins %i tguess %f\n",hlife[ih]->GetName(),hlife[ih]->GetBinWidth(1),hlife[ih]->GetNbinsX(),tguess);
  }
  TF1 *fp[NHIST];  
  TCanvas *can[NHIST];
  //TH1D* hnlife = (TH1D*) hlife[3]->Clone("");
  //for(int ib=0; ib<hnlife->GetNbinsX(); ++ib) hnlife->SetBinContent(ib, -1.0*hnlife->GetBinContent(ib));
  //hlife[3] = hnlife;

  for(int ih=0; ih<NHIST; ++ih) {
    if(hlife[ih]==NULL) continue;
    printf(" fit %i \n",ih);
    hlife[ih]->Rebin(2);
    // fit
    Double_t binwidth = hlife[ih]->GetBinWidth(1);
    Int_t lowBin = hlife[ih]->FindBin(1.2);
    //if(ih==3) lowBin = hlife[ih]->FindBin(1);
    Double_t xlow=hlife[ih]->GetBinLowEdge(lowBin);
    Int_t highBin = hlife[ih]->FindBin(7.0);
    Double_t xhigh=hlife[ih]->GetBinLowEdge(highBin);
    Double_t integral = hlife[ih]->Integral();
  
    fp[ih] = new TF1(Form("fexp-%i",ih),fexp,xlow,xhigh,3);
    fp[ih]->SetNpx(1000); // numb points for function
    fp[ih]->SetParName(0,"lifetime");
    fp[ih]->SetParName(1,"integral");
    fp[ih]->SetParName(2,"binwidth");
    fp[ih]->SetParameter(0,tguess);
    fp[ih]->SetParameter(1,integral);
    fp[ih]->FixParameter(2,binwidth);

    gStyle->SetOptFit();
    TString all0Name;
    all0Name.Form("%s-%s",htag[ih].Data(),tag.Data());
    can[ih]=new TCanvas(all0Name,all0Name);
    can[ih]->SetLogy();
    hlife[ih]->Fit(fp[ih],"RL");
    hlife[ih]->Draw();
    can[ih]->Print(".png");
    printf(" \t number %i  %s nbins %i lowbin %i xlow %E xhigh %E integral %E width %E \n",
        ih,hlife[ih]->GetName(),hlife[ih]->GetNbinsX(),lowBin,xlow,xhigh,integral,binwidth); //,hsum1->GetNa
    fp[ih]->Print();

    printf(" pmt %i tau = %f +- %f micro-sec \n",0,fp[0]->GetParameter(0)*1E6,fp[0]->GetParError(0)*1E6);
  }
  return;

}
