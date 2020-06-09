
using namespace TMath;
static double xstart=1;
static double xstop=10;
static double tiny =.000021; 
/*I fix all of the lifetimes: microseconds*/
static double pmodtS = 7E-3;
static double pmodtT = 1.96;
static double pmodtG = 3.2;
static double pmodtX = 20E-3;


static double kvalFit(double x) {
  double k0 = 3.53;  // offset for undoped argon
  return 0.196*(x+k0)+0.689;
}

static void ratioE(double x, double xe, double y, double ye, double& r,double& re){
  r = x/y;
  re = r * sqrt( pow(xe/x,2.)+pow(ye/y,2.) );
}


// 9 parameter model !
static double pmodel(double *xx, double *par) 
{
  int type = int(par[9]);
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double S0 = par[1];
  double T0 = par[2];
  double G0 = par[3];
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
  T=T/tt;
  S=S/ts;
  X=X/tx;
  G=G/tg;
  if(type ==0) f=T+S;
  else if(type ==1) f=M;
  else if(type ==2) f=X;
  else if(type ==3) f=G;
  else f =  S+T+G+X;  
  return f;
}


// fmodelFit[ifit]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
static double fmodel(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  int type = int(par[7]);
  if(t<0) return tiny;
  double binwidth = par[0];
  double A0 = par[1];
  double ta = par[2];
  double k = par[3];
  double tx = par[4];
  double Am = par[5];
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
  A=A/ta;
  X=X/tx;
  M=M/tm;

  if(type ==0) f=A;
  else if(type ==1) f=D;
  else if(type ==2) f=X;
  else if(type ==3) f=M;
  else f =  A+X+M; 
  return f;
}




void fits() {
  enum {PFITS=5};
  enum {NHIST=5};

  TFile *cfile = new TFile("fits.root","RECREATE");


  TString runTag[PFITS];
  runTag[0]=TString("00_PPM");
  runTag[1]=TString("01_PPM");
  runTag[2]=TString("02_PPM");
  runTag[3]=TString("05_PPM");
  runTag[4]=TString("10_PPM");

  
  double runPPM[PFITS]={0,1,2,5,10};
  TFile *file1 = new TFile("comps.root","readonly");

  //     Life-XeD_10ppm_N2_10ppm_Doping;
  TH1D* hNLife[PFITS];
  TString Nnames[PFITS];
  TString newName;
  Nnames[0]=TString("NLifePureNorm");
  Nnames[1]=TString("NLife1ppmNorm");
  Nnames[2]=TString("NLife2ppmNorm");
  Nnames[3]=TString("NLife5ppmNorm");
  Nnames[4]=TString("NLife10ppmNorm");

  TH1D* hist;
  for(int ifit =0; ifit<PFITS; ++ifit )  {
    file1->GetObject(Nnames[ifit].Data(),hist);
    hNLife[ifit] = (TH1D*) hist->Clone(Form("NLife-%s",runTag[ifit].Data()));
    double total =  hNLife[ifit]->Integral();
    printf(" file %i total %f integral %f bw %f \n",ifit,total,hNLife[ifit]->Integral(),hNLife[ifit]->GetBinWidth(1));
    cfile->Append(hNLife[ifit]);
  }


  Double_t binwidth = hNLife[0]->GetBinWidth(1);
  Double_t integral = hNLife[0]->Integral();

  // pmodel 
  TF1* pmodelFit[NHIST][PFITS];
  int pcolor[PFITS] = {2,kGreen+2,kBlue,kYellow-7,1};
  int pstyle[PFITS] = {21,22,23,24,25};
  double Single0=0*integral;
  double Triplet0=1*integral;
  double Gas0=Triplet0/10.;

  // fmodel 
  double tauX=pmodtX;
  double tauA = pmodtT;
  double tauM = pmodtG; 
  double Amystery = Gas0;
  TF1* fmodelFit[NHIST][PFITS];

  TCanvas *canFmodel[PFITS];


  for(int ihist=0; ihist<NHIST; ++ihist) {
    double k1= kvalFit(runPPM[ihist]);
    if(ihist==0) k1=0;
    for(int ifit=0; ifit<PFITS; ++ifit) {
      pmodelFit[ihist][ifit] = new TF1(Form("pmodelFit-%.f-%i",runPPM[ihist],ifit),pmodel,xstart,xstop,10);
      pmodelFit[ihist][ifit]->SetNpx(1000); // numb points for function
      pmodelFit[ihist][ifit]->SetParNames("binwidth","S0","T0","G0","tauS","tauT","k1","tauG","tauX","type");
      pmodelFit[ihist][ifit]->SetParameters(binwidth,Single0,Triplet0,Gas0,pmodtS,pmodtT,k1,pmodtG,pmodtX,ifit);
      pmodelFit[ihist][ifit]->FixParameter(9,ifit);
      pmodelFit[ihist][ifit]->SetLineColor(pcolor[ifit]);
      pmodelFit[ihist][ifit]->SetMarkerStyle(pstyle[ifit]);
      pmodelFit[ihist][ifit]->Print();
      fmodelFit[ihist][ifit] = new TF1(Form("fmodelFit-%.f-%i",runPPM[ihist],ifit),fmodel,xstart,xstop,8);
      fmodelFit[ihist][ifit]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
      fmodelFit[ihist][ifit]->SetParameters(binwidth,Triplet0, pmodtT, k1 ,pmodtX, Amystery ,pmodtG, ifit );
      fmodelFit[ihist][ifit]->SetNpx(1000); // numb points for function
      fmodelFit[ihist][ifit]->FixParameter(0,binwidth);
      fmodelFit[ihist][ifit]->FixParameter(7,ifit);
      fmodelFit[ihist][ifit]->SetLineColor(pcolor[ifit]);
      fmodelFit[ihist][ifit]->SetLineStyle(10);
      fmodelFit[ihist][ifit]->SetMarkerStyle(pstyle[ifit]);
      fmodelFit[ihist][ifit]->Print();
    }

    canFmodel[ihist] = new TCanvas(Form("fmodel-%.f-PPM",runPPM[ihist]),Form("fmodel-%.f-PPM",runPPM[ihist]));
    canFmodel[ihist]->SetLogy();
    pmodelFit[ihist][4]->Draw("");
    for(int ifit=0; ifit<PFITS; ++ifit) {
      fmodelFit[ihist][ifit]->Draw("same");
      pmodelFit[ihist][ifit]->Draw("same");
    }

  }
   
  return;
  
  TF1* fp[NHIST];
  bool PMODEL=false;
  for(int ihist=0; ihist<NHIST; ++ihist) {
    double k1= kvalFit(runPPM[ihist]);
    if(PMODEL) {
      fp[ihist]= new TF1(Form("Lfit-%.f-%i",runPPM[ihist],ihist),pmodel,xstart,xstop,10);
      fp[ihist]->SetNpx(1000); // numb points for function
      fp[ihist]->SetParNames("binwidth","S0","T0","G0","tauS","tauT","k1","tauG","tauX","type");
      fp[ihist]->SetParameters(binwidth,Single0,Triplet0,Gas0,pmodtS,pmodtT,k1,pmodtG,pmodtX,4);
      fp[ihist]->FixParameter(0,binwidth);
      fp[ihist]->FixParameter(1,0);
      fp[ihist]->SetParLimits(2,0,2*integral);
      fp[ihist]->FixParameter(4,pmodtS);
      fp[ihist]->FixParameter(5,pmodtT);
      //fp[ihist]->FixParameter(7,pmodtG);
      fp[ihist]->FixParameter(8,pmodtX);
      fp[ihist]->FixParameter(9,4);
    } else {
      fp[ihist] = new TF1(Form("LfFit-%.f-%i",runPPM[ihist],ihist),fmodel,xstart,xstop,8);
      fp[ihist]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
      fp[ihist]->SetParameters(binwidth,Triplet0, pmodtT, k1 ,pmodtX, Amystery ,pmodtG, 4 );
      fp[ihist]->SetNpx(1000); // numb points for function
      fp[ihist]->FixParameter(0,binwidth);
      fp[ihist]->FixParameter(2,pmodtT);
      fp[ihist]->FixParameter(4,pmodtX);
      fp[ihist]->FixParameter(4,pmodtX);
      fp[ihist]->FixParameter(5,3.733E4);
      fp[ihist]->FixParameter(7,4);
    }

    fp[ihist]->SetLineColor(pcolor[ihist]);
    fp[ihist]->SetMarkerStyle(pstyle[ihist]);

    printf(" \t number %i  %s nbins %i  xstart %E xstop %E integral %E width %E \n",
        ihist,hNLife[ihist]->GetName(),hNLife[ihist]->GetNbinsX(),xstart,xstop,integral,binwidth); //,hsum1->GetNa
    hNLife[ihist]->Fit(fp[ihist],"RLE0+");
    fp[ihist]->Print();
    for(int ii=0; ii<10; ++ii) {
      printf(" xxxx  %.f PPM \n",runPPM[ihist]);
      printf("\t xxxx  param %i %s %.3f +/- %.3f \n",ii,fp[ihist]->GetParName(ii),fp[ihist]->GetParameter(ii)/integral,
          fp[ihist]->GetParError(ii)/integral);
    }

  }

  TCanvas *canFit[NHIST];
  gStyle->SetOptStat();
  gStyle->SetOptFit();
  for(int ihist=0; ihist<NHIST; ++ihist) {
    canFmodel[ihist] = new TCanvas(Form("LifeFit-%.f-PPM",runPPM[ihist]),Form("LifeFit-%.f-PPM",runPPM[ihist]));
    canFmodel[ihist]->SetLogy();
    hNLife[ihist]->Draw("p");
    fp[ihist]->Draw("sames");
  }



  cfile->Write();
  cfile->ls();


}
