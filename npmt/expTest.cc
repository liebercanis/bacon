using namespace TMath;
/** time unit is microseconds **/
static double xstart=6;
static double xstop=10;
static double tiny=1E-9;


/*I fix all of the lifetimes: microseconds*/
static double pmodtS = 7E-3;
static double pmodtT = 1.96;
static double pmodtG = 3.2;
static double pmodtX = 20E-3;


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
  double ta = par[2];
  double k = par[3];
  double tx = par[4];
  double tg = par[6];
  double tap = 1./( 1./ta + k);
  double cnorm = 1;//Exp(-xstart/tap);
  double gnorm = 1;//Exp(-xstart/tg);
  double A0 = par[1]*binwidth/cnorm;
  double A = A0*Exp(-t/tap)/ta;
  double G0 = par[5]*binwidth/gnorm;
  double t1=1./( 1./tx - k);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta*k;
  double D = C*( Exp(-t*k) - Exp(-t/tap) );
  double X = C/tx*k*(t1*Exp(-t*k) - t2*Exp(-t/tap)) + C/tx*k*(t2-t1)*Exp(-t/tx);
  double G = G0/tg*Exp(-t/tg);
  double f;
  if(type ==0) f=A;
  else if(type ==1) f=D;
  else if(type ==2) f=X;
  else if(type ==3) f=G;
  else if(type ==4) f=A+X;  
  else f=A+X+G; 
  return f;
}



static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double tau = par[0];
  double cnorm = Exp(-xstart/tau) - Exp(-xstop/tau);
  double binwidth = par[2];
  double f = binwidth/cnorm*par[1]/tau*Exp(-x/tau);
  return f; 
}



void expTest(int nsbin=8,unsigned MAXGEN=1000000) {


  enum {MAXSETS=6};
  double runPPM[MAXSETS]={0,1,2,5,10,10};
  TString runTag[MAXSETS];
  runTag[0]=TString("00PPM-MU");
  runTag[1]=TString("01PPM-Ran");
  runTag[2]=TString("02PPM-Ran");
  runTag[3]=TString("05PPM-Ran");
  runTag[4]=TString("10PPM-Ran");
  runTag[5]=TString("10PPM-Mu");


  TFile *fout = new TFile("expTestOut.root","RECREATE");

  Double_t binwidth = double(nsbin)/1000.;
  double maxLife=10.;
  int nbins = int(maxLife/binwidth);
  double lifeBins = maxLife/binwidth; // ns bins
  double tripletT = 1.3;
  double integral = double(MAXGEN);


  // fmodel 
  double tauX=pmodtX;
  double tauA = pmodtT;
  double tauM = pmodtG; 

  double gfrac = 0.1;
  double A0=(1.-gfrac)*double(MAXGEN);
  double G0 = A0*gfrac;
 
  int setColor[MAXSETS]={kRed,kGreen,kBlue,kYellow-6,kGreen-6,1};
  int setStyle[MAXSETS]={21,22,23,24,25,26};

  int iset=5;
  

  TF1* fp = new TF1(Form("fexp-%.2f",tripletT),fexp,xstart,xstop,3);
  fp->SetNpx(1000); // numb points for function
  fp->SetParNames("lifetime","integral","binwidth");
  fp->SetParameters(tripletT,integral,binwidth);
  fp->FixParameter(2,binwidth);

  TF1 * ftail = new TF1(Form("fexp-%i",iset),fexp,xstart,xstop,3);
  ftail->SetNpx(1000); // numb points for function
  ftail->SetParNames("lifetime","integral","binwidth");
  ftail->SetParameters(pmodtG ,double(MAXGEN),binwidth);
  //ftail->FixParameter(0,pmodtG);
  ftail->FixParameter(2,binwidth);
  for(int ii=0; ii<3; ++ii) {
    printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
  }
  double sum = ftail->Integral(xstart,xstop)/binwidth; 
  double cnorm = Exp(-xstart/pmodtG) - Exp(-xstop/pmodtG);

  printf( " gen %.2E cnorm %.2E atarting sum from function ftail %.2E \n", double(MAXGEN),cnorm ,sum);


  //hLifeSum[iset]->Fit(ftail,"RLE0+");
  TH1D* hTail = new TH1D("Tail","",nbins,0,maxLife); 
  for(unsigned ig=0; ig<MAXGEN ; ++ig) hTail->Fill(ftail->GetRandom());
  hTail->Fit(ftail,"RLE0+");

  for(int ii=0; ii<3; ++ii) {
    printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
  }
  
  sum = ftail->Integral(xstart,xstop)/binwidth;
  printf( " sum from function ftail %.2E \n", sum);


  TCanvas *canTail;
  gStyle->SetOptStat();
  gStyle->SetOptFit(1111111);
  canTail = new TCanvas(Form("TailFit-%.f-PPM",runPPM[iset]),Form("TailFit-%.f-PPM",runPPM[iset]));
  canTail->SetLogy();
  hTail->Draw("");
  ftail->Draw("same");


  fout->Write();
  return;

}
