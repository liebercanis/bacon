using namespace TMath;
/** time unit is microseconds **/
static double xstart=1;
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
  double cnorm = Exp(-xstart/tau);
  double binwidth = par[2];
  double f = binwidth/cnorm*par[1]/tau*Exp(-x/tau);
  return f; 
}



void fitTest(int nsbin=8,unsigned MAXGEN=1000000) {


  enum {MAXSETS=6};
  double runPPM[MAXSETS]={0,1,2,5,10,10};
  TString runTag[MAXSETS];
  runTag[0]=TString("00PPM-MU");
  runTag[1]=TString("01PPM-Ran");
  runTag[2]=TString("02PPM-Ran");
  runTag[3]=TString("05PPM-Ran");
  runTag[4]=TString("10PPM-Ran");
  runTag[5]=TString("10PPM-Mu");


  TFile *fout = new TFile("fitTestOut.root","RECREATE");

  Double_t binwidth = double(nsbin)/1000.;
  double maxLife=10.;
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
  double k1= kvalFit(runPPM[iset]);
  double tap = 1./( 1./tauA + k1);

  double cnorm = 1;//Exp(-xstart/tap);
  double Astart = A0*binwidth/cnorm/tap;
  double Gstart = G0*binwidth/cnorm/tauM;


  TH1D* hLifeTest;
  hLifeTest = new TH1D("LifeTest",Form(" lifetime %.2f",tap),lifeBins,0,maxLife);
  hLifeTest->GetXaxis()->SetTitle(" micro-seconds ");
  hLifeTest->SetMarkerColor(kBlack);
  hLifeTest->SetLineColor(kBlack);
  hLifeTest->SetMarkerStyle(22);
  hLifeTest->SetMarkerSize(.2);

  TH1D *hLife128 = new TH1D("Life128",Form(" lifetime 175 %.2f",tripletT),lifeBins,0,maxLife);
  hLife128->GetXaxis()->SetTitle(" micro-seconds ");
  hLife128->SetMarkerColor(kRed);
  hLife128->SetLineColor(kRed);
  hLife128->SetMarkerStyle(23);
  hLife128->SetMarkerSize(.2);

  TH1D *hLife175 = new TH1D("Life175",Form(" lifetime 128 %.2f",pmodtX),lifeBins,0,maxLife);
  hLife175->GetXaxis()->SetTitle(" micro-seconds ");
  hLife175->SetMarkerColor(kBlue);
  hLife175->SetLineColor(kBlue);
  hLife175->SetMarkerStyle(24);
  hLife175->SetMarkerSize(.2);


  TF1* fp = new TF1(Form("fexp-%.2f",tripletT),fexp,xstart,xstop,3);
  fp->SetNpx(1000); // numb points for function
  fp->SetParNames("lifetime","integral","binwidth");
  fp->SetParameters(tripletT,integral,binwidth);
  fp->FixParameter(2,binwidth);




  TF1 *fmodelFit[6];
  for(int jfit=0; jfit<6; ++jfit) {
    fmodelFit[jfit] = new TF1(Form("fmodelFit-%.f-%i",runPPM[iset],jfit),fmodel,xstart,xstop,8);
    fmodelFit[jfit]->SetParNames("binwidth","A0","tauA","k1","tauX","Gfrac","tauM","type");
    fmodelFit[jfit]->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX,G0 ,pmodtG, jfit );
    fmodelFit[jfit]->SetNpx(1000); // numb points for function
    fmodelFit[jfit]->FixParameter(0,binwidth);
    fmodelFit[jfit]->FixParameter(7,jfit);
    fmodelFit[jfit]->SetLineColor(setColor[jfit]);
    fmodelFit[jfit]->SetMarkerStyle(setStyle[jfit]);
    fout->Append(fmodelFit[jfit]);
    if(jfit<5) fmodelFit[jfit]->SetLineStyle(10);

  }

  TF1 *fm=fmodelFit[5];
  for(int ii=0; ii<8; ++ii) {
    if(ii==1) printf("\t  param %i %s %.3E +/- %.3E \n",ii,fm->GetParName(ii),fm->GetParameter(ii),fm->GetParError(ii));
    else printf("\t  param %i %s %.3f +/- %.3f \n",ii,fm->GetParName(ii),fm->GetParameter(ii),fm->GetParError(ii));
  }
  printf("start A %.3E G %.3E \n",Astart,Gstart);

  TCanvas *canFmodelD = new TCanvas("FullModel","FullModel");
  canFmodelD->SetLogy();
  //fmodelFit[5]->GetHistogram()->GetYaxis()->SetRangeUser(tiny,5E12);
  fmodelFit[5]->Draw("");
  for(int jfit=0; jfit<4; ++jfit) {
    fmodelFit[jfit]->Draw("same");
  }


  int ifit=5;

  fm = new TF1(Form("fmodelFit-%i",int(runPPM[iset])),fmodel,xstart,xstop,8);
  fm->SetParNames("binwidth","A0","tauA","k1","tauX","G0","tauG","type");
  fm->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX, G0,pmodtG, ifit );
  fm->SetNpx(1000); // numb points for function
  fm->FixParameter(0,binwidth);
  fm->FixParameter(2,pmodtT);
  fm->FixParameter(4,pmodtX);
  fm->FixParameter(6,pmodtG);
  if(ifit<5) {
    fm->FixParameter(5,G0);
  }
  fm->FixParameter(7,ifit);
  fm->SetLineColor(setColor[ifit]);
  //fm->SetLineStyle(10);
  fm->SetMarkerStyle(setStyle[ifit]);
  fm->Print();


  TCanvas *canFmodel = new TCanvas(Form("fmodel-%.f-PPM",runPPM[iset]),Form("fmodel-%.f-PPM",runPPM[iset]));
  canFmodel->SetLogy();
  fm->Draw("");

  double lightAll = fm->Integral(xstart,xstop)/binwidth;
  double light128 = fmodelFit[0]->Integral(xstart,xstop)/binwidth;
  double light175 = fmodelFit[2]->Integral(xstart,xstop)/binwidth;

  hLifeTest->SetTitle(Form(" tauA %.2f k %.2f tauX %.2f gas frac  %.2f",tripletT,k1,pmodtX,G0/A0));
  int MAX128=  int(double(MAXGEN)*light128/lightAll);
  int MAX175=  int(double(MAXGEN)*light175/lightAll);
  
  printf(" ...... generating all %i 128 %i 175 %i sum %i \n",MAXGEN,MAX128,MAX175,MAX128+MAX175);

  for(unsigned ig=0; ig<MAXGEN ; ++ig) hLifeTest->Fill(fm->GetRandom());
  for(unsigned ig=0; ig<MAX128 ; ++ig) hLife128->Fill(fmodelFit[0]->GetRandom());
  for(unsigned ig=0; ig<MAX175 ; ++ig) hLife175->Fill(fmodelFit[2]->GetRandom());


  TCanvas *can = new TCanvas( Form("Fit-model-%.2f-%i-%.1f-type%i",tap,nsbin,xstart,ifit),Form("Fit-model-%.2f-%i-%.1f",tripletT,nsbin,xstart));
  can->SetLogy();
  gStyle->SetOptFit();
  hLifeTest->Fit(fm,"R");
  hLifeTest->Draw();
  hLife128->Draw("same");
  hLife175->Draw("same");
  hLifeTest->Draw("same");

  double lAll =  hLifeTest->Integral();
  double l128=   hLife128->Integral();
  double l175 =  hLife175->Integral();

  printf(" life %.2f events %i binw %f cnorm %f lightAll %13.3f light128 %13.3f light175 %13.3f 175 ratio =  %.3f \n",tripletT,MAXGEN,binwidth,cnorm,lightAll,light128,light175,
      light175/lightAll);

  printf(" life %.2f events %i binw %f cnorm %f lightAll %13.3f light128 %13.3f light175 %13.3f 175 ratio =  %.3f \n",tripletT,MAXGEN,binwidth,cnorm,lAll,l128,l175,
      l175/lAll);

  for(int ii=0; ii<8; ++ii) {
    if(fm->GetParError(ii)==0) continue;
    if(ii==1) printf("\t  param %i %s %.3E +/- %.3E \n",ii,fm->GetParName(ii),fm->GetParameter(ii),fm->GetParError(ii));
    else printf("\t  param %i %s %.3f +/- %.3f \n",ii,fm->GetParName(ii),fm->GetParameter(ii),fm->GetParError(ii));
  }

  printf(" sum A+G = %E \n",fm->GetParameter(1)+fm->GetParameter(5));
  can->Print(".pdf");

  fout->Write();
  return;

}
