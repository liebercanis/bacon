
using namespace TMath;
static double xstart=1.0;
static double xstop=10;
static double tiny =.000021;

static double afterLow=1.46;
static double afterHigh=1.56;
static double afterLow2=1.95;
static double afterHigh2=2.1;


/*I fix all of the lifetimes: microseconds*/
static double pmodtS = 7E-3;
static double pmodtT = 1.96;
static double pmodtG = 3.2;
static double pmodtX = 20E-3;
static double tailStart=6.0;


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
  T=T/tt;
  S=S/ts;
  G=G/tg;
  X=X/tx;
  double f;
  if(type ==0) f=T+S;
  else if(type ==1) f=M;
  else if(type ==2) f=X;
  else if(type ==3) f=G;
  else f =  S+T+G+X;  
  return f;
}



static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double binwidth = par[2];
  double tau = par[0];
  double cnorm = Exp(-tailStart/tau) - Exp(-xstop/tau);
  double f = par[1]/tau*binwidth/cnorm*Exp(-x/tau);
  return f; 
}



void tbFit1(int iset=0) {

  TFile *inFile = new TFile("tbAnaOut.root","READONLY");
  if(!inFile) return;
  TFile *fout = new TFile("tbFitOut.root","RECREATE");

  enum {MAXSETS=6};

  if(iset>=MAXSETS) {
    printf(" max set is %i \n",MAXSETS); return;
  }

  double runPPM[MAXSETS]={0,1,2,5,10,10};


  TString runTag[MAXSETS];
  runTag[0]=TString("00PPM-MU");
  runTag[1]=TString("01PPM-Ran");
  runTag[2]=TString("02PPM-Ran");
  runTag[3]=TString("05PPM-Ran");
  runTag[4]=TString("10PPM-Ran");
  runTag[5]=TString("10PPM-Mu");

  TString runRange[MAXSETS];
  runRange[0]=TString("03000-03020");
  runRange[1]=TString("20005-20020");
  runRange[2]=TString("20025-20040");
  runRange[3]=TString("20045-20060");
  runRange[4]=TString("20065-20080");
  runRange[5]=TString("20205-20215");

  TH1D *hMpvSet[MAXSETS];
  TH1D *hLifeSum[MAXSETS];
  TH1D *hLifeBlah[MAXSETS];
  TH1D *hMpvSum[MAXSETS];
  int pcolor[MAXSETS] = {kRed,kGreen,kBlue,kYellow,kGreen-6,kBlack};
  int pstyle[MAXSETS] = {21,22,23,24,25,26};
  int setColor[MAXSETS]={1,2,3,4,kYellow-6,kGreen-6};
  int setStyle[MAXSETS]={21,22,23,24,25,26};



  TIter next(inFile->GetListOfKeys());
  TKey *key;
  TH1D* hlife = NULL;
  int iset1=0;
  int iset2=0;
  int iset3=0;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1D *h = (TH1D*)key->ReadObj();
    TString hname(h->GetName());
    cout << hname << endl;
    if(hname.Contains("MpvSet")) hMpvSet[iset1++]=h;
    if(hname.Contains("MpvSum"))  hMpvSum[iset2++]=h;
    if(hname.Contains("LifeSum")) hLifeBlah[iset3++]=h;
  }

  printf(" got %i %i %i \n",iset1,iset2,iset3);

  int ilow = hLifeBlah[0]->FindBin(afterLow);
  int ihigh = hLifeBlah[0]->FindBin(afterHigh);
  int ilow2 = hLifeBlah[0]->FindBin(afterLow2);
  int ihigh2 = hLifeBlah[0]->FindBin(afterHigh2);

   
  printf(" refilling skipping bins from %i (%f) to %i (%f) and from %i (%f) to %i (%f) \n",ilow,afterLow,ihigh,afterHigh,ilow2,afterLow2,ihigh2,afterHigh2);
  for(int i=0; i<MAXSETS; ++i) {
    printf(" set %i events %.0f lifetime integral %.0f \n",i,hMpvSum[i]->GetEntries(),hLifeBlah[i]->Integral());
    hLifeSum[i]= (TH1D*) hLifeBlah[i]->Clone(Form("LifeToFit-%i",i));
    hLifeSum[i]->Reset(); 
    for(int ibin=0; ibin< hLifeBlah[i]->GetNbinsX() ; ++ibin) {
      if(ibin>= ilow && ibin<= ihigh) continue;
      if(ibin>= ilow2 && ibin<= ihigh2) continue;
      hLifeSum[i]->SetBinContent(ibin,hLifeBlah[i]->GetBinContent(ibin));
      hLifeSum[i]->SetBinError(ibin,hLifeBlah[i]->GetBinError(ibin));
    }
  }
 
  // rbin
  for(int iset=0; iset<MAXSETS; ++iset) {printf(" %i %s %s %s \n",iset,hMpvSet[iset]->GetName(),hMpvSum[iset]->GetName(),hLifeSum[iset]->GetName());
    cout << " nbins " << hLifeSum[iset]->GetNbinsX() << " wid  " << hLifeSum[iset]->GetBinWidth(1);
    hLifeSum[iset]->Rebin(5);
    cout << " nbins " << hLifeSum[iset]->GetNbinsX() << " wid  " << hLifeSum[iset]->GetBinWidth(1) << endl;
  }



  TCanvas *canLifeAll=new TCanvas("Life-ALL","Life-ALL");
  canLifeAll->SetLogy();
  hLifeSum[0]->Draw();
  hLifeSum[0]->SetTitle("SPE versus Time all sets ");
  for(int iset=1; iset<MAXSETS; ++iset) {
    hLifeSum[iset]->SetTitle("SPE versus Time all sets ");
    hLifeSum[iset]->Draw("sames");
    hLifeSum[iset]->SetMarkerStyle(setStyle[iset]);
    hLifeSum[iset]->SetMarkerColor(setColor[iset]);
    hLifeSum[iset]->SetLineColor(setColor[iset]);
    hLifeSum[iset]->GetYaxis()->SetRangeUser(1E-1,1E3);
    hLifeSum[iset]->SetMarkerSize(0.3);

  }
  canLifeAll->Print(".pdf");


  TCanvas *canNone = new TCanvas("none","none");


   /***** fitting *****/
  Double_t binwidth = hLifeSum[0]->GetBinWidth(1);

  // pmodel 
  TF1* pmodelFit[MAXSETS][MAXSETS];
   double Single0=0;

  // fmodel 
  double tauX=pmodtX;
  double tauA = pmodtT;
  double tauM = pmodtG; 
  TF1* fmodelFit[MAXSETS][MAXSETS];
  TCanvas *canFmodel[MAXSETS];

  double gfrac = 0.1;


  double k1= kvalFit(runPPM[iset]);

  double integral=0;
  for(int ibin =0; ibin < hLifeSum[iset]->GetNbinsX() ; ++ibin) integral += hLifeSum[iset]->GetBinContent(ibin);
  double integral2 = hLifeSum[iset]->Integral();

  double A0=(1.-gfrac)*integral;
  double G0 = gfrac*integral;
  for(int ifit=0; ifit<MAXSETS; ++ifit) {
    fmodelFit[iset][ifit] = new TF1(Form("fmodelFit-%.f-%i",runPPM[iset],ifit),fmodel,xstart,xstop,8);
    fmodelFit[iset][ifit]->SetParNames("binwidth","A0","tauA","k1","tauX","G0","tauG","type");
    fmodelFit[iset][ifit]->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX, G0 ,pmodtG, ifit );
    fmodelFit[iset][ifit]->SetNpx(1000); // numb points for function
    fmodelFit[iset][ifit]->FixParameter(0,binwidth);
    fmodelFit[iset][ifit]->FixParameter(7,ifit);
    fmodelFit[iset][ifit]->SetLineColor(pcolor[ifit]);
    //fmodelFit[iset][ifit]->SetLineStyle(10);
    fmodelFit[iset][ifit]->SetMarkerStyle(pstyle[ifit]);
  }


  // fit exponential tail
  double tauG=2.89116e+00;
  bool tailFit=true;
  gStyle->SetOptStat();
  gStyle->SetOptFit(1111111);
  double gnorm=0;
  if(tailFit) {
    double tailIntegral=0;
    for(int ibin =  hLifeSum[iset]->FindBin(tailStart) ; ibin <  hLifeSum[iset]->FindBin(xstop) ; ++ibin) tailIntegral += hLifeSum[iset]->GetBinContent(ibin);
    double tailInt2 = hLifeSum[iset]->Integral(hLifeSum[iset]->FindBin(tailStart),hLifeSum[iset]->FindBin(xstop));

    printf(" integral %0.2f (%.2f) tail (%i,%i)  %.2f (%0.2f) \n",integral2,integral,hLifeSum[iset]->FindBin(tailStart),hLifeSum[iset]->FindBin(xstop),
        tailInt2,tailIntegral);

    TF1 * ftail = new TF1(Form("fexp-%i",iset),fexp,tailStart,xstop,3);
    ftail->SetNpx(1000); // numb points for function
    ftail->SetParNames("lifetime","integral","binwidth");
    ftail->SetParameters(3.2,tailIntegral,binwidth);
    //ftail->FixParameter(0,3.2);
    ftail->FixParameter(2,binwidth);
    for(int ii=0; ii<3; ++ii) {
      printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
    }
    double sum = ftail->Integral(tailStart,xstop)/binwidth; 
    cout << " starting sum from function ftail  " << sum << endl;


    //hLifeSum[iset]->Fit(ftail,"RLE0+");
    TH1D* hTail = (TH1D*) hLifeSum[iset]->Clone("Tail");
    //hTail->Reset();
    //for(unsigned ig=0; ig<int(tailIntegral) ; ++ig) hTail->Fill(ftail->GetRandom());
    hTail->Fit(ftail,"RLE0+");

    for(int ii=0; ii<3; ++ii) {
      printf("\t  param %i %s %.3f +/- %.3f \n",ii,ftail->GetParName(ii),ftail->GetParameter(ii),ftail->GetParError(ii));
    }

    sum = ftail->Integral(tailStart,xstop)/binwidth;
    tauG = ftail->GetParameter(0);
    gnorm = sum*(1.0 - Exp(-xstop/tauG))/(Exp(-tailStart/tauG) - Exp(-xstop/tauG));

    cout << " sum from function ftail  " << sum << " tauG =   " << tauG << "  gnorm " << gnorm <<  endl;


    TCanvas *canTail;
    canTail = new TCanvas(Form("TailFit-%.f-PPM",runPPM[iset]),Form("TailFit-%.f-PPM",runPPM[iset]));
    hTail->SetTitle(Form("Tail Fit set %.f PPM %s",runPPM[iset],runTag[iset].Data()));
    hTail->SetLineColor(kBlack);
    hTail->SetMarkerColor(kBlack);

    canTail->SetLogy();
    hTail->Draw("p");
    TPaveStats *st1 =(TPaveStats*) hTail->GetListOfFunctions()->FindObject("stats");
    st1->SetOptFit(1011); //for example
    ftail->Draw("same");
    canTail->Print(".pdf");

  }

  TF1* fp[MAXSETS];
  int pfits=5;
  int npars=0;
  npars=8;
  if(gnorm!=0) G0=gnorm;
  if(iset<3) G0=0;
  fp[iset] = new TF1(Form("LfFit-%.f-%i",runPPM[iset],iset),fmodel,xstart,xstop,npars);
  fp[iset]->SetParNames("binwidth","A0","tauA","k1","tauX","G0","tauG","type");
  fp[iset]->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX,G0,tauG, pfits );
  fp[iset]->SetNpx(1000); // numb points for function
  fp[iset]->FixParameter(0,binwidth);
  fp[iset]->SetParLimits(1,0,2*A0);
  fp[iset]->FixParameter(2,pmodtT);
  fp[iset]->SetParLimits(3,0,10*k1);
  fp[iset]->FixParameter(4,pmodtX);
  //fp[iset]->SetParLimits(5,0,2.*G0);
  fp[iset]->FixParameter(5,G0);
  fp[iset]->FixParameter(6,tauG);
  fp[iset]->FixParameter(7,pfits);

  printf("fit number %i  %s nbins %i  xstart %E xstop %E integral %E width %E \n",
      iset,hLifeSum[iset]->GetName(),hLifeSum[iset]->GetNbinsX(),xstart,xstop,integral,binwidth); //,hsum1->GetNa

  fp[iset]->SetLineColor(pcolor[iset]);
  fp[iset]->SetMarkerStyle(pstyle[iset]);

  hLifeSum[iset]->Fit(fp[iset],"RLE0+");
  fp[iset]->Print();
  printf(" xxxx  %.f PPM \n",runPPM[iset]);
  for(int ii=0; ii<npars; ++ii) {
    printf("\t  param %i %s %.3f +/- %.3f \n",ii,fp[iset]->GetParName(ii),fp[iset]->GetParameter(ii),fp[iset]->GetParError(ii));
  }

  for(int ifit=0; ifit<MAXSETS; ++ifit) {
    fmodelFit[iset][ifit]->SetParameters(binwidth,fp[iset]->GetParameter(1),fp[iset]->GetParameter(2),fp[iset]->GetParameter(3),fp[iset]->GetParameter(4),
        fp[iset]->GetParameter(5),fp[iset]->GetParameter(6),ifit);
  }

  double lightAll = fp[iset]->Integral(xstart,xstop)/binwidth;
  double light128 = fmodelFit[iset][0]->Integral(xstart,xstop)/binwidth;
  double light175 = fmodelFit[iset][2]->Integral(xstart,xstop)/binwidth;
  printf(" set %i lightAll %13.3f light128 %13.3f light175 %13.3f 175 ratio =  %.3f \n",iset,lightAll,light128,light175,light175/lightAll);

  canFmodel[iset] = new TCanvas(Form("fmodel-%.f-PPM",runPPM[iset]),Form("fmodel-%.f-PPM",runPPM[iset]));
  canFmodel[iset]->SetLogy();
  fmodelFit[iset][4]->GetYaxis()->SetRangeUser(1E-1,1E3);
  fmodelFit[iset][4]->Draw("");
  for(int ifit=0; ifit<MAXSETS; ++ifit) {
    fmodelFit[iset][ifit]->SetTitle(Form("fmodel set %i Xe %.f PPM",iset,runPPM[iset]));
    fmodelFit[iset][ifit]->Draw("same");
  }
  canFmodel[iset]->Print(".pdf");


  gStyle->SetOptFit(1111111);
  gStyle->cd();  

  TCanvas *canFit[MAXSETS];
  canFit[iset] = new TCanvas(Form("LifeFit-%.f-PPM-%s",runPPM[iset],runTag[iset].Data()),Form("LifeFit-%.f-PPM-%s",runPPM[iset],runTag[iset].Data()));
  canFit[iset]->SetLogy();
  hLifeSum[iset]->SetTitle(Form("LifeFit set %.f PPM %s",runPPM[iset],runTag[iset].Data()));
  hLifeSum[iset]->GetYaxis()->SetRangeUser(1E-1,1E3);
  hLifeSum[iset]->SetLineColor(kBlack);
  hLifeSum[iset]->SetMarkerColor(kBlack);
  hLifeSum[iset]->Draw("p");
  fp[iset]->SetLineColor(kRed);
  fp[iset]->Draw("sames");
  TPaveStats *st =(TPaveStats*) hLifeSum[iset]->GetListOfFunctions()->FindObject("stats");
  st->SetOptFit(); //for example
  canFit[iset]->Modified();
  canFit[iset]->Update();
  canFit[iset]->Print(".pdf");


  fout->Write();


}
