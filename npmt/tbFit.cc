
using namespace TMath;
static double xstart=1.05;
static double xstop=10;
static double tiny =.000021;

static double afterLow=1.40;
static double afterHigh=2.2;

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


void tbFit() {

  TFile *inFile = new TFile("tbAnaOut.root","READONLY");
  if(!inFile) return;
  TFile *fout = new TFile("tbFitOut.root","RECREATE");

  enum {MAXSETS=6};

  double runPPM[MAXSETS]={0,1,2,5,10,10};

   
  TString runTag[MAXSETS];
  runTag[0]=TString("00PPM-MU");
  runTag[1]=TString("01PPM-Ran");
  runTag[2]=TString("02PPM-Ran");
  runTag[3]=TString("05PPM-Ran");
  runTag[4]=TString("10PPM-Ran");
  runTag[5]=TString("10PPM-Mu");

  TString runRange[MAXSETS];
  runRange[0]=TString("3000-3020");
  runRange[1]=TString("20005-20020");
  runRange[2]=TString("20025-20040");
  runRange[3]=TString("20045-20060");
  runRange[4]=TString("20065-20080");
  runRange[5]=TString("20205-20215");

  TH1D *hMpvSet[MAXSETS];
  TH1D *hLifeSum[MAXSETS];
  TH1D *hMpvSum[MAXSETS];
  int pcolor[MAXSETS] = {2,3,4,45,1};
  int pstyle[MAXSETS] = {21,22,23,24,25};
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
    if(hname.Contains("LifeSum")) hLifeSum[iset3++]=h;
  }

  printf(" got %i %i %i \n",iset1,iset2,iset3);

  
  TCanvas *canSet = new TCanvas("dataSets","dataSets");
  canSet->Divide(2,3);
  for(int i=0; i<MAXSETS; ++i) {
    canSet->cd(i+1);
    cout << i << "  " << hMpvSet[i]->GetName() << endl;
    hMpvSet[i]->SetTitle(Form(" fitted Landau MPV for set %i runs %s  ",i,runRange[i].Data()));
    hMpvSet[i]->GetXaxis()->SetTitle(" fitted Landau MPV for run ");
    hMpvSet[i]->Draw();
  }
  canSet->Print(".pdf");

  return;

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
    hLifeSum[iset]->SetMarkerSize(0.3);

  }
  canLifeAll->Print(".pdf");


  TCanvas *canNone = new TCanvas("none","none");


  /** fit landau **/
  TF1 *flan[MAXSETS];
  for(int iset=0; iset<MAXSETS; ++iset) {
      flan[iset] = new TF1( Form("fitLan-%i",iset),"landau",100,5000);
      hMpvSum[iset]->Fit(flan[iset],"R");
  }


  std::vector<double> vppm;
  std::vector<double> vppmErr;
  std::vector<double> vmpv;
  std::vector<double> vmpvErr;
  std::vector<double> vwid;
  std::vector<double> vwidErr;

  std::vector<double> vppmMu;
  std::vector<double> vppmErrMu;
  std::vector<double> vmpvMu;
  std::vector<double> vmpvErrMu;
  std::vector<double> vwidMu;
  std::vector<double> vwidErrMu;



  for(int iset=0; iset<MAXSETS; ++iset) {
    printf(" landau fit set %i %s \n",iset+1,hMpvSum[iset]->GetName());
    if(iset==0||iset==5) {
      vppmMu.push_back(runPPM[iset]);
      vppmErrMu.push_back(0.05*runPPM[iset]);
      vmpvMu.push_back(flan[iset]->GetParameter(1));
      vmpvErrMu.push_back(flan[iset]->GetParError(1));
      vwidMu.push_back(flan[iset]->GetParameter(2));
      vwidErrMu.push_back(flan[iset]->GetParError(2));
    } else {
      vppm.push_back(runPPM[iset]);
      vppmErr.push_back(0.05*runPPM[iset]);
      vmpv.push_back(flan[iset]->GetParameter(1));
      vmpvErr.push_back(flan[iset]->GetParError(1));
      vwid.push_back(flan[iset]->GetParameter(2));
      vwidErr.push_back(flan[iset]->GetParError(2));

    }

    for(int ii=0; ii<3; ++ii) {
      printf("\t  param %i %s %.3f +/- %.3f \n",ii,flan[iset]->GetParName(ii),flan[iset]->GetParameter(ii),flan[iset]->GetParError(ii));
    }
  }

 
  double ratio,eratio;
  ratioE(flan[5]->GetParameter(1),flan[5]->GetParError(1),flan[0]->GetParameter(1),flan[0]->GetParError(1),ratio,eratio);

  printf(" ... ratio 10ppm %.3f +/- %.3f to 00ppm %.3f +/- %.3f = %.3f +/- %.3f \n",
      flan[5]->GetParameter(1),flan[5]->GetParError(1),flan[0]->GetParameter(1),flan[0]->GetParError(1),ratio,eratio);

  double mean0 = hMpvSet[0]->GetMean(); //double meanErr0 = hMpvSet[0]->GetMeanError();
  double meanErr0 = hMpvSet[0]->GetRMS();
  double mean5 = hMpvSet[5]->GetMean(); //double meanErr5 = hMpvSet[5]->GetMeanError();
  double meanErr5 = hMpvSet[5]->GetRMS();


  double mratio,meratio;
  ratioE(mean5,meanErr5,mean0,meanErr0,mratio,meratio);

  printf(" ... ratio 10ppm %.3f +/- %.3f to 00ppm %.3f +/- %.3f = %.3f +/- %.3f \n",
      mean5,meanErr5,mean0,meanErr0,mratio,meratio);



  TCanvas *canMpv[MAXSETS];
  gStyle->SetOptFit();
  for(int iset=0; iset<MAXSETS; ++iset) {
    canMpv[iset] = new TCanvas(Form("theMpvSet-%i-%s-%s",iset+1,runTag[iset].Data(),runRange[iset].Data()),Form("theMpvSet %i %s %s",iset+1,runTag[iset].Data(),runRange[iset].Data()));
    cout << iset << "  " << hMpvSum[iset]->GetBinContent(hMpvSum[iset]->GetMaximumBin() ) << endl;
    //hMpvSum[iset]->GetYaxis()->SetRangeUser(0, 1.1*hMpvSum[iset]->GetBinContent(hMpvSum[iset]->GetMaximumBin() ) );
    hMpvSum[iset]->GetYaxis()->SetRangeUser(0.000001,0.015);
    //canMpv[iset]->SetLogy();
    hMpvSum[iset]->SetMarkerStyle(setStyle[iset]);
    hMpvSum[iset]->SetMarkerSize(0.1);
    hMpvSum[iset]->SetMarkerColor(setColor[iset]);
    hMpvSum[iset]->SetLineColor(setColor[iset]);
    hMpvSum[iset]->Draw();
    flan[iset]->Draw("sames");
    canMpv[iset]->Print(".pdf");
  }

  TCanvas *canMpvAll=new TCanvas("MPV-ALL","MPV-ALL");
  canMpvAll->Divide(2,3);
  for(int iset=0; iset<MAXSETS; ++iset) {
    canMpvAll->cd(iset+1);
    hMpvSum[iset]->Draw("");
  }
  canMpvAll->Print(".pdf");


  TGraphErrors *gmpv = new TGraphErrors(vppm.size(),&vppm[0],&vmpv[0],&vppmErr[0],&vmpvErr[0]);
  TGraphErrors *gmpvMu = new TGraphErrors(vppmMu.size(),&vppmMu[0],&vmpvMu[0],&vppmErrMu[0],&vmpvErrMu[0]);
  TMultiGraph *mgmp = new TMultiGraph();

  TCanvas *cmpv = new TCanvas("mpvSum","mpvSum");
  cmpv->SetGridx(); cmpv->SetGridy();
  gmpv->SetTitle("Landau MPV");
  gmpv->SetMarkerColor(kRed);
  gmpv->SetMarkerStyle(22);

  gmpvMu->SetMarkerColor(kBlue);
  gmpvMu->SetMarkerStyle(21);
  gmpv->SetMarkerSize(1);
  gmpv->GetXaxis()->SetTitle(" Xe PPM ");
  gmpv->GetYaxis()->SetTitle(" Landau MPV ");
  mgmp->Add(gmpv,"p");
  mgmp->Add(gmpvMu,"p");
  mgmp->SetTitle("Landau MPV; Xe PPM ; Landau MPV ");
  mgmp->Draw("a");
  cmpv->Print(".png");



  TGraphErrors *gwid = new TGraphErrors(vppm.size(),&vppm[0],&vwid[0],&vppmErr[0],&vwidErr[0]);
  TGraphErrors *gwidMu = new TGraphErrors(vppmMu.size(),&vppmMu[0],&vwidMu[0],&vppmErrMu[0],&vwidErrMu[0]);
  TMultiGraph *mgwid = new TMultiGraph();

  TCanvas *cwid = new TCanvas("widSum","widSum");
  cwid->SetGridx(); cwid->SetGridy();
  gwid->SetTitle("Landau width");
  gwid->SetMarkerColor(kRed);
  gwid->SetMarkerStyle(22);
  gwidMu->SetMarkerColor(kBlue);
  gwidMu->SetMarkerStyle(21);

  gwid->SetMarkerSize(1);
  gwid->GetXaxis()->SetTitle(" Xe PPM ");
  gwid->GetYaxis()->SetTitle(" Landau width ");
  mgwid->Add(gwid,"p");
  mgwid->Add(gwidMu,"p");
  mgwid->SetTitle("Landau Width; Xe PPM ; Landau Width ");
  mgwid->Draw("a");
  cwid->Print(".png");

  return;
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


  for(int iset=0; iset<MAXSETS; ++iset) {
    double k1= kvalFit(runPPM[iset]);
    Double_t integral = hLifeSum[iset]->Integral();
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
      fmodelFit[iset][ifit]->Print();
    }

    canFmodel[iset] = new TCanvas(Form("fmodel-%.f-PPM",runPPM[iset]),Form("fmodel-%.f-PPM",runPPM[iset]));
    canFmodel[iset]->SetLogy();
    fmodelFit[iset][4]->Draw("");
    for(int ifit=0; ifit<MAXSETS; ++ifit) {
      fmodelFit[iset][ifit]->Draw("same");
    }
  }


  TF1* fp[MAXSETS];
  bool PMODEL=false;
  int pfits=5;
  int npars=0;
  for(int iset=0; iset<pfits; ++iset) {
    double k1= kvalFit(runPPM[iset]);
    Double_t integral = hLifeSum[iset]->Integral();
    double A0=(1.-gfrac)*integral;
    double G0 = gfrac*integral;
    if(PMODEL) {
      npars=10;
      fp[iset]= new TF1(Form("Lfit-%.f-%i",runPPM[iset],iset),pmodel,xstart,xstop,npars);
      fp[iset]->SetNpx(1000); // numb points for function
      fp[iset]->SetParNames("binwidth","S0","T0","G0","tauS","tauT","k1","tauG","tauX","type");
      fp[iset]->SetParameters(binwidth,Single0,A0,G0,pmodtS,pmodtT,k1,pmodtG,pmodtX,4);
      fp[iset]->FixParameter(0,binwidth);
      fp[iset]->FixParameter(1,0);
      fp[iset]->SetParLimits(2,0,2*integral);
      fp[iset]->FixParameter(4,pmodtS);
      fp[iset]->FixParameter(5,pmodtT);
      //fp[iset]->FixParameter(7,pmodtG);
      fp[iset]->FixParameter(8,pmodtX);
      fp[iset]->FixParameter(9,4);
    } else {
      npars=8;
      fp[iset] = new TF1(Form("LfFit-%.f-%i",runPPM[iset],iset),fmodel,xstart,xstop,npars);
      fp[iset]->SetParNames("binwidth","A0","tauA","k1","tauX","Am","tauM","type");
      fp[iset]->SetParameters(binwidth,A0, pmodtT, k1 ,pmodtX,G0,pmodtG, 1 );
      fp[iset]->SetNpx(1000); // numb points for function
      fp[iset]->FixParameter(0,binwidth);
      fp[iset]->FixParameter(2,pmodtT);
      fp[iset]->SetParLimits(3,0,2*k1);
      //fp[iset]->FixParameter(4,pmodtX);
      //fp[iset]->FixParameter(5,0);
      fp[iset]->FixParameter(6,pmodtG);
      fp[iset]->FixParameter(7,4);
    }

    printf("fit number %i  %s nbins %i  xstart %E xstop %E integral %E width %E \n",
        iset,hLifeSum[iset]->GetName(),hLifeSum[iset]->GetNbinsX(),xstart,xstop,integral,binwidth); //,hsum1->GetNa

    fp[iset]->SetLineColor(pcolor[iset]);
    fp[iset]->SetMarkerStyle(pstyle[iset]);

    //hLifeSum[iset]->Fit(fp[iset],"RLE0+");
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
  }

  TCanvas *canFit[MAXSETS];
  gStyle->SetOptStat();
  gStyle->SetOptFit();
  for(int iset=0; iset<pfits; ++iset) {
    canFmodel[iset] = new TCanvas(Form("LifeFit-%.f-PPM",runPPM[iset]),Form("LifeFit-%.f-PPM",runPPM[iset]));
    canFmodel[iset]->SetLogy();
    hLifeSum[iset]->Draw("p");
    //fp[iset]->Draw("sames");
  }

  fout->Write();

  /*******/
 
  TCanvas *canLife[MAXSETS];
  for(int iset=0; iset<MAXSETS; ++iset) {
    canLife[iset] = new TCanvas(Form("theLifeSet-%i-%s-%s",iset+1,runTag[iset].Data(),runRange[iset].Data()),Form("theLifeSet %i %s %s",iset+1,runTag[iset].Data(),runRange[iset].Data()));
    for(int ibin =0; ibin< hLifeSum[iset]->GetNbinsX() ; ++ibin) 
      if(hLifeSum[iset]->GetBinContent(ibin)<=0) printf(" %s %i %E \n",hLifeSum[iset]->GetName(),ibin,hLifeSum[iset]->GetBinContent(ibin));
    hLifeSum[iset]->SetMarkerStyle(setStyle[iset]);
    hLifeSum[iset]->SetMarkerColor(setColor[iset]);
    hLifeSum[iset]->SetLineColor(setColor[iset]);
    hLifeSum[iset]->SetMarkerSize(0.3);
    canLife[iset]->SetLogy();
    hLifeSum[iset]->GetYaxis()->SetRangeUser(1E-1,1E3);//1.5*hLifeSum[iset]->GetBinContent(hLifeSum[iset]->GetMaximumBin() ) );
    hLifeSum[iset]->GetXaxis()->SetRangeUser(0.0,10.);
    hLifeSum[iset]->Draw();
    canLife[iset]->Print(".pdf");
  }

 

}
