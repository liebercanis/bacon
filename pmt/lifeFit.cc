/* uses TBaconEvent class */
#include "lifeFit.hh"
using namespace TMath;

static double xstart =1.4;

static double fmodel(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double xppm = par[1];
  double A0 = binwidth*par[2];
  double ta = par[3];
  double tdp = par[4]/xppm;
  double tx = par[5];
  double Am = binwidth*par[6];
  double tm = par[7];
  double tap = 1./( 1./ta + 1./tdp);
  double t1 = 1./( 1./tdp - 1./tx);
  double t2 = 1./( 1./tap - 1./tx);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t2*Exp(-t/tap) - t1*Exp(-t/tdp)) - C/tdp*(t1-t2)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f =  A+D+X+M;
  //printf(" f %f A %f D %f X %f M %f\n",f,A,D,X,M);
  return f;
}



static double fmodelA(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double xppm = par[1];
  double A0 = binwidth*par[2];
  double ta = par[3];
  double tdp = par[4]/xppm;
  double tx = par[5];
  double Am = binwidth*par[6];
  double tm = par[7];
  double tap = 1./( 1./ta + 1./tdp);
  double t1 = 1./( 1./tdp - 1./tx);
  double t2 = 1./( 1./tap - 1./tx);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t2*Exp(-t/tap) - t1*Exp(-t/tdp)) - C/tdp*(t1-t2)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);

  return  A;
}

static double fmodelD(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double xppm = par[1];
  double A0 = binwidth*par[2];
  double ta = par[3];
  double tdp = par[4]/xppm;
  double tx = par[5];
  double Am = binwidth*par[6];
  double tm = par[7];
  double tap = 1./( 1./ta + 1./tdp);
  double t1 = 1./( 1./tdp - 1./tx);
  double t2 = 1./( 1./tap - 1./tx);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t2*Exp(-t/tap) - t1*Exp(-t/tdp)) - C/tdp*(t1-t2)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);

  return D;
}

static double fmodelX(double *xx, double *par) 
{
  double  t= xx[0]-xstart;
  double binwidth = par[0];
  double xppm = par[1];
  double A0 = binwidth*par[2];
  double ta = par[3];
  double tdp = par[4]/xppm;
  double tx = par[5];
  double Am = binwidth*par[6];
  double tm = par[7];
  double tap = 1./( 1./ta + 1./tdp);
  double t1 = 1./( 1./tdp - 1./tx);
  double t2 = 1./( 1./tap - 1./tx);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t2*Exp(-t/tap) - t1*Exp(-t/tdp)) - C/tdp*(t1-t2)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);

  return X;
}






static double fexp(double *xx, double *par)
{
  double x = xx[0];
  double binwidth = par[2];
  double tau = par[0];
  double f = par[1]/tau*binwidth*Exp(-x/tau);
  return f; 
}



static double ftriple(double *xx, double *par)
{
  double x= xx[0];
  double bw = par[0];
  double A1 = par[1];
  double A2 = par[2];
  double A3 = par[3];
  double ts = par[4];
  double td = par[5];
  double tl = par[6];
  double f =  bw*( A1/ts*Exp(-x/ts) - A2/td*Exp(-x/td) + A3/tl*Exp(-x/tl) );
  return f; 
}



lifeFit::lifeFit(Long64_t nloop)
{
  TString tag;
  tag = TString("XenonDoping10ppm_1ppmN2_30451");
  //tag = TString("XenonDoping10ppmGasTest_40000");

  TString fileName ; fileName.Form("TBacon_%s_Ev_0_derivative_3.50.root",tag.Data());

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(!_file0) {
    cout << " not finding file " << fileName << endl;
    return;
  }
  TTree *tree;
  _file0->GetObject("TBacon",tree);
  TBaconEvent *bEvent = new TBaconEvent();
  tree->SetBranchAddress("bevent",&bEvent);


  cout << " TBacon has " << tree->GetEntries() << endl;

  TFile *fout = new TFile(Form("life-%s.root",tag.Data()),"RECREATE");

  int ipmt=0;
  double maxLife=10.0;
  double qCutLow=20;
  double qCut=100;

  TH1D *hLife[3];

  TH2D* hQTime = new TH2D("QTime"," charge by time ",200,0,10,200,0,20);
  hQTime->GetXaxis()->SetTitle(" micro-seconds ");
  hQTime->GetYaxis()->SetTitle(" charge  ");

  hLife[0] = new TH1D("LifeAll"," lifetime PMT cut q<20 ",500,0,maxLife);
  hLife[0]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[0]->SetMarkerColor(kBlack);
  hLife[0]->SetMarkerStyle(22);
  hLife[0]->SetMarkerSize(.2);


  hLife[1] = new TH1D("LifePass",Form(" lifetime PMT cut q < %.0f",qCut),500,0,maxLife);
  hLife[1]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[1]->SetMarkerColor(kBlue);
  hLife[1]->SetMarkerStyle(21);
  hLife[1]->SetMarkerSize(.2);

  
  hLife[2] = new TH1D("LifeFail",Form(" lifetime PMT cut > %.0f",qCutLow),500,0,maxLife);
  hLife[2]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[2]->SetMarkerColor(kBlue);
  hLife[2]->SetMarkerStyle(23);
  hLife[2]->SetMarkerSize(.2);


  Double_t binwidth = hLife[0]->GetBinWidth(1);


  double tauD=0.5;
  TF1* modelFit = new TF1("modelFit",fmodel,xstart,10,8);
  modelFit->SetNpx(1000); // numb points for function
  modelFit->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFit->SetParameters(binwidth,10, 1.0 , 1.8, tauD , 0.02,1.0, 3 );
  modelFit->FixParameter(0,binwidth);
  modelFit->FixParameter(1,10);

  TF1* modelFitA = new TF1("modelFitA",fmodelA,xstart,10,8);
  modelFitA->SetNpx(1000); // numb points for function
  modelFitA->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFitA->SetParameters(binwidth,10, 1.0 , 1.8, tauD , 0.02,1.0, 3  );
  modelFitA->FixParameter(0,binwidth);
  modelFitA->FixParameter(1,10);

  TF1* modelFitD = new TF1("modelFitD",fmodelD,xstart,10,8);
  modelFitD->SetNpx(1000); // numb points for function
  modelFitD->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFitD->SetParameters(binwidth,10, 1.0 , 1.8, tauD , 0.02,1.0, 3  );
  modelFitD->FixParameter(0,binwidth);
  modelFitD->FixParameter(1,10);

  TF1* modelFitX = new TF1("modelFitX",fmodelX,xstart,10,8);
  modelFitX->SetNpx(1000); // numb points for function
  modelFitX->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFitX->SetParameters(binwidth,10, 1.0 , 1.8, tauD , .02,1.0, 3  );
  modelFitX->FixParameter(0,binwidth);
  modelFitX->FixParameter(1,10);


 /* 
  for(int it=0; it<100; ++it) {
    double tt = xstart + double(it)/10.0;
    printf(" modelFit t=%0.2f model %f  A %f D %f X %f\n",tt,modelFit->Eval(tt),modelFitA->Eval(tt),modelFitD->Eval(tt),modelFitX->Eval(tt));
  }
  */


  modelFit->SetLineColor(kBlack);
  modelFit->GetHistogram()->GetYaxis()->SetRangeUser(1E-5,0.1);
  modelFit->GetHistogram()->GetXaxis()->SetRangeUser(0,4);
  modelFitA->SetLineColor(kBlue);
  modelFitD->SetLineColor(kRed);
  modelFitX->SetLineColor(kGreen);

  TCanvas *canFullFit = new TCanvas("modelFit","modelFit");
  canFullFit->SetLogy();
  modelFit->Draw("");
  modelFitA->Draw("same");
  modelFitD->Draw("same");
  modelFitX->Draw("same");
  

  double timeCut=5;
  TH1D *hQSpe = new TH1D("QSpe"," late pulse charge  ",100,0,2);
  hQSpe->GetXaxis()->SetTitle(Form(" late (>%.1f) microsec) pulse  charge ",timeCut));
  //ntEvent->Draw("qspe0>>QSpe","nspe0>0");
  //
  TH1D *hQSpeCut = new TH1D("QSpeCut"," SPE charge XenonDoping10ppm_1ppmN2_30452 ",100,0,2);
  hQSpeCut->GetXaxis()->SetTitle(Form(" q>0.2 and  late (>%.1f) microsec) pulse  charge ",timeCut));
  //hQSpeCut->GetXaxis()->SetTitle(" charge q SPELED_500 ");
  //ntEvent->Draw("qspe0>>QSpe","qspe0>0");


  /*
  TH1D *hQSumE = new TH1D("QSumE"," total event charge  ",150,0,300);
  hQSumE->GetXaxis()->SetTitle(" sum pulse Q for event ");
  hQSumE->SetLineColor(kRed);
  */

  TH2D* hQsumVsN = new TH2D("QsumVsN"," total pulse charge VS # pulses ",100,0,100,500,0,100);
  hQsumVsN->GetXaxis()->SetTitle(" number pulses in event ");
  hQsumVsN->GetYaxis()->SetTitle(" charge "); 


  TH1D* hQRatio = new TH1D("QRatio"," q900/qtot ",100,0,1);
  hQRatio->GetXaxis()->SetTitle(" q900/qtot ");
  hQRatio->GetYaxis()->SetTitle(" events/bin ");  
  
  TH1D *hQSum = new TH1D("QSum"," total event charge  ",60,0,300);
  hQSum->GetXaxis()->SetTitle(" sum pulse Q (x10^9) for event ");

  TH1D *hTriggerTime = new TH1D("TriggerTime"," 3 Qspe Trigger time ",1000,.5,1.5);
  TH1D *hTriggerTime0 = new TH1D("TriggerTime0"," first hit Trigger time ",1000,.5,1.5);


  if(nloop==0) nloop = tree->GetEntries();

  double qspe=0.64;
  for(Long64_t entry =0; entry< nloop ; ++ entry) {
    tree->GetEntry(entry);
    if(bEvent->hits.size()==0) continue;
    double evn = bEvent->npulse;
    hQsumVsN->Fill(bEvent->qsum,evn);
    hQRatio->Fill(bEvent->q900/bEvent->qsum);
    hQSum->Fill(bEvent->qsum);
    double hitqsum = bEvent->qsum;
    if(entry%100==0) printf(" ... %lld nhits %lu \n",entry,bEvent->hits.size());

    if(1) continue;
    // find start time 
    double triggerTime=0;
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      if(bEvent->hits[ip].q>3*qspe) {
        triggerTime = bEvent->hits[ip].time; 
        break;
      }
    }
    hTriggerTime->Fill(triggerTime);
    if( bEvent->hits.size()>0 ) hTriggerTime0->Fill(bEvent->hits[0].time);


    // loop over pulses
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      //cout << entry << "  " << hitq << endl;
      //if(theHit.order==0) triggerTime = theHit.time;
      double hitTime = bEvent->hits[ip].time - triggerTime + 1;
      double hitq = bEvent->hits[ip].q;
      double hitqerr = bEvent->hits[ip].qerr;
      if(hitTime>timeCut) hQSpe->Fill(hitq);
      if(hitTime>timeCut) hQSpeCut->Fill(hitq);
      int hitBin =  hLife[0]->FindBin(hitTime); 
      hQTime->Fill(hitTime,hitq);
      if(hitqsum>20) {
        hLife[0]->SetBinContent(hitBin, hLife[0]->GetBinContent(hitBin)+hitq);
        hLife[0]->SetBinError(hitBin, sqrt( pow(hLife[0]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      }
      if(hitqsum>qCut) {
        hLife[1]->SetBinContent(hitBin, hLife[1]->GetBinContent(hitBin)+hitq);
        hLife[1]->SetBinError(hitBin, sqrt( pow(hLife[1]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      } 
      if(hitqsum<20) {
        hLife[2]->SetBinContent(hitBin, hLife[2]->GetBinContent(hitBin)+hitq);
        hLife[2]->SetBinError(hitBin, sqrt( pow(hLife[2]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      }
    }
  }

  TCanvas *canqn = new TCanvas("qsumVn","qsumVn");
  hQsumVsN->Draw();
  TCanvas *canqratio = new TCanvas("qratio","qratio");
  hQRatio->Draw();

  /*
     TF1 *flanSum = new TF1("flanSum", "expo(0)+landau(2)",20, 300);
     hQSum->Fit("flanSum","R");
     */

  //gStyle->SetOptFit();
  TCanvas *canqsum = new TCanvas("qsum","qsum");
  hQSum->Draw();



  qspe = hQSpeCut->GetMean();
  printf(" spe mean %f error %f \n",qspe,hQSpeCut->GetMeanError());

  /*
  TH1D *hQSumE = new TH1D("QSumE"," total event charge  ",150,0,300);
  hQSumE->GetXaxis()->SetTitle(" sum pulse Q for event ");
  hQSumE->SetLineColor(kRed);
  for(Long64_t entry =0; entry< nloop ; ++ entry) {
    treeHit->GetEntry(entry);
    //cout << entry << "  " << hitq << endl;
    if(theHit.order!=0) continue;
    hQSum->Fill(hitqsum/qspe);
  }
  */

  
  TCanvas *canspe = new TCanvas("qspe","qpe");
  gStyle->SetOptFit();
  TF1 *flan = new TF1("flan", "gaus", 0.2, 2);
  //TF1 *flan = new TF1("flan", "gaus", 0, .4);
  hQSpeCut->Fit("flan","R");
  hQSpeCut->Draw();
 
 
  TF1 *fp[3];  
 
  for(int ih=0; ih<3; ++ih) {
    if(hLife[ih]==NULL) continue;
    Double_t tguess  = (hLife[ih]->GetBinWidth(1))*Double_t(hLife[ih]->GetNbinsX())/10;

    printf(" fit %i \n",ih);
    //hLife[ih]->Rebin(2);
    // fit
    Int_t lowBin = hLife[ih]->FindBin(xstart);
    //if(ih==3) lowBin = hLife[ih]->FindBin(1);
    Double_t xlow= xstart; //hLife[ih]->GetBinLowEdge(lowBin);
    Int_t highBin = hLife[ih]->FindBin(10.0);
    Double_t xhigh=hLife[ih]->GetBinLowEdge(highBin);
    Double_t integral = hLife[ih]->Integral();

    printf(" fit range xlow %f bin %i xhigh %f bin %i \n",xlow,lowBin,xhigh,highBin);
    
    fp[ih] = new TF1("modelFit",fmodel,xlow,xhigh,8);
    fp[ih]->SetNpx(1000); // numb points for function
    fp[ih]->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
    fp[ih]->SetParameters(binwidth,10, integral ,1.8, 2 , 0.02, 1 , 3 );
    fp[ih]->FixParameter(0,binwidth);
    fp[ih]->FixParameter(3,1.787);
    fp[ih]->SetParLimits(6,0,integral);
    //fp[ih]->SetParLimits(7,0,100);


    /*
    fp[ih] = new TF1(Form("fexp-%i",ih),fexp,xlow,xhigh,3);
    fp[ih]->SetNpx(1000); // numb points for function
    fp[ih]->SetParNames("lifetime","integral","binwidth");
    fp[ih]->SetParameters(tguess,integral,binwidth);
    fp[ih]->FixParameter(2,binwidth);
    */

    /*
    fp[ih] = new TF1(Form("ftriple-%i",ih),ftriple,xlow,xhigh,7);
    fp[ih]->SetNpx(1000); // numb points for function
    fp[ih]->SetParNames("BW","AS","AD","AL","TS","TD","TL");
    fp[ih]->SetParameters(binwidth,0.6*integral,0.1*integral,0.3*integral,tguess,tguess,tguess);
    //fp[ih]->FixParameter(4,0.45);
    */

    printf(" \t number %i  %s nbins %i lowbin %i xlow %E xhigh %E integral %E width %E \n",
        ih,hLife[ih]->GetName(),hLife[ih]->GetNbinsX(),lowBin,xlow,xhigh,integral,binwidth); //,hsum1->GetNa
    printf(" pmt %i tau = %f +- %f micro-sec \n",0,fp[0]->GetParameter(0)*1E6,fp[0]->GetParError(0)*1E6);
    hLife[ih]->Fit(fp[ih],"RLE");
    fp[ih]->Print();
  }

  TCanvas *can[3];
  gStyle->SetOptFit();

  hLife[0]->GetYaxis()->SetRangeUser(0.1, hLife[0]->GetBinContent( hLife[0]->GetMaximumBin())*1.1);
  hLife[1]->GetYaxis()->SetRangeUser(0.1,hLife[1]->GetBinContent( hLife[1]->GetMaximumBin())*1.1);
  hLife[2]->GetYaxis()->SetRangeUser(0.1,hLife[2]->GetBinContent( hLife[2]->GetMaximumBin())*1.1);


  can[0] = new TCanvas(Form("life-%s-pass%.0f",tag.Data(),qCutLow),Form("life-%s-pass%.0f",tag.Data(),qCutLow));
  hLife[0]->SetTitle( Form("life-%s-pass%.0f",tag.Data(),qCutLow) );
  can[0]->SetLogy();
  hLife[0]->Draw("E1");
  can[0]->Print(".png");


  can[1] = new TCanvas(Form("life-%s-pass%.0f",tag.Data(),qCut),Form("life-%s-pass%.0f",tag.Data(),qCut));
  can[1]->SetLogy();
  hLife[1]->SetTitle( Form("life-%s-pass%.0f",tag.Data(),qCut) );
  hLife[1]->Draw("E1");


  can[2] = new TCanvas(Form("life-%s-fail%.0f",tag.Data(),qCut),Form("life-%s-fail%.0f",tag.Data(),qCut));
  hLife[2]->SetTitle( Form("life-%s-fail%.0f",tag.Data(),qCut) );
  can[2]->SetLogy();
  hLife[2]->Draw("E1");


  fout->ls();
  fout->Write();

}
