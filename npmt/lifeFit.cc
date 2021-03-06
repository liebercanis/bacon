/* uses TBaconEvent class */
#include "lifeFit.hh"
using namespace TMath;

static double xstart =1;


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
  double t1=1./( 1./tx- 1./tdp);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t1*Exp(-t/tdp) - t2*Exp(-t/tap)) + C/tdp*(t2-t1)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f =  A+X+M;  
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
  double t1=1./( 1./tx- 1./tdp);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t1*Exp(-t/tdp) - t2*Exp(-t/tap)) + C/tdp*(t2-t1)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f =  A+X+M; // D dark
  //printf(" f %f A %f D %f X %f M %f\n",f,A,D,X,M);
  return A;
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
  double t1=1./( 1./tx- 1./tdp);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t1*Exp(-t/tdp) - t2*Exp(-t/tap)) + C/tdp*(t2-t1)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f =  A+X+M; // D dark
  //printf(" f %f A %f D %f X %f M %f\n",f,A,D,X,M);
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
  double t1=1./( 1./tx- 1./tdp);
  double t2=1./( 1./tx- 1./tap);
  double C = A0*ta/tdp;
  double A = A0 * Exp(-t/tap);
  double D = C*( Exp(-t/tdp) - Exp(-t/tap) );
  double X = C/tdp*(t1*Exp(-t/tdp) - t2*Exp(-t/tap)) + C/tdp*(t2-t1)*Exp(-t/tx);
  double M = Am*Exp(-t/tm);
  double f =  A+X+M; // D dark
  //printf(" f %f A %f D %f X %f M %f\n",f,A,D,X,M);
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



lifeFit::lifeFit(Int_t runstart, Int_t runstop )
{
  bool VmaxCut = false;
  TString tag;
  if(!VmaxCut) tag.Form("%i_%i_DS_2-All",runstart,runstop);
  else tag.Form("%i_%i_DS_2-Mu",runstart,runstop);
  double lifeBins = 400; //1.0E4/8.0; // ns bins
  TChain *tree = new TChain("TBacon");
  TBaconEvent *bEvent = new TBaconEvent();

  TString fileName1("TBAcon_10000_20000_DS_2.root");
  TString fileName2("TBAcon_20001_20100_DS_2.root");
  TString fileName3("TBAcon_20101_20222_DS_2.root");
  TString fileName4("TBAcon_30000_30452_DS_2.root");
  TString fileName5("TBAcon_30440_30443_DS_2.root");

  TString fileName;
  if(runstart>=10000&&runstop<20001) fileName = fileName1;
  else if(runstart>=20001&&runstop<20100) fileName = fileName2;
  else if(runstart>=20101&&runstop<20222) fileName = fileName3;
  else if(runstart>=30000&&runstop<30453) fileName = fileName4;
  else fileName = fileName5;

  printf(" run start %i runstop %i file %s \n",runstart,runstop,fileName.Data());

  tree->Add(fileName);

  tree->SetBranchAddress("bevent",&bEvent);

  tree->GetListOfFiles()->ls();

  cout << " TBacon has " << tree->GetEntries() << endl;

  TFile *fout = new TFile(Form("life-%s-bins-%0.f.root",tag.Data(),lifeBins),"RECREATE");
  TTree *oTree = new TTree("outBacon","output bevent");
  TBaconEvent *oEvent = new TBaconEvent();
  oTree->Branch("outEvent",&oEvent);

  int ipmt=0;
  double maxLife=10.0;
  double qCutLow=5;
  double qCut=100;


  TH1D *hLife[3];

  TH2D* hQTime = new TH2D("QTime"," charge by time ",200,0,10,200,0,20);
  hQTime->GetXaxis()->SetTitle(" micro-seconds ");
  hQTime->GetYaxis()->SetTitle(" charge  ");

  
  TH1D* hLifePeak;
  hLifePeak = new TH1D("LifePeak"," lifetime PMT cut q<20 ",lifeBins,0,maxLife);
  hLifePeak->GetXaxis()->SetTitle(" micro-seconds ");
  hLifePeak->SetMarkerColor(kBlack);
  hLifePeak->SetMarkerStyle(22);
  hLifePeak->SetMarkerSize(.2);


  hLife[0] = new TH1D("LifeAll"," lifetime PMT >10 hits ",lifeBins,0,maxLife);
  hLife[0]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[0]->SetMarkerColor(kBlack);
  hLife[0]->SetMarkerStyle(22);
  hLife[0]->SetMarkerSize(.2);

  TH1D* hLifeResiduals = new TH1D("LifeAllResiduals"," residuals ",100,-5,5);
  TNtuple *nResidual = new TNtuple("nResidual","","bin:x:res:fit:val:error");


  hLife[1] = new TH1D("LifePass",Form(" lifetime PMT cut q < %.0f",qCut),lifeBins,0,maxLife);
  hLife[1]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[1]->SetMarkerColor(kBlue);
  hLife[1]->SetMarkerStyle(21);
  hLife[1]->SetMarkerSize(.2);

  
  hLife[2] = new TH1D("LifeFail",Form(" lifetime PMT cut > %.0f",qCut),lifeBins,0,maxLife);
  hLife[2]->GetXaxis()->SetTitle(" micro-seconds ");
  hLife[2]->SetMarkerColor(kBlue);
  hLife[2]->SetMarkerStyle(23);
  hLife[2]->SetMarkerSize(.2);


 

 /* 
  for(int it=0; it<100; ++it) {
    double tt = xstart + double(it)/10.0;
    printf(" modelFit t=%0.2f model %f  A %f D %f X %f\n",tt,modelFit->Eval(tt),modelFitA->Eval(tt),modelFitD->Eval(tt),modelFitX->Eval(tt));
  }
  */


  double timeCut=5;
  TH1D *hQSpe = new TH1D("QSpe"," late pulse charge  ",1000,0,10);
  hQSpe->GetXaxis()->SetTitle(Form(" late (>%.1f) microsec) pulse  charge ",timeCut));
  //ntEvent->Draw("qspe0>>QSpe","nspe0>0");
  //
  TH1D *hQSpeCut = new TH1D("QSpeCut"," SPE charge XenonDoping10ppm_1ppmN2_30452 ",1000,0,10);
  hQSpeCut->GetXaxis()->SetTitle(Form(" q>0.2 and  late (>%.1f) microsec) pulse  charge ",timeCut));
  //hQSpeCut->GetXaxis()->SetTitle(" charge q SPELED_500 ");
  //ntEvent->Draw("qspe0>>QSpe","qspe0>0");
  //
  TH1D *hMuVmax = new TH1D("MuVmax"," muVmax ",1000,0,1);
  hQSpeCut->GetXaxis()->SetTitle("muon Vmax" );



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
  
  TH1D *hQSum = new TH1D("QSum"," total event charge  ",1000,0,100);
  hQSum->GetXaxis()->SetTitle(" event total dt-Q (x10^9) for event ");
  
  TH1D *hHitResidual = new TH1D("HitResidual"," residual  ",1000,-5,5);

  TH1D *hQSum10 = new TH1D("QSum10"," total event charge > 10 pulses  ",1000,0,100);
  hQSum10->GetXaxis()->SetTitle(" event total dt-Q (x10^9)for event ");

  TH1D *hEventSumCut = new TH1D("EventSumCut"," good summed event charge > 10 pulses  ",1000,0,10000);
  hEventSumCut->GetXaxis()->SetTitle(" event sum dt-Q (x10^9) ");

  
  TH1D *hRunSumCut = new TH1D("RunSumCut"," good summed event charge > 10 pulses  ",1000,0,10000);
  hRunSumCut->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");

  TH1D *hEventSumNeil = new TH1D("EventSumNeil"," good summed event charge > 10 pulses  ",1000,0,10000);
  hEventSumNeil->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");

  TH1D *hRunSumNeil = new TH1D("RunSumNeil"," good summed event charge > 10 pulses  ",1000,0,10000);
  hRunSumNeil->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");




  TH1D *hTriggerTime = new TH1D("TriggerTime"," 3 Qspe Trigger time ",100,.9,1.1);
  TH1D *hTriggerTimeMax = new TH1D("TriggerTimeMax"," first hit Trigger time ",100,.9,1.1);

  TH2D *hTriggerQTime = new TH2D("TriggerQTime"," charge vs time  ",100,.9,1.1,1000,0,200);
  TH1D *hTriggerQ = new TH1D("TriggerQ"," trig charge ",1000,.9,1.1);
  TH1D *hTrigEvent = new TH1D("TrigEvent"," trig charge ",1000,.9,1.1);
  TH1D *hTrigEventQhit = new TH1D("TrigEventQhit"," trig charge ",1000,.9,1.1);
  TH1D *hTrigEventQsum = new TH1D("TrigEventQsum"," trig charge ",1000,.9,1.1);
  TNtuple *nTrig =  new TNtuple("nTrig","","event:triggerTime:qtime:qmax:qtrigSum");
  TNtuple *ntEvent = new TNtuple("ntEvent","","event:trigTime:nhits:maxq:peakQ:peakTime:pre:sum:early:late");
  TNtuple *ntRun = new TNtuple("ntRun","","run:nev:mean:meanErr:sig:sigErr:mpv:mpvErr:wid:widErr:sumSPE");


  Long64_t nloop = tree->GetEntries();
  Long64_t goodTrigger=0;
  //double qspe=0.64;
  double qspe=0.1;
  
  double lowVal=0.9;
  double highVal=1.10;
  double afterLow=1.45;
  double afterHigh=1.6;

  Int_t currentRun=runstart;

  TDirectory *byRun = fout->mkdir("ByRun");
  byRun->cd();
  TH1D* hLifeRun=NULL;
  TH1D* hChargeRun=NULL;

  Long64_t treeNumber = -1;
  for(Long64_t entry =0; entry< nloop ; ++ entry) {
    Long64_t itree = tree->LoadTree(entry);
    if(itree!=treeNumber) {
      tree->SetBranchAddress("bevent",&bEvent);
      treeNumber= itree;
    }
    tree->GetEntry(entry);
    
    if(bEvent->npmt!=0) continue; // only use PMT zero
    qspe = 1E9*bEvent->spe; 

    // throw away low hit number events
    //if(  bEvent->hits.size()<10) continue;

    // find run range
    if(bEvent->run< runstart) continue;
    if(bEvent->run> runstop)  break;
    //if(bEvent->run<30000&&bEvent->run>20199) break;
    //

    // if first run
    if(!hLifeRun) {
      hLifeRun   = (TH1D*) hLife[0]->Clone(Form("LifeRun%i",bEvent->run));
      hLifeRun->SetTitle(Form("Life Run %i",bEvent->run));
      hLifeRun->Reset();
      hChargeRun = (TH1D*) hRunSumCut->Clone(Form("ChargeRun%i",bEvent->run));
      hChargeRun->SetTitle(Form("Summed Event Charge Run %i",bEvent->run));
      hChargeRun->Reset();
    }
    // Landau fits
    if(bEvent->run!=currentRun) {
      if(hChargeRun) printf(" end of run %i events %E totSPE %E \n" , currentRun, hChargeRun->GetEntries(),hLifeRun->Integral());
      hRunSumCut->Fit("landau");
      TF1 *lfit =  hRunSumCut->GetFunction("landau");
      hRunSumNeil->Fit("landau");
      TF1 *lfit2 =  hRunSumNeil->GetFunction("landau");
      if(lfit&&lfit2) ntRun->Fill(float(currentRun),float(hRunSumCut->GetEntries()),lfit->GetParameter(1),lfit->GetParError(1),lfit->GetParameter(2),lfit->GetParError(2),
          lfit2->GetParameter(1),lfit2->GetParError(1),lfit2->GetParameter(2),lfit2->GetParError(2),hLifeRun->Integral());
      // prepare for new run
      currentRun=bEvent->run;
      hRunSumCut->Reset();
      hRunSumNeil->Reset();
      hLifeRun   = (TH1D*) hLife[0]->Clone(Form("LifeRun%i",bEvent->run));
      hLifeRun->SetTitle(Form("Life Run %i",bEvent->run));
      hLifeRun->Reset();
      hChargeRun = (TH1D*) hRunSumCut->Clone(Form("ChargeRun%i",bEvent->run));
      hChargeRun->SetTitle(Form("Summed Event Charge Run %i",bEvent->run));
      hChargeRun->Reset();
    }

    if(entry%1000==0) printf(" ... %lld run %d  nhits %lu \n",entry,bEvent->run,bEvent->hits.size());
    
    //if(bEvent->run!=30440) continue;
    double evn = bEvent->npulse;

    // find start time 
    double qtrigSum=0;
    double triggerTime=0;
    double qmax=0;
    double qtime=0;
    double hitLast=0;
    TH1D *hTrigQEvent=NULL;
    TH1D *hTrigQSumEvent=NULL;
    if(entry<10) {
      hTrigQSumEvent = (TH1D*) hTrigEvent->Clone(Form("TrigQSumEvent%i",int(entry)));
      hTrigQEvent = (TH1D*) hTrigEvent->Clone(Form("TrigQEvent%i",int(entry)));
      hTrigQEvent->SetMarkerStyle(4);
      hTrigQSumEvent->SetMarkerStyle(4);
    }
    // loop to find trigger time 
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      double hitTime = bEvent->hits[ip].time*1E6;
      if(bEvent->hits.size()<10) continue;
      if(hitTime<0.9) continue;
      if(hitTime>1.5) continue;
      double qhit = 1.0E9*bEvent->hits[ip].q/qspe;  // spe
      qtrigSum+= qhit;
      int ibin = hTrigEvent->FindBin(hitTime);
      hTrigEventQhit->SetBinContent(ibin,qhit);
      hTrigEventQsum->SetBinContent(ibin,qtrigSum);
      if(hitLast>hitTime) break;
      if(hTrigQEvent) {
        hTrigQEvent->SetBinContent(ibin,qhit);
        hTrigQSumEvent->SetBinContent(ibin,qtrigSum);
        //printf("\t xxxx event %i hit %i ibin %i time %f q %f sum %f  \n",entry, ip, ibin, hitTime, qhit, qtrigSum);
      }
      if(hitTime>0.90&&hitTime<1.1) {
        if(qhit>qmax) {
          qmax=qhit;
          qtime = hitTime;
        }
        hTriggerQTime->Fill(hitTime,qtrigSum);
        hTriggerQ->Fill(hitTime,qhit);
      }
      // trigger time found as first hit in time window
      if(qhit>1&&hitTime>0.9&&hitTime<1.1&&triggerTime==0 ) triggerTime = bEvent->hits[ip].time*1E6; 
      hitLast=hitTime;
    }
    hTriggerTime->Fill(triggerTime);
    hTriggerTimeMax->Fill(qtime);
    nTrig->Fill(float(entry),triggerTime,qtime,qmax,qtrigSum);

    //printf(" ... xxxx trigger %lld run %d  trigtime %f %f  \n",entry,bEvent->run,triggerTime,qtime);
   
    //if(qtrigSum<20) continue;
    // loop to find max hit 
    double hitPre=0;
    double hitSum=0;
    double hitEarly=0;
    double hitLate=0;
    double hitMax=0;
    hLifePeak->Reset(0);
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      //cout << entry << "  " << hitq << endl;
      //if(theHit.order==0) triggerTime = theHit.time;
      double hitTime = 1.0E6*bEvent->hits[ip].time - triggerTime + 1;
      double hitq = 1.0E9*bEvent->hits[ip].q/qspe;  // spe
      hQSum->Fill(hitq);
      if(bEvent->hits.size()>10)  hQSum10->Fill(hitq);
      bool peakCut = (bEvent->hits[ip].peak>0.795&&bEvent->hits[ip].peak<0.804)||(bEvent->hits[ip].peak>1.604&&bEvent->hits[ip].peak<1.611);
      if(peakCut) continue; 
      bool afterCut = hitTime>afterLow&&hitTime<afterHigh;
      if(afterCut) continue;
      if(hitq>hitMax) hitMax=hitq;
      hitSum+= hitq;
      if(hitTime<triggerTime) hitPre += hitq;
      else if(hitTime>triggerTime&&hitTime<triggerTime+.5) hitEarly += hitq;
      else hitLate += hitq;
      int hitBin =  hLifePeak->FindBin(hitTime); 
      hLifePeak->SetBinContent(hitBin, hLifePeak->GetBinContent(hitBin)+hitq);
    } //
    int binMax = hLifePeak->GetMaximumBin();
    double hitPeak = hLifePeak->GetBinContent(binMax);
    double hitPeakTime = hLifePeak->GetBinLowEdge(binMax);
    ntEvent->Fill(float(entry),triggerTime,float(bEvent->hits.size()),hitMax,hitPeak,hitPeakTime,hitPre,hitSum,hitEarly,hitLate);
    oEvent->run=bEvent->run;  
    oEvent->event=bEvent->event;     
    oEvent->npulse=bEvent->npulse;    

    // min hits cut
    if(bEvent->hits.size()<10) continue;

    // throw out events with bad trigger
    if(triggerTime<0.9||triggerTime>1.1) continue;

    //run cut
    //if(bEvent->run>20090) continue;

    ++goodTrigger;

    // muonVmax cut 
    hMuVmax->Fill(bEvent->muVmax);
    if(VmaxCut&&bEvent->muVmax  < 0.05) continue;
    //if(!VmaxCut&&bEvent->muVmax > 0.05) continue;

    hRunSumNeil->Fill(bEvent->totQ/bEvent->spe);
    hEventSumNeil->Fill(bEvent->totQ/bEvent->spe);

    double QSumCut=0;
    // loop over pulses
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      //cout << entry << "  " << hitq << endl;
      //if(theHit.order==0) triggerTime = theHit.time;
      double hitTime = 1.0E6*bEvent->hits[ip].time - triggerTime + 1;
      double hitq = 1.0E9*bEvent->hits[ip].q/qspe;  // spe

      bool peakCut = (bEvent->hits[ip].peak>0.795&&bEvent->hits[ip].peak<0.804)||(bEvent->hits[ip].peak>1.604&&bEvent->hits[ip].peak<1.611);
      bool afterCut = hitTime>afterLow&&hitTime<afterHigh;
      if(hitTime>0.9&&!peakCut&&!afterCut) QSumCut+=hitq; 
      //QSumCut+=hitq; 

      // guess  at noise of tenth qspe
      //  ***** setting time unit 1.000000E-09 maxLife 10.000000 # digis 10000 
      double width = (qspe*.5)*1.0E9*bEvent->hits[ip].pwidth;
      //double hitqerr= sqrt(2*(abs(hitq)+ width*width));
      //gain is 4.36e6 electrons per photon
      double hitqerr= sqrt(abs(hitq)+ width*width);
      TPulse thePulse;
      thePulse.istart=bEvent->hits[ip].istart;
      thePulse.time= hitTime;
      thePulse.peak= bEvent->hits[ip].peak;
      thePulse.tpeak= bEvent->hits[ip].tpeak;
      thePulse.q = hitq;
      double hitres = hitq/hitqerr;
      thePulse.pwidth= width; 
      thePulse.qerr=hitqerr;
      oEvent->hits.push_back(thePulse);

      if(hitTime>7) hHitResidual->Fill(hitres);

      if(hitTime>timeCut) hQSpe->Fill(hitq);
      if(hitTime>timeCut) hQSpeCut->Fill(hitq);
      int hitBin =  hLife[0]->FindBin(hitTime); 
      hQTime->Fill(hitTime,hitq);
      double hitqerr2 = pow(hitqerr,2);
      if(isnan(hitqerr2)) printf (" XXXXXXX %f %f \n",hitq,hitqerr);
      // cut bad hits
      //if(peakCut||afterCut) continue;
      if(bEvent->hits.size()>10 ) {
        hLife[0]->SetBinContent(hitBin, hLife[0]->GetBinContent(hitBin)+hitq);
        hLife[0]->SetBinError(hitBin, sqrt( pow(hLife[0]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
        hLifeRun->SetBinContent(hitBin, hLifeRun->GetBinContent(hitBin)+hitq);
        hLifeRun->SetBinError(hitBin, sqrt( pow(hLifeRun->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      }
      if( hitSum>qCutLow&&hitSum<qCut) {
        hLife[1]->SetBinContent(hitBin, hLife[1]->GetBinContent(hitBin)+hitq);
        hLife[1]->SetBinError(hitBin, sqrt( pow(hLife[1]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      } 
      if(hitSum<qCutLow) {
        hLife[2]->SetBinContent(hitBin, hLife[2]->GetBinContent(hitBin)+hitq);
        hLife[2]->SetBinError(hitBin, sqrt( pow(hLife[2]->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      }
    }
    if(entry<1000) oTree->Fill();
    hEventSumCut->Fill(QSumCut);
    hRunSumCut->Fill(QSumCut);
    hChargeRun->Fill(QSumCut);
  }

  // run landau fit for last run
  hRunSumCut->Fit("landau");
  TF1 *lfit =  hRunSumCut->GetFunction("landau");
  hRunSumNeil->Fit("landau");
  TF1 *lfit2 =  hRunSumNeil->GetFunction("landau");
  if(lfit&&lfit2) ntRun->Fill(float(currentRun),float(hRunSumCut->GetEntries()),lfit->GetParameter(1),lfit->GetParError(1),lfit->GetParameter(2),lfit->GetParError(2),
          lfit2->GetParameter(1),lfit2->GetParError(1),lfit2->GetParameter(2),lfit2->GetParError(2),hLifeRun->Integral());

  // scale buy number of good triggers
  // for(int ih=0; ih<3; ++ih) hLife[ih]->Scale(1./goodTrigger);

  printf("  total events %lld goodTriggers %lld\n",nloop,goodTrigger);
  cout<< "  oTree size = " << oTree->GetEntries() << " total events " << hQSum10->GetEntries()  << endl;
  cout<< "  ntRun size = " << ntRun->GetEntries() << endl;
  //fout->ls();
  fout->Write();

  return;



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

  
  //gStyle->SetOptFit();
  TCanvas *canqsum10 = new TCanvas("qsumTen","qsumTen");
  hQSum10->Draw();


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
  Double_t binwidth = hLife[0]->GetBinWidth(1);
  Double_t integral = hLife[0]->Integral();

  double tauD=0.389;
  double tauX=0.02;
  double xppm=1.0;
  double tauA = 1.787;
  double tauM = 3.0;
  double Amystery = .001*integral;

  TF1* modelFit = new TF1("modelFit",fmodel,xstart,10,8);
  modelFit->SetNpx(1000); // numb points for function
  modelFit->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFit->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX, Amystery , tauM);
  modelFit->FixParameter(0,binwidth);
  modelFit->FixParameter(1,10);

  TF1* modelFitA = new TF1("modelFitA",fmodelA,xstart,10,8);
  modelFitA->SetNpx(1000); // numb points for function
  modelFitA->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX, Amystery , tauM);
  modelFitA->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX,1.0, 3 );
  modelFitA->FixParameter(0,binwidth);
  modelFitA->FixParameter(1,10);

  TF1* modelFitD = new TF1("modelFitD",fmodelD,xstart,10,8);
  modelFitD->SetNpx(1000); // numb points for function
  modelFitD->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX, Amystery , tauM);
  modelFitD->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX,1.0, 3 );
  modelFitD->FixParameter(0,binwidth);
  modelFitD->FixParameter(1,10);

  TF1* modelFitX = new TF1("modelFitX",fmodelX,xstart,10,8);
  modelFitX->SetNpx(1000); // numb points for function
  modelFitX->SetParameters(binwidth,xppm, 1.0 , tauA, tauD , tauX, Amystery , tauM);
  modelFitX->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
  modelFitX->FixParameter(0,binwidth);
  modelFitX->FixParameter(1,10);


  //modelFit->SetParameters(binwidth,10, integral ,1.787, .2 , 0.02, 1 , 3 );

  modelFit->SetLineColor(kBlack);
  //modelFit->GetHistogram()->GetYaxis()->SetRangeUser(1E-5,0.1);
  modelFit->GetHistogram()->GetXaxis()->SetRangeUser(0,4);
  modelFitA->SetLineColor(kRed);
  modelFitD->SetLineColor(kGreen);
  modelFitD->SetLineStyle(2);
  modelFitX->SetLineColor(kBlue);
  modelFit->SetLineColor(kBlack);


  modelFit->SetParameter(2,integral);
  modelFitA->SetParameter(2,integral);
  modelFitD->SetParameter(2,integral);
  modelFitX->SetParameter(2,integral);


  
  modelFit->Print();
  for(int ii=0; ii<8 ; ++ii) {
    printf(" param %i %s %.3f \n",ii,modelFit->GetParName(ii),modelFit->GetParameter(ii));
  }

  double A0 =  modelFitD->GetParameter(2);
  double ta = modelFitD->GetParameter(3);
  double tdp = modelFitD->GetParameter(4)/xppm;
  double tx = modelFitD->GetParameter(5);
  double tap = 1./( 1./ta + 1./tdp);
  double t1inv = ( -1./tdp + 1./tx);
  double t2inv = ( -1./tap + 1./tx);


  printf(" tapinv %f tdpinv %f txinv %f 1/t1=taudpinv-tauxinv %f 1/t2=taudpinv-tauxinv %f \n",1./tap,1./tdp,1./tx,t1inv,t2inv);

  /*
  for(unsigned j=0; j<100; ++j ) {
    double x = xstart+double(j)/10.;
    double xx = double(j)/10.;
    double DD =  A0*ta/tdp*(Exp(-xx/tdp) - Exp(-xx/tap)) ;
    printf(" x=%f A %.3E D %.3E (%.3E) X %.3E T %.3E \n",x, modelFitA->Eval(x),modelFitD->Eval(x),DD,modelFitX->Eval(x),modelFit->Eval(x));
  }
  */


  TCanvas *canFullFit = new TCanvas(Form("modelFit-tauD-%.1f",tauD),Form("modelFit-tauD-%.1f",tauD));
  canFullFit->SetLogy();
  modelFit->Draw("");
  //modelFit->Draw("same");
  modelFitA->Draw("same");
  modelFitD->Draw("same");
  modelFitX->Draw("same");


  
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
    integral = hLife[ih]->Integral();
    Amystery = .001*integral;

    xlow = xstart;
    xhigh= 10.0;

    printf(" fit range xlow %f bin %i xhigh %f bin %i %f \n",xlow,lowBin,xhigh,highBin,integral);
    
    fp[ih] = new TF1(Form("modelFit%i",ih),fmodel,xlow,xhigh,8);
    //fp[ih]->SetNpx(1000); // numb points for function
    fp[ih]->SetParNames("binwidth","xppm","A0","tauA","tauD","tauX","Am","tauM");
    fp[ih]->FixParameter(0,binwidth);
    fp[ih]->SetParameters(binwidth,10, integral ,tauA, tauD , tauX, Amystery , tauM );
    //fp[ih]->FixParameter(3,tauA);
    fp[ih]->FixParameter(5,tauX);
    fp[ih]->SetParLimits(6,0,0.1*integral);
    fp[ih]->SetParLimits(7,1,10);
    //if(ih==2) fp[ih]->FixParameter(7,10000);

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


  for(int ibin =1; ibin< hLife[0]->GetNbinsX() -1 ; ++ibin) {
    double x =  hLife[0]->GetBinLowEdge(ibin);
    if(x<xstart) continue;
    double fitval = fp[0]->Eval(x);
    double residual = 0;
    if(hLife[0]->GetBinError(ibin)>0) (hLife[0]->GetBinContent(ibin)-fitval)/hLife[0]->GetBinError(ibin);
    nResidual->Fill(float(ibin),x,residual,fitval,hLife[0]->GetBinContent(ibin), hLife[0]->GetBinError(ibin));
    hLifeResiduals->Fill(residual);
  }



  TCanvas *can[3];
  gStyle->SetOptFit();

  hLife[0]->GetYaxis()->SetRangeUser(0.1, hLife[0]->GetBinContent( hLife[0]->GetMaximumBin())*1.1);
  hLife[1]->GetYaxis()->SetRangeUser(0.1,hLife[1]->GetBinContent( hLife[1]->GetMaximumBin())*1.1);
  hLife[2]->GetYaxis()->SetRangeUser(0.1,hLife[2]->GetBinContent( hLife[2]->GetMaximumBin())*1.1);


  can[0] = new TCanvas(Form("life-%s-pass-%i-hits",tag.Data(),10),Form("life-%s-passs-%i-hits",tag.Data(),10));
  hLife[0]->SetTitle( Form("life-%s-pass-%i-hits",tag.Data(),10) );
  can[0]->SetLogy();
  hLife[0]->Draw("E1");
  can[0]->Print(".png");


  can[1] = new TCanvas(Form("life-%s-Q-%.0f-%.0f",tag.Data(),qCutLow,qCut),Form("life-%s-%.0f-%0.f",tag.Data(),qCutLow,qCut));
  can[1]->SetLogy();
  hLife[1]->SetTitle( Form("life-%s-Q-%.0f-%.0f",tag.Data(),qCutLow,qCut) );
  hLife[1]->Draw("E1");


  can[2] = new
  TCanvas(Form("life-%s-Q-%.0f",tag.Data(),qCutLow),Form("life-%s-Q-%.0f",tag.Data(),qCutLow));
  hLife[2]->SetTitle( Form("life-%s-Q-%.0f",tag.Data(),qCutLow) );
  can[2]->SetLogy();
  hLife[2]->Draw("E1");

}
