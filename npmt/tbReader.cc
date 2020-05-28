/* uses TBaconEvent Class */
//////////////////////////////////////////////////////////
//  M.Gold May 2020 
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <complex>//includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort
#include "TSpectrum.h"
#include "TRandom3.h"

#include "TBaconEvent.hxx"

class tbReader {
  public :
    tbReader(Int_t runstart=20000, Int_t runstop=20009);
    virtual ~tbReader(){;}
    void newRun();
    Int_t currentRun;
    Long64_t treeNumber;
    Long64_t goodTrigger;
    Long64_t totTrigger;
    Long64_t nloop;
    TH1D *hLife;
    TH1D *hChargeSum;

    TH1D *hLifeRun;
    TH1D *hChargeSumRun;
    TH1D *hChargeCutRun;
    TH1D *hLifeQCutRun;
    TH1D *hChargeQCutRun;
    
    TChain *tree;
    TBaconEvent *bEvent;
    TNtuple *ntRun;
    TNtuple *ntTrig;
};


void tbReader::newRun() {
  printf(" xxxxxxxxx at entry %lld end of run %i tot(good) %lld(%lld) triggers\n" , nloop, currentRun, totTrigger, goodTrigger) ;

  if(hChargeCutRun) { 
    hChargeCutRun->Fit("landau");
    TF1 *lfit =  hChargeCutRun->GetFunction("landau");
    hChargeSumRun->Fit("landau");
    TF1 *lfit2 =  hChargeSumRun->GetFunction("landau");
    hChargeQCutRun->Fit("landau");
    TF1 *lfit3 = hChargeQCutRun->GetFunction("landau");

    if(lfit&&lfit2&&lfit3) ntRun->Fill(float(currentRun),float(hChargeSumRun->GetEntries()),
        lfit->GetParameter(1), lfit->GetParError(1),  lfit->GetParameter(2),lfit->GetParError(2),
        lfit2->GetParameter(1),lfit2->GetParError(1),lfit2->GetParameter(2),lfit2->GetParError(2),
        lfit3->GetParameter(1),lfit3->GetParError(1),lfit3->GetParameter(2),lfit3->GetParError(2));

    printf(" xxxxxxxxx ntRun size %lld events %E totSPE %E \n" , ntRun->GetEntries(),hChargeSumRun->GetEntries(),hLifeRun->Integral());

  }

  totTrigger=0;
  goodTrigger=0;
  currentRun = bEvent->run;
  hLifeRun = (TH1D*) hLife->Clone(Form("LifeRun%5i",int(currentRun)));
  hChargeSumRun = (TH1D*) hChargeSum->Clone(Form("ChargeSumRun%5i",int(currentRun)));
  hChargeCutRun = (TH1D*) hChargeSum->Clone(Form("ChargeCutRun%5i",int(currentRun)));
  hLifeQCutRun = (TH1D*) hLife->Clone(Form("LifeQCutRun%5i",int(currentRun)));
  hChargeQCutRun = (TH1D*) hChargeSum->Clone(Form("ChargeQCutRun%5i",int(currentRun)));
}


tbReader::tbReader(Int_t runstart, Int_t runstop)
{
  currentRun=-1;
  treeNumber-1;
  goodTrigger=0;
  totTrigger=0;
  hLifeRun=NULL;
  hChargeSumRun=NULL;
  hChargeCutRun=NULL;
  hLifeQCutRun=NULL;
  hChargeQCutRun=NULL;

  double maxLife=10.0;
  int lifeBins = int(1.0E4/8.0); // ns bins
  double QSumCutValue = 200;
  //int lifeBins=400;

  //bool VmaxCut = false;
  
  // files must be in order by run
  TString nfile[6];
  nfile[0]=TString("TBAcon_10000_20000_DS_2.root");
  nfile[1]=TString("TBAcon_20001_20100_DS_2.root");
  nfile[2]=TString("TBAcon_20101_20222_DS_2.root");
  nfile[3]=TString("TBAcon_30000_30452_DS_2.root");
  nfile[4]=TString("TBAcon_30440_30443_DS_2.root");
  nfile[5]=TString("TBAcon_30453_40000_DS_2.root");
  tree = new TChain("TBacon");
  bEvent = new TBaconEvent();
  tree->SetBranchAddress("bevent",&bEvent);
  tree->GetListOfFiles()->ls();
  for(int ifile=0; ifile<6; ++ifile) tree->Add(nfile[ifile]);
  cout << " TBacon has " << tree->GetEntries() << endl;

  TFile *fout = new TFile(Form("tbReader-%i-%i-%i-bins.root",runstart,runstop,lifeBins),"RECREATE");
  cout << " output file is  " << fout->GetName() << endl;
  TDirectory *trigDir = fout->mkdir("trigDir");
  fout->cd();
  ntRun = new TNtuple("ntRun","","run:nev:mpvc:mpvcErr:widc:widcErr:mpvn:mpvnErr:widn:widnErr:mpvq:mpvqErr:widq:widqErr");

  TH2D *hNHitTotQ = new TH2D("NHitByTotQ"," nhit by tot Q ",50,0,50,50,0,50);
  TH1D *hNHits = new TH1D("NHits"," hits in event  ",350,0,350);
  TH1D *hTotQ  = new TH1D("TotQ"," totQ in event  after 10 cut ",350,0,350);


  // trigger
  std::map<double,double> mapDer;
  std::map<double,double> mapMax;
  ntTrig =  new TNtuple("ntTrig","","event:sigTime:sigQ:derTime:derQ:maxTime:maxQ");
  

  TH1D *hTrigQhit = new TH1D("TrigQhit"," trig hit charge ",1000,.9,1.1);
  TH1D *hTrigQsum = new TH1D("TrigQsum"," trig sum charge ",1000,.9,1.1);
  TH2D *hTriqQsumvsTime = new TH2D("TrigQsumvsTime"," charge vs time  ",100,.9,1.1,1000,0,200);

  TH1D *hTrigTimeQsig = new TH1D("TrigTimeQsig"," 3 Qspe Trigger time ",240,.9,1.5);
  TH1D *hTrigTimeQder = new TH1D("TrigTimeQder"," max der Trigger time ",240,.9,1.5);
  TH1D *hTrigTimeQmax = new TH1D("TrigTimeQmax"," max Q   Trigger time ",240,.9,1.5);

  
  // event 
  TH1D *hMuVmax = new TH1D("MuVmax"," muVmax ",1000,0,1);
  hMuVmax->GetXaxis()->SetTitle("muon Vmax" );
  TNtuple *ntEvent = new TNtuple("ntEvent"," event ","entry:run:ntrig:nhits:muVmax:QSumCut");

  //hit
  TH1D *hHitResidual = new TH1D("HitResidual7microsec"," residual  ",1000,0,5);

  
  // prototypes 
  hLife = new TH1D("LifeCut"," lifetime PMT >10 hits ",lifeBins,0,maxLife);
  hLife->GetXaxis()->SetTitle(" micro-seconds ");
  hLife->SetMarkerColor(kBlack);
  hLife->SetMarkerStyle(22);
  hLife->SetMarkerSize(.2);

  hChargeSum = new TH1D("ChargeSum"," good summed event charge > 10 pulses  ",1000,0,10000);
  hChargeSum->GetXaxis()->SetTitle(" run sum dt-Q (x10^9) ");

   
  nloop = tree->GetEntries();
  for(Long64_t entry =0; entry< nloop ; ++ entry) {
    Long64_t itree = tree->LoadTree(entry);
    if(itree!=treeNumber) {
      tree->SetBranchAddress("bevent",&bEvent);
      treeNumber= itree;
    }
    tree->GetEntry(entry);
    // find run range
    if(bEvent->run < runstart) continue;
    if(bEvent->run > runstop)  break;
    
    // if new run
    if(bEvent->run!=currentRun) newRun();     

    if(bEvent->npmt!=0) continue; // only use PMT zero

    ++totTrigger;
    if(entry%1000==0) printf(" ... %lld run %d  nhits %lu \n",entry,bEvent->run,bEvent->hits.size());

    hNHitTotQ->Fill(bEvent->totQ,bEvent->hits.size());
    hNHits->Fill(bEvent->hits.size());
    if(bEvent->hits.size()<10) continue;
    hTotQ->Fill(bEvent->totQ);

    double qspe = 1E9*bEvent->spe; 
    // early trigger data
    TH1D *hTrigQhitEvent =  NULL;
    TH1D *hTrigQsumEvent = NULL;
    if(totTrigger<10) {
      trigDir->cd();
      hTrigQhitEvent =  (TH1D*) hTrigQhit->Clone(Form("TrigQRun%iEvent%i",int(bEvent->run),int(entry)));
      hTrigQsumEvent =  (TH1D*) hTrigQsum->Clone(Form("TrigQSumRun%iEvent%i",int(bEvent->run),int(entry)));
      fout->cd();
    } 
    
    // find start time
    double lastQhit=0;
    double qtrigSum=0;
    double triggerTime=-1;
    double trigQ=0;
    double trigSig=3;
    mapDer.clear();
    mapMax.clear();
    // hit loop to find trigger time 
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      double hitTime = bEvent->hits[ip].time*1E6;
      if(hitTime<0.9) continue;
      if(hitTime>1.5) continue;
      double hitq = 1.0E9*bEvent->hits[ip].q/qspe;  // spe
      mapDer.insert( pair<double,double>(hitq-lastQhit,hitTime) );
      lastQhit=hitq;
      mapMax.insert( pair<double,double>(hitq,hitTime) );
      qtrigSum+= hitq;
      int ibin = hTrigQhit->FindBin(hitTime);
      hTrigQhit->SetBinContent(ibin,hitq);
      hTrigQsum->SetBinContent(ibin,qtrigSum);
      hTriqQsumvsTime->Fill(hitTime,qtrigSum);

      if(hitq>trigSig&&triggerTime==-1) { 
        triggerTime=hitTime;
        trigQ=hitq;
      }

      if(hTrigQhitEvent) {
        hTrigQhitEvent->SetBinContent(ibin,hitq);
        hTrigQsumEvent->SetBinContent(ibin,qtrigSum);
        //printf("\t xxxx event %i hit %i ibin %i time %f q %f sum %f  \n",entry, ip, ibin, hitTime, hitq, qtrigSum);
      }
    }
    // hit loop to find QSum 
    double QSumCut=0;
    double QSum=0;
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      double hitTime = bEvent->hits[ip].time*1E6;
      double hitq = 1.0E9*bEvent->hits[ip].q/qspe;  // spe
      if(hitq<=0) continue;
      bool peakCut = (bEvent->hits[ip].peak>0.795&&bEvent->hits[ip].peak<0.804)||(bEvent->hits[ip].peak>1.604&&bEvent->hits[ip].peak<1.611);
      QSum += hitq;
      if(!peakCut) QSumCut += hitq;
    }

    hChargeSumRun->Fill(QSum);
    hChargeCutRun->Fill(QSumCut);
    if(QSumCut>QSumCutValue) hChargeQCutRun->Fill(QSumCut);


    std::map<double,double>::iterator derIter;
    std::map<double,double>::iterator maxIter;

    derIter =  mapDer.begin();
    maxIter =  mapMax.begin();


    hTrigTimeQsig->Fill(triggerTime);
    ntTrig->Fill(float(entry),triggerTime,trigQ,derIter->second,derIter->first,maxIter->second,maxIter->first);

    // throw out events with bad trigger
    if(triggerTime<0.9||triggerTime>1.1) continue;

    ++goodTrigger;

    // muonVmax cut 
    hMuVmax->Fill(bEvent->muVmax);
    //if(VmaxCut&&bEvent->muVmax  < 0.05) continue;

       // loop over pulses
    for(unsigned ip=0; ip< bEvent->hits.size(); ++ip) {
      //cout << entry << "  " << hitq << endl;
      //if(theHit.order==0) triggerTime = theHit.time;
      double hitTime = 1.0E6*bEvent->hits[ip].time - triggerTime + 1;
      double hitq = 1.0E9*bEvent->hits[ip].q/qspe;  // spe

      // throw out negative hits
      if(hitq<=0) continue;

      // guess  at noise of tenth qspe
      //  ***** setting time unit 1.000000E-09 maxLife 10.000000 # digis 10000 
      double width = (qspe*.5)*1.0E9*bEvent->hits[ip].pwidth;
      //double hitqerr= sqrt(2*(abs(hitq)+ width*width));
      //gain is 4.36e6 electrons per photon
      double hitqerr= sqrt(abs(hitq)+ width*width);
      double hitres = hitq/hitqerr;

      if(hitTime>7) hHitResidual->Fill(hitres);

      int hitBin =  hLifeRun->FindBin(hitTime); 
      double hitqerr2 = pow(hitqerr,2);
      if(isnan(hitqerr2)) printf (" XXXXXXX %f %f \n",hitq,hitqerr);
      // cut bad hits
      //if(peakCut||afterCut) continue;
      hLifeRun->SetBinContent(hitBin, hLifeRun->GetBinContent(hitBin)+hitq);
      hLifeRun->SetBinError(hitBin, sqrt( pow(hLifeRun->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      hLifeRun->SetBinContent(hitBin, hLifeRun->GetBinContent(hitBin)+hitq);
      hLifeRun->SetBinError(hitBin, sqrt( pow(hLifeRun->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      if(QSumCut>QSumCutValue) {
        hLifeQCutRun->SetBinContent(hitBin, hLifeQCutRun->GetBinContent(hitBin)+hitq);
        hLifeQCutRun->SetBinError(hitBin, sqrt( pow(hLifeQCutRun->GetBinError(hitBin),2)+pow(hitqerr,2) ));
        hLifeQCutRun->SetBinContent(hitBin, hLifeQCutRun->GetBinContent(hitBin)+hitq);
        hLifeQCutRun->SetBinError(hitBin, sqrt( pow(hLifeQCutRun->GetBinError(hitBin),2)+pow(hitqerr,2) ));
      }
    }

    ntEvent->Fill(float(entry),float(bEvent->run),float(goodTrigger),float(bEvent->hits.size()),float(bEvent->muVmax),float(QSumCut));
  }
  // run landau fit for last run

  //fout->ls();
  fout->Write();

  return;
}