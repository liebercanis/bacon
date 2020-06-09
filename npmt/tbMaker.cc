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
#include "TBaconRun.hxx"

class tbMaker {
  public :
    tbMaker(Int_t runStart=3000, Int_t runStop=3001);
    virtual ~tbMaker(){;}
    Long64_t nloop;
    TTree *tbTree;
    TTree *tPulse;
    TBaconEvent *bEvent;
    //TBaconRun *bRun;  
 };



tbMaker::tbMaker(Int_t runStart=3000, Int_t runStop)
{
  

  TFile *fout = new TFile(Form("tbMaker-%i-%i.root",runStart,runStop),"RECREATE");
  cout << " output file is  " << fout->GetName() << endl;
  fout->cd();

  double meanSPE  = 4.60e-11;
  double sigmaSPE = 1.5e-11;

  TString fileName;
  TString fileDir = TString("processedData/DS1/");

   
  /* ntuplePulse structure */
  Float_t         ientry;
  Float_t         pmtNum;
  Float_t         nhits;
  Float_t         charge;
  Float_t         startTime;
  Float_t         peakWidth;
  Float_t         T0;
  Float_t         vMax;
  Float_t         vMaxTime;
  Float_t         Sdev;
  Float_t         baseline;
  tPulse = new TTree("ntuplePulse","ntuplePulse");
  tPulse->Branch("ientry", &ientry);
  tPulse->Branch("pmtNum", &pmtNum);
  tPulse->Branch("nhits", &nhits);
  tPulse->Branch("charge", &charge);
  tPulse->Branch("startTime", &startTime);
  tPulse->Branch("peakWidth", &peakWidth);
  tPulse->Branch("T0", &T0);
  tPulse->Branch("vMax", &vMax);
  tPulse->Branch("vMaxTime", &vMaxTime);
  tPulse->Branch("Sdev", &Sdev);
  tPulse->Branch("baseline", &baseline);

  tPulse->GetListOfBranches()->ls();



  TNtuple *nt;
  vector<int> runEvents;
  vector<int> vrun;
  TList* ntList = new TList;
  for(Int_t irun=runStart ; irun<=runStop; ++irun) {
    fileName = fileDir+TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(irun)+TString(".root");
    TFile *fin  = new TFile(fileName,"READONLY");
    fin->GetObject("ntuplePulse",nt);
    if(nt) {
      ntList->Add(nt);
      runEvents.push_back(int(nt->GetEntries()));
      vrun.push_back(irun);
    }
  }


  vector<int> runEventSum(vrun.size());
  int runsum=0;
  for(unsigned ir=0; ir<vrun.size(); ++ir) {
    runsum += runEvents[ir];
    runEventSum[ir]=runsum;
    printf(" run %i events %i sum %i \n",vrun[ir],runEvents[ir],runEventSum[ir]);
  }

  tPulse->Merge(ntList);
 

    
  cout << " ntpulse has " << tPulse->GetEntries() << endl;

  tbTree = new TTree("TBacon","TBacon data");
  bEvent = new TBaconEvent();
  tbTree->Branch("bevent",&bEvent);
  
  cout << " TBacon has " << tbTree->GetEntries() << endl;


  // fill TBacon 
  int irun=0;
  for(Long64_t entry=0; entry< tPulse->GetEntries() ; ++entry) {

    tPulse->GetEntry(entry);
    if(entry>=runEventSum[irun]) ++irun;

    // ientry is event index //

    //cout << "run  " << vrun[irun] <<  " entry " << entry << "  charge " << charge*1E9 << endl;

    /**

        baconEvent->event=ientry;
        baconEvent->run=z;
        baconEvent->npulse=nPulses;
        baconEvent->npmt=j;
        baconEvent->totQ = pmtRun->totalCharge[j];
        baconEvent->spe = meanSPE[j];
        baconEvent->muVmax=muonVmax;


        baconEvent->T0=pmtRun->T0[j];
        baconEvent->totalCharge=pmtRun->totalCharge[j];
        baconEvent->tMax=       pmtRun->tMax[j];
        baconEvent->vMax=       pmtRun->vMax[j];
        baconEvent->cMax=       pmtRun->cMax[j];
        baconEvent->baseline=   pmtRun->baseline[j];
        baconEvent->sDev=       pmtRun->sDev[j]; 


        if(ientry%1000 == 0 || ientry == NEvents - 1) printf(" ev %i  pmt %i totQ %E \n ",ientry, j,baconEvent->totQ);

        for(int i = 0; i < nPulses;i++){
          if(debug && i == 0) cout<<"Looping over Pulses"<<endl;
          Double_t charge     = pmtRun->charge[j][i];
          Double_t startTime  = pmtRun->startTimes[j][i];
          Double_t peakWidth  = pmtRun->peakWidths[j][i];
          Double_t peakHeight = pmtRun->peakHeights[j][i];
          Double_t peakTime   = pmtRun->peakTimes[j][i];
          TPulse thePulse;
          thePulse.time=startTime;
          thePulse.tpeak=peakTime;
          thePulse.pwidth=peakWidth;
          thePulse.peak=peakHeight;
          thePulse.q=charge;
          //thePulse.qerr=phitQErr;


          //SPE Fill
          if(peakTime > 7e-6) hSPE[j]->Fill(charge);

          //cut on small peaks
          debug = false;
          if(peakHeight < vMinCut && j == 0){
            if(debug) cout<<"vMinCut "<<vMinCut<<", pulse value "<<peakHeight<<", PMT "<<j<<endl;
            continue;
          }
          if(peakWidth < peakWidthCut && j == 0) {
            if(debug) cout<<"peakWidthCut "<<peakWidthCut<<", pulse value "<<peakWidth<<", PMT "<<j<<endl;
            continue;
          }

          //avoid events that start early
          if(startTime < tMinCut && j == 0){
            if(debug) cout<<"tMinCut "<<tMinCut<<", pulse value "<<startTime<<", PMT "<<j<<endl;
            continue;
          }
          //avoid events with large maximum
          //i.e. events that go out of range
          if(vMax >= vMaxEventCut && j == 0){
            if(debug) cout<<"vMaxEventCut "<<vMaxEventCut<<", pulse value "<<vMax<<", PMT "<<j<<endl;
            continue;
          }
          //Fill Triplet
          //if(peakTime > 6e-6 && peakTime < 7e-6) cout<<"charge "<<charge<<endl;
          Int_t peakTimeBin = hTripletChargeWeightedSummed->FindBin(peakTime);
          Int_t nPhotons = charge/meanSPE[j];
          //triplet fit for both channels
          hTripletChargeWeightedSummed->Fill(peakTime,charge/meanSPE[j]);
          //hTripletChargeWeightedSummed->Fill(peakTime,charge);
          //hTripletChargeWeightedPMT[j]->Fill(peakTime,charge);
          if(j == 0){
            hTripletChargeWeightedPMT[0]->Fill(peakTime,charge/meanSPE[j]);
            hTripletChargeWeightedPMT[1]->Fill(peakTime,charge/meanSPE[j]);
          }
          baconEvent->hits.push_back(thePulse);
        }
        //if(debug) cout<<"F40 cut "<<endl;
        //pulse finding loop
        //F40 cut time window
        Double_t F40 = 0;
        Double_t F40Window = 40e-9;
        Double_t F40Start = 1e-6,F40Stop = F40Start+ F40Window;
        for(int i = F40Start/deltaT; i <= F40Stop/deltaT; i++){
          F40 += -volts[i]*deltaT;  
        }
        //if(debug) cout<<" Filling final histogram"<<endl;
        hF40Summed->Fill(F40/totalCharge);
        baconEvent->F40=F40;
        TBacon->Fill();
      }

**/

   }


  fout->ls();
  fout->Write();

  return;
}
