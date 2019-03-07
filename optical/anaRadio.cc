#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
#include "MGTMCEventHeader.hh"
#include "io/MGOutputMCOpticalRun.hh"
#include "MGTMCRun.hh"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TString.h"
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>
#include <TDatime.h>
#include <TNtuple.h>
#include <unordered_map>
#include <map>

using namespace std;


int main(int argc, char *argv[])
{ 
  if(argc<2) {
    std::cout << "Usage: " << argv[0] << " <input file name without .root extension in directory ../ >" << '\n';
    return 0;
  }
  TString fileName = argv[1];

  // setup input 
  TString dir = "../";
  TString fullFileName = dir + fileName + TString(".root");
  TFile* infile = TFile::Open(fullFileName,"READONLY");
  if(!infile){
    cout << " cannot find file " << fullFileName << endl;
    return 0;
  }

  //TFile *infile = new TFile(dir+fileName,"READONLY");
  TChain *fTree = new TChain("fTree");
  printf(" adding to tree %s \n",fullFileName.Data());
  Int_t nfiles = fTree->Add(fullFileName);
  //MGOutputMCOpticalRun* fMCOpticalRun;	
  fTree->GetListOfBranches()->ls();
  //fTree->SetBranchAddress("fMCOpticalRun",&fMCOpticalRun);

  Long64_t nentries = fTree->GetEntries();
  cout<<"TTree with  "<< nfiles <<" file and  "<<nentries<<" entries"<<endl;
  fTree->GetListOfFiles()->ls();

  MGTMCEventSteps*  eventSteps=NULL;      // MGDO encapsulation of steps 
  MGTMCEventSteps*  eventPrimaries=NULL;  // MGDO encapsulation of steps 
  MGTMCEventHeader* eventHeader=NULL;
  fTree->SetBranchStatus("fMCOpticalRun",false);
  fTree->SetBranchAddress("eventHeader",&eventHeader);
  fTree->SetBranchAddress("eventSteps",&eventSteps);
  fTree->SetBranchAddress("eventPrimaries",&eventPrimaries);

  // setup outpuf file 
  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  TString outFileName = TString("Radio.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(outFileName,"recreate");
  cout << " opened output file " << outFileMaGe->GetName() << endl;
  outFileMaGe->cd();
  TH1D* hEThresh = new TH1D("EThresh","EThresh",1000,0,1.);
  TNtuple *ntStepLAr = new TNtuple("stepLAR","LAr Core","event:nstep:edep:ke:x:y:z:r:t");
  TNtuple *ntLAr = new TNtuple("primaryLAr","primary LAr","x:y:z:r:theta:px:py:pz:nPhotons:edep:mapThresh:pmtCounter:totalPhotons");
  TNtuple *ntGe = new TNtuple("primaryGe","primary Ge","x:y:z:r:theta:px:py:pz:nPhotons:edep:mapThresh:pmtCounter:totalPhotons");

  TString mapFileName = "OpticalMapBACONArgon.1e10.5mm";
  TFile* mapFile = TFile::Open(mapFileName+TString(".root"),"READONLY");
  if(!mapFile) {
    cout << " failed to open map file " << mapFileName << endl;
    return 0;
  }
  //mapFile->ls();
  TH3D* hMap = NULL;
  mapFile->GetObject("EnergyMap",hMap);
  TH2D *hYZ = NULL;
  mapFile->GetObject("2DOpticalMap_YZ",hYZ);
  if(!hMap) { printf("!!!!! missing map hMap \n") ;return 0;}
  if(!hYZ ) { printf("!!!!! missing map hYZ \n") ;return 0;}
  printf("...... maps %s %s \n",hMap->GetTitle(),hYZ->GetTitle());

  printf("starting run over events %lld \n",nentries);

  Int_t noHitcounter = 0;
  for(Int_t i = 0; i < nentries ; i++){
    if((i+1)%100 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
    fTree->GetEntry(i);
    TString physName;
    int nsteps = eventPrimaries->GetNSteps();
    printf(" event %i NSteps %i \n",i,nsteps);    
    if(nsteps<1) { 
      printf("  event %i NSteps %i MISSING PRIMARY! \n",i,nsteps); 
      continue;
    }
    const MGTMCStepData* primaries = eventPrimaries->GetStep(0);
    if(primaries == NULL){
      cout<<"null primary" << " nsteps " << nsteps <<endl;
      continue;
    }
    //virtual inline TVector3 GetMomentumVector() const { return TVector3(fPx, fPy, fPz); }
    //virtual inline TVector3 GetPositionVector() const { return TVector3(fX, fY, fZ); }
    //virtual inline TVector3 GetLocalPositionVector() const { return TVector3(fLocalX, fLocalY, fLocalZ); }
    
    Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
    Double_t px = primaries->GetPx(),py = primaries->GetPy(),pz = primaries->GetPz();
    Double_t r = sqrt(x*x+y*y);
    Double_t theta = std::acos(x/r);
    if( y < 0) theta += TMath::Pi(); 
    Double_t nPhotons = 0,eDepLAr = 0,eDepGe = 0,mapThresh = 0.;
    Int_t stepCounter = 0,pmtCounter = 0,nTotalPhotons = 0;
    bool hitLAr = false,hitGe = false;
    for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
      const MGTMCStepData *step = eventSteps->GetStep(j);
      physName = step->GetPhysVolName();
      TString procName = step->GetProcessName();
      if(step->GetTrackWeight() == 1 && step->GetStepNumber() == 1) {
        nTotalPhotons++;
      }
      //eDep +=step->GetEdep();
      if(physName.Contains("physicalPMT")&& procName.Contains("WLS") && step->GetParticleID() == 0 && step->GetEdep() != 0) pmtCounter++;
      //if(physName.Contains("SiPM")&& step->GetParticleID() == 0 && step->GetEdep() != 0) pmtCounter++;
      if(physName =="Detector"  ){
        if(step->GetParticleID() == 0) continue;
        Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
        Double_t eThresh = (hMap->GetBinContent(bin)/1000.); //map in keV, Geant4 is in MeV
        if(eThresh == 0){
          eThresh = 100.; //100 MeV as a large threshold
        }
        hEThresh->Fill(eThresh);
        nPhotons += step->GetEdep()/eThresh;
        mapThresh+= eThresh;
        stepCounter++;
        eDepLAr +=step->GetEdep();
        hitLAr = true;
        ntStepLAr->Fill(i,j,
            step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),
            sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT());
      }
      else if(physName.Contains("ActiveDet")){
        eDepGe += step->GetEdep();
        hitGe = true;
      }
    }
    if(hitLAr){
      ntLAr->Fill(x,y,z,r,theta,px,py,pz,nPhotons,eDepLAr,mapThresh/stepCounter,pmtCounter,nTotalPhotons);
    }
    if(hitGe){ 
      ntGe->Fill(x,y,z,r,theta,px,py,pz,nPhotons,eDepGe,mapThresh/stepCounter,pmtCounter,nTotalPhotons);
    }
    if(!hitGe && !hitLAr){
      noHitcounter++;
    }
  }
  delete hMap;
  delete hYZ;
  outFileMaGe->cd();
  outFileMaGe->ls();
  outFileMaGe->Write();
  //printf(" closing output file \n");
  outFileMaGe->Close();
  //printf(" closing infile \n");
  infile->Close();
  //printf(" closing mapfile \n");
  mapFile->Close();
  printf("Total events = %lld LAR %i GE %i  no hit %i \n ",nentries,int(ntLAr->GetEntriesFast()),int(ntGe->GetEntriesFast()),noHitcounter);
  cout<<"end of file  "<<fileName<<endl;
  return 0;
}
