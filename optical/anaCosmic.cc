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

inline bool fileExist (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main(int argc, char *argv[])
{ 
  std::vector<TString> knownName;
  knownName.push_back("lab");
  knownName.push_back("cryostat");
  knownName.push_back("airSpace");
  knownName.push_back("argonGasPhysical");
  knownName.push_back("Detector");
  knownName.push_back("physicalPMT_0");
  knownName.push_back("physicalPMT_1");
  knownName.push_back("physicalHousingPMT_0");
  knownName.push_back("physicalHousingPMT_1");
  knownName.push_back("physicalWLS_0");
  knownName.push_back("physicalWLS_1");
  knownName.push_back("World");
  printf(" known physical names %lu \n",knownName.size());
  for(unsigned iname=0; iname<knownName.size(); ++ iname) printf("%u %s \n",iname,knownName[iname].Data());

  if(argc<2) {
    std::cout << "Usage: " << argv[0] << " <input file name without .root extension in directory ../ >" << '\n';
    return 0;
  }
  int nprocess=0;
  if(argc>2) nprocess = atoi(argv[2]);
  TString fileName = argv[1];
  TString dir = "/home/gold/XenonDoping/";
  TString fullFileName= dir+fileName+TString(".root");
  TFile* infile = TFile::Open(fullFileName,"readonly");
  if(!infile){
    printf(" cannot find input file %s\n",fullFileName.Data()); 
    return 0;
  }
  printf(" reading input file %s\n",fullFileName.Data()); 

  TChain *fTree = new TChain("fTree");
  fTree->Add(dir+fileName+TString(".root"));
  Long64_t nentries = (Long64_t)fTree->GetEntries();

  MGTMCEventSteps *eventSteps = 0;
  MGTMCEventSteps *eventPrimaries = 0;
  if(fTree != NULL){
    fTree->SetBranchAddress("eventSteps",&eventSteps);
    fTree->SetBranchAddress("eventPrimaries",&eventPrimaries);
  }
  else{
    cout<<"NULL fTree"<<endl;
    return 0;
  }

  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  TString outFileName = TString("anaCosmic-")+fileName+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opened output %s \n",outFile->GetName());
  TNtuple *ntPrim = new TNtuple("ntPrim","","nsteps:e:z:r:theta:edep:ltot:eloss:zlast:rlast");
  TNtuple *ntStepMu = new TNtuple("ntStepMu","","nstep:pid:vol:edep:length:ltot:eloss:zstep:rstep");
  TNtuple *ntStepOpt = new TNtuple("ntStepOpt","",  "nstep:pid:hit:edep:length:ltot:eloss:zstep:rstep");
  TNtuple *ntStepOth = new TNtuple("ntStepOth","","nstep:pid:vol:edep:length:ltot:eloss:zstep:rstep");


  if(nprocess==0) nprocess=nentries;
  cout<<">> Processing "<<nprocess<< " of " <<nentries<<" entries"<<endl;
  for(Long64_t ientry = 0; ientry < nprocess ; ++ientry) {
    //if(ientry%10 == 0 ) cout<<"\tprocessed "<< ientry <<" events"<<endl;
    fTree->GetEntry(ientry);
    const MGTMCStepData *prim0 = eventPrimaries->GetStep(0);
    if(prim0 == NULL){
      cout<<"null prim0"<<endl;
      continue;
    }
    TVector3 pvec = prim0->GetMomentumVector(); 
    TVector3 rvec = prim0->GetPositionVector(); 
    Double_t r = rvec.Perp();
    Double_t theta = pvec.Theta()*180.0/TMath::Pi();
    Double_t energy = prim0->GetKineticE();
    if(ientry%10 == 0 ) { 
      cout<<"processed "<< ientry << " pid " << prim0->GetParticleID() << " primary steps " << eventPrimaries->GetNSteps() 
        << " length " << prim0->GetStepLength() << " energy " << energy << " theta " << theta  
        << " event steps  " << eventSteps->GetNSteps() <<  endl;
    }
    TVector3 rstep;
    TVector3 rMuStep;
    Double_t eloss=0;
    Double_t edep=0;
    Double_t totLength=0;
    Double_t edepMu=0;
    Double_t totLengthMu=0;
    Double_t elossMu=0;
    for(int jstep=0; jstep<eventSteps->GetNSteps(); ++ jstep) {
      const MGTMCStepData *evstep = eventSteps->GetStep(jstep);
      int pId = evstep->GetParticleID(); 
      rstep  = evstep->GetPositionVector();  
      if(rstep.z()<-500) continue;
      totLength = evstep->GetTotalTrackLength();
      edep = evstep->GetEdep();
      eloss+=edep;
      int volId = evstep-> GetSensitiveVolumeID();
      TString physName = TString(evstep->GetPhysVolName());
      bool unknown=true;
      // get volume code
      int vol=-1;
      for(unsigned iname=0; iname<knownName.size(); ++ iname) if( physName==knownName[iname]) {
        vol=int(iname);
        unknown=false;
      }
      if(unknown) printf(" \n unknown name %s !!!! \n",physName.Data() );
      TString procName = TString(evstep->GetProcessName());
      
      if(abs(pId)==13) {
        ntStepMu->Fill(jstep,pId,vol,edep,evstep->GetStepLength(),totLength,eloss,rstep.z(),rstep.Perp());
        edepMu=edep;
        totLengthMu=totLength;
        elossMu=eloss;
        rMuStep = rstep;
      }
      else if(pId==0) {
      //optical photons don't have a pdgID so mage sets them to 0
      // photon hits the PMT
        int hit=0;
        if(physName.Contains("physicalPMT")&& procName.Contains("WLS") && pId == 0 && edep != 0) hit=1;
        ntStepOpt->Fill(jstep,pId,hit,edep,evstep->GetStepLength(),totLength,eloss,rstep.z(),rstep.Perp());
      }
      else ntStepOth->Fill(jstep,pId,vol,edep,evstep->GetStepLength(),totLength,eloss,rstep.z(),rstep.Perp());
    }
    ntPrim->Fill(eventSteps->GetNSteps(),energy,rvec.z(),r,theta,edepMu,totLengthMu,elossMu,rMuStep.z(),rMuStep.Perp());
  }
  outFile->Write();
  outFile->Close();
  cout<<" anaCosmic finished  "<<fileName<<endl;
  return 0;
}
