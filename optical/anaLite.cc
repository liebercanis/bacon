#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
#include "MGTMCEventHeader.hh"
#include "io/MGOutputMCOpticalRun.hh"
#include "G4OpBoundaryProcess.hh"
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

//static double rMaximum = 240.0;

int main(int argc, char *argv[])
{
  Double_t e128 = 2.*TMath::Pi()*200/128.0*1E-6; 
  Double_t e442 = 2.*TMath::Pi()*200/442.0*1E-6; 
  printf(" constants e128 %.2E e442 MeV %.2E MeV \n",e128,e442);
  
  
  if(argc<2) {
    std::cout << "Usage: " << argv[0] << " <input file name without .root extension in directory ../ >" << '\n';
    return 0;
  }
  Long64_t nprocess=0;
  if(argc>2) nprocess = atoi(argv[2]);
  TString fileName = argv[1];
  TString dir = "/data2/mgold/MaGe_data/";
  TString fullFileName= dir+fileName+TString(".root");
  TFile* inFile = TFile::Open(fullFileName,"readonly");
  if(!inFile){
    printf(" cannot find input file %s\n",fullFileName.Data()); 
    return 0;
  }
  printf(" reading input file %s\n",fullFileName.Data()); 


  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  TString outFileName = TString("anaLight-")+fileName+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opened output %s \n",outFile->GetName());

  TString mapFileName = "OpticalMapBACoN.1e10";
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


  outFile->cd();
  TH1D* hYield = new TH1D("Yield","Yield of photons ",1000,0,10);
  TH1D* hTrkWeight = new TH1D("TrkWeight"," track weight ",10,0,10);
  TH1D* hPhotonsHit = new TH1D("PhotonsHit"," photon count hit ",200,0,200);
  TH1D* hPhotonsAll = new TH1D("PhotonsAll"," photon  count  ",200,0,200);
  TNtuple *ntStepLAr = new TNtuple("stepLAr","LAr","event:nstep:pid:liquid:edep:ke:x:y:z:r:yield:nph");
  TNtuple *ntLAr = new TNtuple("primaryLAr","primary LAr","nsteps:liquid:eprim:eion:x:y:z:r:theta:px:py:pz:phSum:edepSum:aveYield");

  if(nprocess==0) nprocess=nentries;
  cout<<">> Processing "<<nprocess<< " of " << nentries <<" entries"<<endl;

  Int_t noHitCounter = 0;

  for(Long64_t i = 0; i < nprocess ; ++i) {
    //if(i%10 == 0 ) cout<<"\tprocessed "<< i <<" events"<<endl;
    fTree->GetEntry(i);
    const MGTMCStepData *prim0 = eventPrimaries->GetStep(0);
    if(prim0 == NULL){
      cout<<"null prim0"<<endl;
      continue;
    }
    TVector3 pprim = prim0->GetMomentumVector(); 
    TVector3 rprim = prim0->GetPositionVector(); 
    Double_t r = rprim.Perp();
    Double_t x = prim0->GetX(),y = prim0->GetY(),z = prim0->GetZ();//,time = prim0->GetT();
    Double_t px = prim0->GetPx(),py = prim0->GetPy(),pz = prim0->GetPz();
    
    Double_t theta = pprim.Theta()*180.0/TMath::Pi();
    Double_t eprim = prim0->GetKineticE();
    if(i%1 == 0 ) { 
      cout<<"processed "<< i << " pid " << prim0->GetParticleID() << " primary steps " << eventPrimaries->GetNSteps() 
        << " length " << prim0->GetStepLength() << " eprim " << eprim << " theta " << theta  
        << " event steps  " << eventSteps->GetNSteps() << endl;
    }
    // steps loop 
    Double_t sumPhotons = 0, eDepLAr = 0, sumYield = 0.;
    Int_t stepCounter = 0, pmtCounter = 0, nTotalPhotons = 0;
    double eIon=0;
    bool hitLAr=false;
    bool hitGas=false;

    for(int j=0; j<eventSteps->GetNSteps(); ++ j) {
      const MGTMCStepData *step = eventSteps->GetStep(j);
      hTrkWeight->Fill( step->GetTrackWeight() );

      if(step->GetTrackWeight() == 1 && step->GetStepNumber() == 1) nTotalPhotons++;
  
      TVector3 rstep = step->GetPositionVector();
      TString physName = step->GetPhysVolName();
      TString procName = step->GetProcessName();
      
      if(physName.Contains("physicalPMT")&& procName.Contains("WLS") && step->GetParticleID() == 0 && step->GetEdep() != 0) pmtCounter++;
      //cout << i << "  " << j << " " << physName << "  " << procName << endl;
      // map assumes 40000 photons/MeV
      bool inLAr = physName =="Detector";
      if(inLAr) hitLAr=true;  
      bool inGas = physName == "argonGasPhysical";
      if(inGas) hitGas=true;
      if(inLAr||inGas){
        if(abs(step->GetParticleID())==11||abs(step->GetParticleID())==13) eIon+=step->GetEdep(); 
        if(step->GetParticleID() == 0) continue;
        Int_t bin = hMap->FindBin(step->GetX(),step->GetY(),step->GetZ());
        // map stores keV/photon scale nominal 40000 by actual value 
        // yield is photons/keV
        Double_t binVal = hMap->GetBinContent(bin); //map in keV,        
        if(binVal==0) {
          printf(" \t !!!! %f \n",binVal);
          binVal=1.0E9;
        }
        Double_t yield = 1.0/binVal;
        hYield->Fill(yield);
        Double_t stepNpho  = step->GetEdep()*yield;
        sumPhotons += step->GetEdep()*yield;
        sumYield+= yield;
        stepCounter++;
        eDepLAr +=step->GetEdep();
        ntStepLAr->Fill(i,j,step->GetParticleID(),double(inLAr),step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),rstep.Perp(),yield,stepNpho);
      } 
    }

    hPhotonsAll->Fill(sumPhotons);
    if(sumPhotons>0) hPhotonsHit->Fill(sumPhotons);
    if(sumPhotons>1000) printf(" sumPhotons %f eprim %f eIon %f z %f \n",sumPhotons,eprim,eIon,z);
    ntLAr->Fill(eventSteps->GetNSteps(),hitLAr,eprim,eIon,x,y,z,r,theta,px,py,pz,sumPhotons,eDepLAr,sumYield/double(stepCounter));
    if(!hitLAr&&!hitGas) ++noHitCounter;
  }

  mapFile->Close();
  inFile->Close();
  outFile->cd();
  printf("Total events = %lld LAR %lld no hit events %i \n ",nprocess,ntLAr->GetEntriesFast(),noHitCounter);
  printf(" events %lld no hit events %f  mean photons per event with at least one hit %f \n",nprocess,hPhotonsAll->GetBinContent(1), hPhotonsAll->GetMean());
  outFile->ls();
  outFile->Write();
  outFile->Close();
  cout<<" anaLite finished  "<<fileName<<endl;
  return 0;
}
