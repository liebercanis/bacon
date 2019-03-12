#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
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


int main(){

  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  TString fileName = TString("BACoN.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
//  TFile *outFileOpticalMap = new TFile(TString("OpticalMap.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root"),"recreate");
  
  TNtuple *ntpleStepPMT = new TNtuple("stepPMT","stepPMT","event:nstep:edep:ke:x:y:z:r:t:eventN");
  TNtuple *ntplePMT = new TNtuple("primaryPMT","primaryPMT","x:y:z:r:theta:edep:px:py:pz:nsteps:eventN");
  TNtuple *ntpleStepLAr = new TNtuple("stepLAR","LArCore","event:nstep:edep:ke:x:y:z:r:t:eventN");
  TNtuple *ntpleLAr = new TNtuple("primaryLAr","primaryLAr","x:y:z:r:theta:edep:px:py:pz:nsteps:eventN");

  std::map<std::vector<Int_t>,Double_t> prob_map;
  Double_t gridSpacing = 5;//*mm
  Double_t maxX = (19./2)*2.54*10,maxY=(19./2)*2.54*10,maxZ=23.*2.54*10/2.,minZ =-23.*2.54*10./2.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
  TH3D* hMap = new TH3D("OpticalMap","OpticalMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapUnscaled = new TH3D("OpticalMap_unScaled","OpticalMap_unScaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  TH2D* h2DMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZUnscaled = new TH2D("2DOpticalMap_RZUnscaled","2DOpticalMap_RZUnscaled",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);

  TH2D* h2DMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* h2DMapXYUnscaled = new TH2D("2DOpticalMap_XYUnscaled","2DOpticalMap_XYUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);

  TH2D* h2DMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DMapYZUnscaled = new TH2D("2DOpticalMap_YZUnscaled","2DOpticalMap_YZUnscaled",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  Double_t totalEvents = 0;

  cout<<"starting run"<<endl; 
  
  for(int k = 0; k < 1000 ; k++){
   
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/BACoN/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/LGND_200Orig100M1.1mAttenuation/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
    TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrBACON/";
    //TString fileName = "SensitiveVolumesLGND_200Alt1" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig1.1mAttenuation" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig" +to_string(k);
    //TString fileName = "SensitiveVolumes" +to_string(k);
    //TString fileName = "PMTCountingTest"+to_string(k);
    //TString fileName = "PMTSensitiveVolume";
    //TString fileName = "XenonDopedTest"+to_string(k);
    TString fileName = "Argon"+to_string(k);
    ///*
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      break;
    }
    //*/
    TFile* infile = TFile::Open(dir+fileName+TString(".root"));
    string neventsString = infile->Get("NumberOfEvents")->GetTitle();
    totalEvents += std::stoi(neventsString,nullptr,10);
    infile->Close();
    outFileMaGe->cd();
    
    //TFile *infile = new TFile(dir+fileName,"READONLY");
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));

    
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    
    cout<<"File location at: "<<dir+fileName+TString(".root")<<". Processing "<<nentries<<" entries"<<endl;
    MGTMCEventSteps *eventSteps = 0;
    MGTMCEventSteps *eventPrimaries = 0;
    if(fTree != NULL){
      fTree->SetBranchAddress("eventSteps",&eventSteps);
      fTree->SetBranchAddress("eventPrimaries",&eventPrimaries);
    }
    else{
      cout<<"NULL fTree"<<endl;
      break;
    }
       
    const MGTMCStepData *step,*primaries;
    
    for(Int_t i = 0; i < nentries ; i++){
      if((i+1)%100000 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      fTree->GetEntry(i);
      TString physName;
      TString procName;
      primaries = eventPrimaries->GetStep(0);
      if(primaries == NULL){
        cout<<"null primary"<<endl;
        continue;
      }
      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t px = primaries->GetPx(),py = primaries->GetPy(),pz = primaries->GetPz();
      Double_t r = sqrt(x*x+y*y);
      Double_t theta = std::acos(x/r);
      Double_t edep = 0;
      bool hitPMT = false,hitLAr = false;
      Int_t hitPMTCounter = 0,hitLArCounter = 0;
      if( y < 0) theta += 3.14159265359;
      for (Int_t j = 0; j < eventSteps->GetNSteps();j++){
        step = eventSteps->GetStep(j);
        physName = step->GetPhysVolName();
        procName = step->GetProcessName();
        //optical photons don't have a pdgID so mage sets them to 0
        if(physName.Contains("physicalPMT")&& procName.Contains("WLS") && step->GetParticleID() == 0 && step->GetEdep() != 0){
          hitPMT = true;
          hitPMTCounter++;
//          edep = step->GetEdep();
          ntpleStepPMT->Fill(i,j,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT(),i);
          ntplePMT->Fill(x,y,z,r,theta,step->GetEdep(),px,py,pz,hitPMTCounter,i);
          hMapUnscaled->Fill(x,y,z);
          h2DMapRZUnscaled->Fill(r,z);
          h2DMapXYUnscaled->Fill(x,y);
          h2DMapYZUnscaled->Fill(y,z);
        }
        if(physName.Contains("Detector")){
          hitLAr = true;
          hitLArCounter++;
          edep += step->GetEdep();
          ntpleStepLAr->Fill(i,j,step->GetEdep(),step->GetKineticE(),step->GetX(),step->GetY(),step->GetZ(),sqrt(step->GetX()*step->GetX() +step->GetY()*step->GetY()),step->GetT(),i);
        }
      }
      /*
      if(hitPMT){
        ntplePMT->Fill(x,y,z,r,theta,edep,px,py,pz,hitPMTCounter,i);
        if(edep > 4.e-6 && edep < 7.e-6)
          cout<<"edep = "<<edep<<", "<<procName<<endl;
        //hMap->Fill(x,y,z,1./(100*1000000.));
        hMapUnscaled->Fill(x,y,z);
        h2DMapRZUnscaled->Fill(r,z);
        h2DMapXYUnscaled->Fill(x,y);
        h2DMapYZUnscaled->Fill(y,z);
      }
      */
      if(hitLAr){
        ntpleLAr->Fill(x,y,z,r,theta,edep,px,py,pz,hitLArCounter,i);
      }
    }
  }
  //weight 3Dhistogram
  for(int i = 0; i <= hMap->GetNbinsX();i++){
    for(int j = 0; j <= hMap->GetNbinsY();j++){
      for(int k = 0; k <= hMap->GetNbinsZ();k++){
        Double_t binVal = hMapUnscaled->GetBinContent(i,j,k);
        Double_t x = gridSpacing*(i)-maxX;
        Double_t y = gridSpacing*(j)-maxY;
        Double_t z = gridSpacing*(k)+minZ;
        if(binVal == 0) binVal = hMapUnscaled->Interpolate(x,y,z);
        hMap->SetBinContent(i,j,k,binVal/totalEvents);
      }
    }
  }
  for(int i = 0; i <= h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZUnscaled->GetBinContent(i,j);
       Double_t r = gridSpacing*(i);
       Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0) binVal = h2DMapRZUnscaled->Interpolate(r,z);
      h2DMapRZ->SetBinContent(i,j,binVal/totalEvents);
    }
  }
  for(int i = 0; i <= h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXYUnscaled->GetBinContent(i,j);
      Double_t x = gridSpacing*(i)-maxX;
      Double_t y = gridSpacing*(j)-maxY;
      if(binVal == 0) binVal = h2DMapXYUnscaled->Interpolate(x,y);
      h2DMapXY->SetBinContent(i,j,binVal/totalEvents);
    }
  }

  for(int i = 0; i <= h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZUnscaled->GetBinContent(i,j);
      Double_t y = gridSpacing*(i)-maxY;
      Double_t z = gridSpacing*(j)+minZ;
      if(binVal == 0) binVal = h2DMapYZUnscaled->Interpolate(y,z);
      h2DMapYZ->SetBinContent(i,j,binVal/totalEvents);
    }
  }

//  outFileMaGe->cd();
  //ntpleStep->Write();
  //ntplePMT->Write();
  //ntpleLAr->Write();
  outFileMaGe->Write();
  outFileMaGe->Close();
/*
  outFileOpticalMap->cd();
  hMap->Write();
  hMapUnscaled->Write();
  outFileOpticalMap->Close();
*/
  cout<<"Total events = "<<totalEvents<<endl;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;
}
