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
  //TString fileName = TString("BACoN.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TString fileName = TString("MaGe.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(fileName,"recreate");
  //TFile *outFileOpticalMap = new TFile(TString("OpticalMap.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root"),"recreate");
/*
  Double_t gridSpacing = 5;//mm
  Double_t maxX = (19./2)*2.54*10,maxY=(19./2)*2.54*10,maxZ=23.*2.54*10/2.,minZ =-23.*2.54*10./2.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
*/
  Double_t gridSpacing = 50;//mm
  //Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.;
  Double_t maxX = 1000.,maxY=1000,maxZ=1200,minZ = -800.,minR = 320., maxR = 1000.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;

  //TH3D* hMap = new TH3D("OpticalMap","OpticalMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH3D* hMapDistribution = new TH3D("OpticalMap_Distribution","OpticalMap_unScaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  //TH2D* h2DMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);
  TH2D* h2DMapRZDistribution = new TH2D("2DOpticalMap_RZDistribution","2DOpticalMap_RZDistribution",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);

  //TH2D* h2DMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* h2DMapXYDistribution = new TH2D("2DOpticalMap_XYDistribution","2DOpticalMap_XYDistribution",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);

  //TH2D* h2DMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DMapYZDistribution = new TH2D("2DOpticalMap_YZDistribution","2DOpticalMap_YZDistribution",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  cout<<"starting run"<<endl; 
  Int_t totalEvents = 0; 
  for(int k = 0; k < 1 ; k++){
   
    TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/BACoN/bin/Linux-g++/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrOpticalDist/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/LGND_200Orig100M1.1mAttenuation/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
    //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/arrBACON/";
    //TString fileName = "SensitiveVolumesLGND_200Alt1" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig1.1mAttenuation" +to_string(k);
    //TString fileName = "SensitiveVolumesLGND_200Orig" +to_string(k);
    //TString fileName = "SensitiveVolumes" +to_string(k);
    //TString fileName = "PMTCountingTest"+to_string(k);
    //TString fileName = "PMTSensitiveVolume";
    //TString fileName = "OpticalDistribution_"+to_string(k);
    //TString fileName = "SpeedTest";
    TString fileName = "Exterior.Distirbution";
    ///*
    if(!fileExist(string(dir+fileName+TString(".root")))){
      cout<<"processed "<<k<<" files"<<endl;
      break;
    }
    //*/
    TFile* infile = TFile::Open(dir+fileName+TString(".root"));
    //string neventsString = infile->Get("NumberOfEvents")->GetTitle();
    //totalEvents += std::stoi(neventsString,nullptr,10);
    infile->Close();
    outFileMaGe->cd();
    
    TChain *fTree = new TChain("fTree");
    fTree->Add(dir+fileName+TString(".root"));

    
    Long64_t nentries = (Long64_t)fTree->GetEntries();
    //CHANGE!!!TODO TODO
    totalEvents = nentries;
    
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
      if((i+1)%10000 == 0 || i == nentries - 1 ) cout<<"\tprocessed "<<i+1<<" events"<<endl;
      fTree->GetEntry(i);
      TString physName;
      primaries = eventPrimaries->GetStep(0);
      if(primaries == NULL){
        cout<<"null primary"<<endl;
        continue;
      }
      Double_t x = primaries->GetX(),y = primaries->GetY(),z = primaries->GetZ();//,time = primaries->GetT();
      Double_t r = sqrt(x*x+y*y); 
      hMapDistribution->Fill(x,y,z);
      h2DMapRZDistribution->Fill(r,z);
      h2DMapXYDistribution->Fill(x,y);
      h2DMapYZDistribution->Fill(y,z);

    }
  }
  for(int i = 0; i <= hMapDistribution->GetNbinsX();i++){
    for(int j = 0; j <= hMapDistribution->GetNbinsY();j++){
      for(int k = 0; k <= hMapDistribution->GetNbinsZ();k++){
        Double_t binVal = hMapDistribution->GetBinContent(i,j,k);
        hMapDistribution->SetBinContent(i,j,k,binVal/totalEvents);
      }
    }
  }

  for(int i = 0; i <= h2DMapRZDistribution->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapRZDistribution->GetNbinsY();j++){
      Double_t binVal = h2DMapRZDistribution->GetBinContent(i,j);
      h2DMapRZDistribution->SetBinContent(i,j,binVal/totalEvents);
    }
  }

  for(int i = 0; i <= h2DMapXYDistribution->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapXYDistribution->GetNbinsY();j++){
      Double_t binVal = h2DMapXYDistribution->GetBinContent(i,j);
      h2DMapXYDistribution->SetBinContent(i,j,binVal/totalEvents);
    }
  }

  for(int i = 0; i <= h2DMapYZDistribution->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapYZDistribution->GetNbinsY();j++){
      Double_t binVal = h2DMapYZDistribution->GetBinContent(i,j);
      h2DMapYZDistribution->SetBinContent(i,j,binVal/totalEvents);
    }
  }

  outFileMaGe->Write();
  outFileMaGe->Close();

//  cout<<"Total events = "<<totalEvents<<endl;
  cout<<"Histogram has "<<nbinsX<<" "<<nbinsY<<" "<<nbinsZ<<" "<<nbinsX*nbinsY*nbinsZ <<" bins"<<endl;
  cout<<"root -l "<<fileName<<endl;
  return 0;

}
