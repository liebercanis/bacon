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
  TString outFileName = TString("ntupleReader.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(outFileName,"recreate");

  TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString fileName = "Radon222.86.Cryostat.2e5";

  if(!fileExist(string(dir+fileName+TString(".root")))){
    cout<<"file does not exist"<<endl;
    return 0;
  }

  TFile* infile = TFile::Open(dir+fileName+TString(".root"));

  TTree *ntuple = (TTree*)infile->Get("stepLAR");
  Float_t edep,x,y,z,r,nPhotons;
  
  ntuple->SetBranchAddress("x",&x);
  ntuple->SetBranchAddress("y",&y);
  ntuple->SetBranchAddress("z",&z);
  ntuple->SetBranchAddress("r",&r);
  ntuple->SetBranchAddress("edep",&edep);
  //ntuple->SetBranchAddress("nPhotons",&nPhotons);

  Double_t gridSpacing = 1;//*mm
  Double_t maxX = (19./2)*2.54*10,maxY=(19./2)*2.54*10,maxZ=23.*2.54*10/2.,minZ =-23.*2.54*10./2.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
  TH2D* hEdepYZWeighted = new TH2D("eDep_yzWeighted","eDep_yzWeighted",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* hEdepXYWeighted = new TH2D("eDep_xyWeighted","eDep_xyWeighted",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* hEdepRZWeighted = new TH2D("eDep_rzWeighted","eDep_rzWeighted",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);
  TH2D* hEdepYZ = new TH2D("eDep_yz","eDep_yz",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* hEdepXY = new TH2D("eDep_xy","eDep_xy",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* hEdepRZ = new TH2D("eDep_rz","eDep_rz",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);



  Int_t entries = (Int_t)ntuple->GetEntries();
  for(Int_t i = 0; i < .1*entries; i++){
    if((i+1)%10000000 == 0 || i == entries - 1 ) cout<<"\tprocessed "<<i+1<<" events out of "<<entries<<endl;
    ntuple->GetEntry(i);
    hEdepYZ->Fill(y,z);
    hEdepXY->Fill(x,y);
    hEdepRZ->Fill(r,z);
    hEdepYZWeighted->Fill(y,z,edep);
    hEdepXYWeighted->Fill(x,y,edep);
    hEdepRZWeighted->Fill(r,z,edep);
  }
  for(int i = 0; i < hEdepYZ->GetNbinsX();i++){
    for(int j = 0; j < hEdepYZ->GetNbinsY();j++){
      Double_t weight = hEdepYZ->GetBinContent(i,j);
      if(weight == 0) continue;
      Double_t val = hEdepYZWeighted->GetBinContent(i,j);
      hEdepYZWeighted->SetBinContent(i,j,val/weight);
    }
  }
  for(int i = 0; i < hEdepXY->GetNbinsX();i++){
    for(int j = 0; j < hEdepXY->GetNbinsY();j++){
      Double_t weight = hEdepXY->GetBinContent(i,j);
      if(weight == 0) continue;
      Double_t val = hEdepXYWeighted->GetBinContent(i,j);
      hEdepXYWeighted->SetBinContent(i,j,val/weight);
    }
  }
  for(int i = 0; i < hEdepRZ->GetNbinsX();i++){
    for(int j = 0; j < hEdepRZ->GetNbinsY();j++){
      Double_t weight = hEdepRZ->GetBinContent(i,j);
      if(weight == 0) continue;
      Double_t val = hEdepRZWeighted->GetBinContent(i,j);
      hEdepRZWeighted->SetBinContent(i,j,val/weight);
    }
  }
  outFileMaGe->cd();
  hEdepYZWeighted->Write();
  hEdepXYWeighted->Write();
  hEdepRZWeighted->Write();
  infile->Close();
  outFileMaGe->Write();
  outFileMaGe->Close();
  cout<<"root -l "<<outFileName<<endl;
  return 0;
}
