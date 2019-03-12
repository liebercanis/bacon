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
  TString fileName = "OpticalMapL200.root";

  TFile *outFileMaGe = new TFile(fileName,"recreate");
  
  TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  
  TString inFileNameTH3D = "OpticalMapLEGEND200.4e9_2mm.TwoFiber.5Rows_RawProb.root";//"OpticalMapLEGEND200.4e9_10mm.5Rows_RawProb.root";//OpticalMapLEGEND200.4e9_2mm.TwoFiber.5Rows_RawProb.root
  TString inFileNameTH2D = "OpticalMapLEGEND200.4e9_50mm.1e7.Exterior.root";

  /*
  Double_t gridSpacing = 10;
  Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;
  TH3D * h3D = new TH3D("","",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);

  gridSpacing = 50;
  maxX = 1000.;maxY=1000;maxZ=1200;minZ = -800.;minR = 320.; maxR = 1000.;
  nbinsX = 2*maxX/gridSpacing;
  TH2D* h2D = new TH2D("","",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);
  */
  TH3D * h3D;
  TH2D* h2D;
  TFile* fileTH3D = TFile::Open(dir+inFileNameTH3D);
  fileTH3D->GetObject("EnergyMapUnscaled",h3D);
  h3D->SetNameTitle("ProbMapInterior","ProbMapInterior");
  TFile* fileTH2D = TFile::Open(dir+inFileNameTH2D);
  fileTH2D->GetObject("2DOpticalMap_RZ_Unscaled",h2D);
  h2D->SetNameTitle("ProbMapExterior","ProbMapExterior");
  outFileMaGe->cd();
  h3D->Write();
  h2D->Write();
  outFileMaGe->Write();
  outFileMaGe->Close();
  cout<<"root -l "<<fileName<<endl;
}

