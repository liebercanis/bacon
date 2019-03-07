#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
#include "MGTMCRun.hh"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TString.h"
#include "TRandom2.h" 
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
  TString sFileName = TString("Detect.")+to_string(time.GetTime())+TString("-")+to_string(time.GetDate())+TString(".root");
  TFile *outFileMaGe = new TFile(sFileName,"recreate");
  //TNtuple *ntupleMap = new TNtuple("ntupleMap","ntupleMap","nphotons:x:y:z:r:theta:ke:eDep");

  TH3D* hMap;// = new TH3D("OpticalMap","OpticalMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);
  TH2D* h2DMapRZ;
  TH2D* h2DMapXY;
  TH2D* h2DMapYZ;

  TH3D* h3DMapDistribution;
  TH2D* h2DMapRZDistribution;
  TH2D* h2DMaXYDistribution;
  TH2D* h2DMapYZDistribution;

  Double_t gridSpacing = 2;//*mm
  Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.,minR=0,maxR=320;
  //Double_t maxX = 1000.,maxY=1000,maxZ=1200,minZ = -800.,minR = 320., maxR = 1000.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;

  TH3D* hOutMapUnscaled = new TH3D("EnergyMapUnscaled","EnergyMapUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);
  
  TH3D* hOutMap = new TH3D("EnergyMap","EnergyMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);
  TH2D* h2DOutMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);
  TH2D* h2DOutMapRZUnscaled = new TH2D("2DOpticalMap_RZ_Unscaled","2DOpticalMap_RZ_Unscaled",nbinsX/2.,minR,maxR,nbinsZ,minZ,maxZ);
  

  TH2D* h2DOutMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* h2DOutMapXYUnscaled = new TH2D("2DOpticalMap_XYUnscaled","2DOpticalMap_XYUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);

  TH2D* h2DOutMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DOutMapYZUnscaled = new TH2D("2DOpticalMap_YZUnscaled","2DOpticalMap_YZUnscaled",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/";
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
  TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
   
  //TString fileName = "SensitiveVolumes_NoOptical-1000";
  //TString fileName = "OpticalMap10mm400MEventsOriginalGeometry";
  //TString fileName = "MaGe.155758-20181001";
  //TString fileName = "MaGe.82655-20181030";
  //TString fileName = "RawLEGEND200.4e9.10mm";
  TString fileName = "RawLEGEND200.4e9.2mm.TwoFiber";
  //TString fileName = "RawLEGEND200.1e7.50mm.Exterior";
 
  if(!fileExist(string(dir+fileName+TString(".root")))){
    cout<<"no file, no cry, strong boy, good hear, look, find, file. "<<endl;
  }
  TFile* mapfile = TFile::Open(dir+fileName+TString(".root"));
  
  Double_t nEvents = 4e9;//1.e7;//4000000000;
  //get Map
  mapfile->ls();
  mapfile->GetObject("OpticalMap_unScaled",hMap);
  mapfile->GetObject("2DOpticalMap_RZUnscaled",h2DMapRZ);
  mapfile->GetObject("2DOpticalMap_XYUnscaled",h2DMapXY);
  mapfile->GetObject("2DOpticalMap_YZUnscaled",h2DMapYZ);

  /*
  mapfile->ls();
  mapfile->GetObject("OpticalMap",hMap);
  mapfile->GetObject("2DOpticalMap_RZ",h2DMapRZ);
  mapfile->GetObject("2DOpticalMap_XY",h2DMapXY);
  mapfile->GetObject("2DOpticalMap_YZ",h2DMapYZ);
  */

  
  //TFile* distFile = TFile::Open(dir+TString("distributionMapLGND200.2mm.TwoFiber.root"));
  TFile* distFile = TFile::Open(dir+TString("distributionMapLGND200.2mm.TwoFiber.root"));
  //TFile* distFile = TFile::Open(dir+TString("distributionMapLGND200.10mm.root"));
  //TFile* distFile = TFile::Open(dir+TString("distributionMapLGND200.50mm.Exterior.root"));

  distFile->GetObject("OpticalMap_Distribution",h3DMapDistribution);    
  distFile->GetObject("2DOpticalMap_RZDistribution",h2DMapRZDistribution);    
  distFile->GetObject("2DOpticalMap_XYDistribution",h2DMaXYDistribution);    
  distFile->GetObject("2DOpticalMap_YZDistribution",h2DMapYZDistribution);
  //Double_t Ndist = h3DMapDistribution->GetEntries();


  outFileMaGe->cd();

  Double_t scintYield = 40.;//40 photons/KeV
  for(int i = 0; i <= hMap->GetNbinsX();i++){
    for(int j = 0; j <= hMap->GetNbinsY();j++){
      for(int k = 0; k <= hMap->GetNbinsZ();k++){
        Double_t binVal = hMap->GetBinContent(i,j,k);
        Double_t weight = h3DMapDistribution->GetBinContent(i,j,k);
/*        
        Double_t z = k*gridSpacing+minZ;
        Double_t y = j*gridSpacing-maxY;
        Double_t x = i*gridSpacing-maxX;
        
        Int_t newBin;
        //cout<<"("<<y<<","<<z<<")"<<endl;
        Double_t top = 203,bot = 19,shift =(top-bot)/4;
        if(z>(top-shift)){//195//107
          z -= shift;
          newBin = hMap->FindBin(x,y,z)-1;
        }  
        else if(z<(bot+shift)){//20//107
          z += shift;
          newBin = hMap->FindBin(x,y,z)-1;
        }
        else continue;
        //      cout<<"\t new ("<<y<<","<<z<<")"<<endl;
        Double_t prob = binVal/(weight*nEvents);
        Double_t eThresh = 1./(prob*scintYield);
        if(binVal == 0)continue;
        if(weight == 0)continue;
        //hOutMap->SetBinContent(newBin,eThresh);
        hOutMapUnscaled->SetBinContent(newBin,prob);
*/
        ///*
        if(binVal == 0)continue;
        if(weight == 0)continue;

        //Double_t prob = binVal/weight;      
        Double_t prob = binVal/(weight*nEvents);      
        Double_t eThresh = 1./(prob*scintYield);        
        //hOutMap->SetBinContent(i,j,k,eThresh);
        hOutMapUnscaled->SetBinContent(i,j,k,prob);
        //*/
      }
    }
  }

///*
  for(int i = 0; i <= h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZ->GetBinContent(i,j);
      Double_t weight = h2DMapRZDistribution->GetBinContent(i,j);
      Double_t z = j*gridSpacing+minZ;
      Double_t r = i*gridSpacing+minR;
      ///*
      Int_t newBin;
      //cout<<"("<<y<<","<<z<<")"<<endl;
      Double_t delta = 1,top = 203+delta,bot = 19-delta,shift =(top-bot)/4;
      if(z>(top-shift)){//195//107
        z -= shift;
        newBin = h2DMapRZ->FindBin(r,z)-1;
      }  
      else if(z<(bot+shift)){//20//107
        z += shift;
        newBin = h2DMapRZ->FindBin(r,z)-1;
      }
      else continue;
      Double_t prob = binVal/(weight*nEvents);
      Double_t eThresh = 1./(prob*scintYield);
      if(binVal == 0)continue;
      if(weight == 0)continue;
      h2DOutMapRZ->SetBinContent(newBin,eThresh);
      h2DOutMapRZUnscaled->SetBinContent(newBin,prob);
      //*/
      /*
      if(binVal == 0) binVal = h2DMapRZ->Interpolate(r,z);
      Double_t prob = binVal/(weight*nEvents);
      Double_t eThresh = 1./(prob*scintYield);
      //if(binVal == 0)continue;
      if(binVal == 0) prob = 1e-4;
      if(weight == 0)continue;
      h2DOutMapRZ->SetBinContent(i,j,eThresh);
      h2DOutMapRZUnscaled->SetBinContent(i,j,prob);
      */
    }
  }
  for(int i = 0; i <= h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXY->GetBinContent(i,j);
      Double_t weight = h2DMaXYDistribution->GetBinContent(i,j);

      if(binVal == 0)continue;
      if(weight == 0)continue;

      //Double_t prob = binVal/weight;
      Double_t prob = binVal/(weight*nEvents);
      Double_t eThresh = 1./(prob*scintYield);
      
      h2DOutMapXY->SetBinContent(i,j,eThresh);
      h2DOutMapXYUnscaled->SetBinContent(i,j,prob);
    }
  }

  for(int i = 0; i <= h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZ->GetBinContent(i,j);
      Double_t weight = h2DMapYZDistribution->GetBinContent(i,j);

      if(binVal == 0)continue;
      if(weight == 0)continue;
      //Double_t prob = binVal/weight;
      Double_t prob = binVal/(weight*nEvents);
      Double_t eThresh = 1./(prob*scintYield);      
      //if(eThresh > 1e4) eThresh = 1.e4;
      h2DOutMapYZ->SetBinContent(i,j,eThresh);
      h2DOutMapYZUnscaled->SetBinContent(i,j,prob);
    }
  }
//*/
  outFileMaGe->Write();
  outFileMaGe->Close();
cout<<"root -l "<<sFileName<<endl;
  return 0;
}
