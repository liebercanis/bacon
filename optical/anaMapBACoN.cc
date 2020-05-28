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

/* 
  Double_t gridSpacing = 10;//mm
  Double_t maxX = 320.,maxY=320.0,maxZ=600,minZ = -200.;
*/
  Double_t gridSpacing = 5;//mm
  Double_t maxX = (19./2)*2.54*10,maxY=(19./2)*2.54*10,maxZ=23.*2.54*10/2.,minZ =-23.*2.54*10./2.;
  Int_t nbinsX = 2*maxX/gridSpacing,nbinsY = 2*maxY/gridSpacing,nbinsZ = (maxZ-minZ)/gridSpacing;

  TH3D* hOutMap = new TH3D("EnergyMap","EnergyMap",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);
  TH3D* hOutMapUnscaled = new TH3D("EnergyMapUnscaled","EnergyMapUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY,nbinsZ,-maxZ,maxZ);
  
  TH2D* h2DOutMapRZ = new TH2D("2DOpticalMap_RZ","2DOpticalMap_RZ",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);
  TH2D* h2DOutMapRZUnscaled = new TH2D("2DOpticalMap_RZ_Unscaled","2DOpticalMap_RZ_Unscaled",nbinsX/2.,0,maxX,nbinsZ,minZ,maxZ);
  

  TH2D* h2DOutMapXY = new TH2D("2DOpticalMap_XY","2DOpticalMap_XY",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);
  TH2D* h2DOutMapXYUnscaled = new TH2D("2DOpticalMap_XYUnscaled","2DOpticalMap_XYUnscaled",nbinsX,-maxX,maxX,nbinsY,-maxY,maxY);

  TH2D* h2DOutMapYZ = new TH2D("2DOpticalMap_YZ","2DOpticalMap_YZ",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);
  TH2D* h2DOutMapYZUnscaled = new TH2D("2DOpticalMap_YZUnscaled","2DOpticalMap_YZUnscaled",nbinsY,-maxY,maxY,nbinsZ,minZ,maxZ);

  TH1D* hWeight = new TH1D("weight","weight",nbinsX/2.,0.,maxX);
  TH1D* hWeightYZ = new TH1D("weightYZ","weightYZ",nbinsY+1,-maxY,maxY);
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/MaGe/bin/Linux-g++/";
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/";
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/array/";
  //TString dir = "/mnt/mjdDisk1/Majorana/users/nmcfadden/RooT/";
  TString dir = "/home/nmcfadden/RooT/mage/";
   
  //TString fileName = "SensitiveVolumes_NoOptical-1000";
  //TString fileName = "OpticalMap10mm400MEventsOriginalGeometry";
  //TString fileName = "MaGe.155758-20181001";
  //TString fileName = "BACoN.81128-20181116";
  TString fileName = "RawBACoN_1PMT.4e9"; 
cout<<"nbinsX "<<nbinsX<<", nbinsY "<<nbinsY<<", nbinsZ "<<nbinsZ<<endl;
  if(!fileExist(string(dir+fileName+TString(".root")))){
    cout<<"no file, no cry, strong boy, good hear, look, find, file. "<<endl;
  }
  TFile* mapfile = TFile::Open(dir+fileName+TString(".root"));
  //get Map
  mapfile->ls();
  mapfile->GetObject("OpticalMap",hMap);
  mapfile->GetObject("2DOpticalMap_RZ",h2DMapRZ);
  mapfile->GetObject("2DOpticalMap_XY",h2DMapXY);
  mapfile->GetObject("2DOpticalMap_YZ",h2DMapYZ);
  
  TFile* distFile = TFile::Open(dir+TString("ProbDistBACON1PMT1.5*mm.112455-20190305.root"));
  distFile->GetObject("OpticalMap_Distribution",h3DMapDistribution);    
  distFile->GetObject("2DOpticalMap_RZDistribution",h2DMapRZDistribution);    
  distFile->GetObject("2DOpticalMap_XYDistribution",h2DMaXYDistribution);    
  distFile->GetObject("2DOpticalMap_YZDistribution",h2DMapYZDistribution);
  Double_t Ndist = h3DMapDistribution->GetEntries();

  outFileMaGe->cd();
  
  TRandom2 rand;
  Double_t Nrand = 1.e6;
  //for(int i = 0; i < h2DMapRZ->GetNbinsX();i++){
  for(int i = 0; i < Nrand;i++){
    Double_t weight = rand.Rndm();
    weight = sqrt(weight);
    hWeight->Fill(weight*maxX,1./Nrand);
  }
  for(int i = 1; i <=hWeightYZ->GetNbinsX();i++){
    Double_t y = hWeightYZ->GetBinCenter(i);//gridSpacing*(i)-maxY;
    Double_t val = sqrt(1-(y/maxY)*(y/maxY));
    hWeightYZ->SetBinContent(i,val );
//    cout<<y<<" "<<val<<endl;
  }
  ///*
  Double_t maxWeight = hWeightYZ->Integral();
  for(int i = 0; i <= hWeightYZ->GetNbinsX();i++){
    hWeightYZ->SetBinContent(i,hWeightYZ->GetBinContent(i)/maxWeight);
  }
  //*/
  /*
  Double_t maxWeight = hWeight->GetMaximum();
  for(int i = 0; i < hWeight->GetNbinsX();i++){
    hWeight->SetBinContent(i+1,hWeight->GetBinContent(i+1)/maxWeight);
  }
  */

  Double_t scintYield = 40.;//40 photons/KeV
  for(int i = 0; i <= hMap->GetNbinsX();i++){
    for(int j = 0; j <= hMap->GetNbinsY();j++){
      for(int k = 0; k <= hMap->GetNbinsZ();k++){
        Double_t binVal = hMap->GetBinContent(i,j,k);
        Double_t x = gridSpacing*(i)-maxX;
        Double_t y = gridSpacing*(j)-maxY;
        Int_t r = sqrt(x*x+y*y);
        
        //Double_t weight = hWeight->GetBinContent(hWeight->FindBin(r) );
        Double_t weight = h3DMapDistribution->GetBinContent(i,j,k)/Ndist;

        if(binVal == 0)continue;
        if(weight == 0)continue;

        Double_t prob = binVal/(weight);      
        Double_t eThresh = 1./(prob*scintYield);        
        hOutMap->SetBinContent(i,j,k,eThresh);
        hOutMapUnscaled->SetBinContent(i,j,k,prob);
      }
    }
  }

  for(int i = 0; i <= h2DMapRZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapRZ->GetNbinsY();j++){
      Double_t binVal = h2DMapRZ->GetBinContent(i,j);
      //Double_t weight = hWeight->GetBinContent(i);
      Double_t weight = h2DMapRZDistribution->GetBinContent(i,j)/Ndist;
      
      //Double_t prob = nbinsZ*binVal/(weight);
      Double_t prob = binVal/weight;
      Double_t eThresh = 1./(prob*scintYield);
      if(binVal == 0)continue;
      if(weight == 0)continue;
      h2DOutMapRZ->SetBinContent(i,j,eThresh);
      h2DOutMapRZUnscaled->SetBinContent(i,j,prob);
    }
  }

  for(int i = 0; i <= h2DMapXY->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapXY->GetNbinsY();j++){
      Double_t binVal = h2DMapXY->GetBinContent(i,j);
      Double_t x = gridSpacing*(i)-maxX;//hMap->ProjectionX()->GetBinCenter(i);
      Double_t y = gridSpacing*(j)-maxY;//hMap->ProjectionY()->GetBinCenter(j);
      Int_t r = sqrt(x*x+y*y);

      //Double_t weight = hWeight->GetBinContent(hWeight->FindBin(r)+1 );
      Double_t weight = h2DMaXYDistribution->GetBinContent(i,j)/Ndist;

      if(binVal == 0)continue;
      if(weight == 0)continue;

      //cout<<binVal<<" "<<1./(scintYield*binVal)<<" "<<1-weight<<endl;
      //binVal = (weight)/(scintYield*binVal);
      //Double_t prob = nbinsZ*binVal/(weight);
      Double_t prob = binVal/weight;
      Double_t eThresh = 1./(prob*scintYield);
      //binVal = binVal/(weight);
      //if(binVal > 1e4) binVal = 1.e4;
      //if(binVal < 10) binVal = 10;
      h2DOutMapXY->SetBinContent(i,j,eThresh);
      h2DOutMapXYUnscaled->SetBinContent(i,j,prob);
    }
  }

  for(int i = 0; i <= h2DMapYZ->GetNbinsX();i++){
    for(int j = 0; j <= h2DMapYZ->GetNbinsY();j++){
      Double_t binVal = h2DMapYZ->GetBinContent(i,j);
      Double_t y = hWeightYZ->GetBinCenter(i);//gridSpacing*(i)-maxX;//hMap->ProjectionX()->GetBinCenter(i);

      //Double_t weight = hWeightYZ->GetBinContent(hWeightYZ->FindBin(y) );
      Double_t weight = h2DMapYZDistribution->GetBinContent(i,j)/Ndist;

      if(binVal == 0)continue;
      if(weight == 0)continue;
      //cout<<binVal<<" "<<1./(scintYield*binVal)<<" "<<1-weight<<endl;
      //binVal = (weight)/(scintYield*binVal);
      //binVal = binVal/(weight);
//      Double_t eThresh = (weight)/(scintYield*binVal);
      //Double_t prob = nbinsZ*binVal/(weight);
      Double_t prob = binVal/weight;
      Double_t eThresh = 1./(prob*scintYield);      
     // cout<<"binVal = "<<binVal<<", weight = "<<weight<<endl;
     // cout<<"\t eThresh*prob "<<eThresh*prob<<", prob = "<<prob<<", eThresh = "<<eThresh<<", scintYield/eThresh = "<<scintYield/eThresh<<endl;
      if(eThresh > 1e4) eThresh = 1.e4;
      h2DOutMapYZ->SetBinContent(i,j,eThresh);
      h2DOutMapYZUnscaled->SetBinContent(i,j,prob);
    }
  }

  outFileMaGe->Write();
  outFileMaGe->Close();
cout<<"root -l "<<sFileName<<endl;
  return 0;
}
