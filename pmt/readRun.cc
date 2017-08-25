////////////////////////////////////////////////////////
#include<stdio.h>
#include <string.h>
#include <dirent.h> 
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
// local 
#include "TPmtEvent.hxx"
TTree *pmtTree;
TPmtEvent *pmtEvent;


void getEventsInDirectory(std::string directory, std::vector<Int_t> &eventNumbers)
{
  DIR *dir;
  struct dirent *dp;
  //class stat st;
    string tag("event");
    printf(" looking in directory %s \n",directory.c_str());
    dir = opendir(directory.c_str());
    while ((dp = readdir(dir)) != NULL) {
        string file_name = dp->d_name;
        if ( strstr(file_name.c_str(), "event" )==NULL ) continue;
        if ( strstr(file_name.c_str(), ".txt" )== NULL ) continue;
        TString fname(file_name.c_str());
        Int_t ievent = TString(fname( fname.First("_")+1, fname.First(".") - fname.First("_"))).Atoi();
        eventNumbers.push_back(ievent);
    }
    closedir(dir);
    printf(" found %i event files \n",Int_t(eventNumbers.size()));
} // GetFilesInDirectory

int readEvent(Int_t ievent, TString fileName)
{
  //cout << " reading file " << fileName << endl;
  pmtEvent->clear();
  pmtEvent->event=ievent;

  //printf(" looking for file %s \n",fileName.Data());

  Int_t nlines = 0;
  Double_t time,volt1,volt2;
  ifstream in;
  in.open(fileName);
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return 0;
  }
  while (1) {
    in >> time >> volt1 >> volt2 ;
    if(!in.good()) break;
    pmtEvent->time.push_back(time);
    pmtEvent->volt1.push_back(volt1);
    pmtEvent->volt2.push_back(volt2);
    nlines++;
  }
  //printf(" have read %i lines and %i entries \n",nlines,pmtEvent->time.size());
  in.close();
  
  pmtTree->Fill();
  return nlines;
    
}

void readRun(Int_t irun=0)
{
  // open ouput file and make some histograms
  TString outFileName; outFileName.Form("baconRun_%i.root",irun);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);

  // get list of files
  TString dirName;
  dirName.Form("data/run_%i",irun);
  std::string directory(dirName.Data());

  std::vector<Int_t> eventNumbers;
  getEventsInDirectory(directory,eventNumbers);
  Int_t nlines=0;
  for( unsigned ievent =0; ievent < eventNumbers.size() ; ++ievent ) {
    TString fname;
    fname.Form("event_%i.txt",ievent);
    //int ievent = TString(fname( fname.First("_")+1, fname.First(".") - fname.First("_"))).Atoi();
    TString fullFileName = dirName + TString("/")+fname;
    if(ievent%100==0) cout << ievent << " tree size  " << pmtTree->GetEntries()  << endl;
    nlines += readEvent(ievent,fullFileName);
  }
  printf(" total of lines %i total number of events is %i \n",nlines,int(pmtTree->GetEntries()));
  outfile->Write();
}
