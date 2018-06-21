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

TString tag;

bool getEventsInDirectory(std::string directory, std::vector<Int_t> &eventNumbers)
{
  DIR *dir;
  struct dirent *dp;
  //class stat st;
    printf(" looking in directory %s \n",directory.c_str());
    dir = opendir(directory.c_str());
    while ((dp = readdir(dir)) != NULL) {
        string file_name = dp->d_name;
        //if ( strstr(file_name.c_str(), tag.c_str() )== NULL ) continue;
        if ( strstr(file_name.c_str(), ".txt" )== NULL ) continue;
        TString fname(file_name.c_str());
        Int_t ievent = TString(fname(0,fname.First("_"))).Atoi();
        tag = TString( fname(0,fname.First("_") ));
        eventNumbers.push_back(ievent);
    }
    closedir(dir);
    printf(" In directory %s found %i event files tag %s \n",directory.c_str(),Int_t(eventNumbers.size()),tag.Data());
    if(eventNumbers.size()==0) return false;
    return true;
} // GetFilesInDirectory

int readEvent(Int_t ievent, TString fileName)
{
  //cout << " reading file " << fileName << endl;

  //pmtEvent->clear();
  //pmtEvent->event=ievent;
  //printf(" looking for file %s \n",fileName.Data());

  string line;
  Int_t nlines = 0;
  Double_t time,volt1,volt2;
  std::vector<Double_t> timeVec,volt1Vec,volt2Vec;
  ifstream in;
  in.open(fileName);
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return 0;
  }
  while (in.good()) {
    in >> time >> volt1 >> volt2 ;
    if(in.eof()||in.fail()||in.bad()) break;
    //printf("%.11f %.8f %.8f \n",time,volt1,volt2);
    timeVec.push_back(time);
    volt1Vec.push_back(volt1);
    volt2Vec.push_back(volt2);
    /*
    pmtEvent->time.push_back(time);
    pmtEvent->volt1.push_back(volt1);
    pmtEvent->volt2.push_back(volt2);
    */
    nlines++;
  }
  in.close();


  pmtEvent->volt1 = volt1Vec;
  pmtEvent->volt2 = volt2Vec;
  pmtEvent->time = timeVec;
  pmtEvent->event=ievent;
  pmtTree->Fill();

  pmtEvent->clear();
  //printf(" have read %i lines and %llu entries \n",nlines,pmtTree->GetEntries());

  return nlines;
    
}

void readRun(Int_t firstRun=0, unsigned maxEvents=0)
{
  Int_t lastRun = firstRun;
  printf("readRun tag %s  maxEvents %u first run %i last run %i \n",tag.Data(),maxEvents,firstRun,lastRun);
  for(Int_t irun  = firstRun; irun <=lastRun;irun++){
  // open ouput file and make some histograms
  TString outFileName; outFileName.Form("rootData/baconRun_%i_%u.root",irun,maxEvents);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outFileName.Data());

  // ttree
  pmtTree = new TTree("pmtTree","pmtTree");
  pmtEvent  = new TPmtEvent();
  pmtTree->Branch("pmtEvent",&pmtEvent);

  // get list of files
  TString dirName;
  dirName.Form("rawData/run_%i",irun);
  std::string directory(dirName.Data());

  std::vector<Int_t> eventNumbers;
  if(!getEventsInDirectory(directory,eventNumbers)) return;
  if(maxEvents==0) maxEvents = eventNumbers.size();
  Int_t nlines=0;
  for( unsigned ievent = eventNumbers[0]; ievent < maxEvents ; ++ievent ) {
    TString fname;
    fname.Form("%s_%i.txt",tag.Data(),ievent);
    TString fullFileName = dirName + TString("/")+fname;
    if(ievent%100==0) cout << ievent << "  " << fname << " tree size  " << pmtTree->GetEntries()  << endl;
    nlines += readEvent(ievent,fullFileName);
  }
  printf(" total of lines %i total number of events is %i \n",nlines,int(pmtTree->GetEntries()));
  outfile->Write();
  }
}
