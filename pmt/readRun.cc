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

// file scope variables 
TTree *pmtTree;
TPmtEvent *pmtEvent;
bool TwoPmts;
int nDuplicates =0;
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
    if(TwoPmts) volt2Vec.push_back(volt2);
    nlines++;
  }
  in.close();

  bool duplicate=false;
  if(pmtTree->GetEntries() !=0 &&
      pmtEvent->volt1[0] == volt1Vec[0]     && pmtEvent->volt1[10] == volt1Vec[10]     && 
      pmtEvent->volt1[100] == volt1Vec[100] && pmtEvent->volt1[110] == volt1Vec[110] &&
      pmtEvent->volt1[150] == volt1Vec[150] && pmtEvent->volt1[160] == volt1Vec[160] &&
      pmtEvent->volt1[200] == volt1Vec[200] && pmtEvent->volt1[210] == volt1Vec[210] 
    ) duplicate = true;

  if(!duplicate) {
    pmtEvent->volt1 = volt1Vec;
    pmtEvent->volt2 = volt2Vec;
    pmtEvent->time = timeVec;
    pmtEvent->event=ievent;
    pmtTree->Fill();
  } else ++nDuplicates;

  // dont want to claer pmtEvent->clear();
  //printf(" have read %i lines and %llu entries \n",nlines,pmtTree->GetEntries());

  return nlines;
    
}

void readRun(TString dirTag="runNov2018_1", unsigned maxEvents=0)
{
  TwoPmts = false;
  Int_t firstRun=0;
  Int_t lastRun = firstRun;
  printf("readRun tag %s  maxEvents %u first \n",dirTag.Data(),maxEvents);
  if(TwoPmts) printf("\t reading 2 pmts \n");
  else  printf("\t WARNING... reading 1 pmt only \n");
  for(Int_t irun  = firstRun; irun <=lastRun;irun++){
    // open ouput file and make some histograms
    TString outFileName; outFileName.Form("rootData/baconRun_%s_%u.root",dirTag.Data(),maxEvents);
    TFile *outfile = new TFile(outFileName,"recreate");
    outfile->cd();
    printf(" opening output file %s \n",outFileName.Data());

    // ttree
    pmtTree = new TTree("pmtTree","pmtTree");
    pmtEvent  = new TPmtEvent();
    pmtTree->Branch("pmtEvent",&pmtEvent);

    // get list of files
    TString dirName;
    dirName.Form("rawData/%s",dirTag.Data());
    std::string directory(dirName.Data());

    std::vector<Int_t> eventNumbers;
    if(!getEventsInDirectory(directory,eventNumbers)) return;
    if(maxEvents==0) maxEvents = eventNumbers.size();
    Int_t nlines=0;
    for( unsigned ievent = eventNumbers[0]; ievent < maxEvents ; ++ievent ) {
      TString fname;
      fname.Form("%s_%i.txt",tag.Data(),ievent);
      TString fullFileName = dirName + TString("/")+fname;
      if(ievent%100==0) cout << ievent << "  " << fname << " tree size  " << pmtTree->GetEntries()  << " duplicates " << nDuplicates << endl;
      nlines += readEvent(ievent,fullFileName);
    }
    printf(" total of lines %i total number of events is %i number duplicates %i \n",nlines,int(pmtTree->GetEntries()),nDuplicates);
    outfile->Write();
  }
}
