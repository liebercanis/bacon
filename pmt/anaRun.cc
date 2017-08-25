////////////////////////////////////////////////////////
#include<stdio.h>
#include <string.h>
#include <sstream>
#include <iostream>
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
TTree *tree;


void anaRun(Int_t irun=0)
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/baconRun_%i.root",irun);
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  // get tree from file 
  fin->GetObject("pmtTree",tree);
  Long64_t nentries = tree->GetEntries();
  cout << " number of entries is " << nentries << endl;

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  tree->SetBranchAddress("pmtEvent", &pmtEvent);
  
  // open file for histograms
  TString outFileName; outFileName.Form("baconRunAna_%i.root",irun);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outFileName.Data());
  // loop over entries 
  for (Long64_t ientry=0; ientry<nentries; ientry++) {
        tree->GetEntry(ientry);
        cout << pmtEvent->event << endl;
  }
  
  outfile->Write();
}
