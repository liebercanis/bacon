#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>//includes std::pair, std::make_pair
#include <valarray>
//
#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort
#include "TSpectrum.h"
#include <TMath.h>
#include <TDatime.h>
#include "TRandom2.h" 

//
#include "TPmtEvent.hxx"
#include "TPmtSimulation.hxx"
#include "TPmtRun.hxx"
#include "TBaconEvent.hxx"
#include "TPulse.hxx"



class anaNRun
{
  public:
    anaNRun(Int_t dsNum,Int_t runStart,Int_t runStop);
    ~anaNRun() {;};
    std::vector<Double_t> monthVector = {0,31,59,90,120,151,181,212,243,273,304,334,365};
    std::vector<Double_t> runTime;
    
    std::vector<Double_t> tripletLifetime;
    std::vector<Double_t> tripletLifetimeError;
    
    std::vector<Double_t> tripletLifetimePMT0;
    std::vector<Double_t> tripletLifetimeErrorPMT0;
    std::vector<Double_t> tripletLifetimePMT1;
    std::vector<Double_t> tripletLifetimeErrorPMT1;
    
    std::vector<Double_t> tripletLifetimePMTSum0;
    std::vector<Double_t> tripletLifetimeErrorPMTSum0;
    std::vector<Double_t> tripletLifetimePMTSum1;
    std::vector<Double_t> tripletLifetimeErrorPMTSum1;

    
    TPmtRun *    pmtRun;
    TPmtEvent *  pmtEvent;
    TTree* TBacon;
    TBaconEvent *baconEvent;
    

    TH1D * hTripletvTime;
    TH1D * hTripletvTimeSum0;
    TH1D * hTripletvTimeSum1;
    TH1D * hTripletvTimePMT0;
    TH1D * hTripletvTimePMT1;

};
