#include <sstream>
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
#include <math.h>
#include <TDatime.h>
#include "TRandom2.h" 

//
#include "TPmtEvent.hxx"
#include "TPmtSimulation.hxx"

class genPulses{
  public:
    genPulses() {;};
    genPulses(Int_t maxEvents);
    ~genPulses() {;};
    TTree* genTree;
    TH1D * hSignal[2];
    TPmtSimulation* pmtSimulation;
    TPmtEvent* pmtEvent;
    std::vector<Double_t> PulseStartTime(Int_t nPhotons,Double_t randSeed);
    Double_t convolve(Double_t *x, Double_t *par);
    void conv(double s=1.e-9,double t1=2.e-9,double t2=60e-9,double t12=0.75,double mean = 10.0e-9);

};
