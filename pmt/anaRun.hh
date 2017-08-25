//////////////////////////////////////////////////////////
//  M.Gold August 2017
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
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
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
//
#include "TPmtEvent.hxx"

typedef std::complex<double> Complex;


class anaRun {
public :
  enum {NPMT=2};
  enum {MAXSAMPLES=10000};
  
  anaRun(Int_t irun=0);
  virtual ~anaRun(){;}

  TTree* pmtTree;
  TPmtEvent* pmtEvent;
   
  //std::vector<Int_t> findMaxPeak(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
  //std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); 
  //Int_t findHits(Int_t ipmt, Double_t sum,  std::vector<Int_t> peakTime, std::vector<Double_t> ddigi); 

  double getBaseline(int ipmt ) { return hBase->GetBinContent(ipmt+1); }

  Int_t nFFTSize;
  /// The fft class to take the fourier transform.
  TVirtualFFT *fFFT;
  /// The fft class to take the inverse fourier transform.
  TVirtualFFT *fInverseFFT;
  TH1D* FFTFilter(Int_t ipmt);

  // histogram pointers
  TH1D* hSamples[NPMT];  // include RF
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hNHits[NPMT];
  TH1D* hQMax[NPMT];

  TH1D* hBaseline[NPMT];
  TH1D* hNoise;
  TH1D* hBase;   
};
