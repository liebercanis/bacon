//////////////////////////////////////////////////////////
//  M.Gold August 2017
//////////////////////////////////////////////////////////
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
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
//
#include "TPmtEvent.hxx"
#include "TPmtHit.hxx"

typedef std::complex<double> Complex;
typedef std::map<Double_t,TPmtHit,std::greater<Double_t> >  hitMap;
typedef std::map<Double_t,TPmtHit,std::greater<Double_t> >::iterator  hitMapIter;


class anaRun {
public :
  enum {NPMT=2};
  enum {MAXSAMPLES=10000};
  enum {MAXHIT=3};
  enum {minLength=10,maxHalfLength=100};
  double fsigma;
  int aveWidth;
  Double_t deltaT = 0,deltaV = 999,singlePhoton = 0;
  Double_t baseline[NPMT];
  Double_t sDev[NPMT];


 
  anaRun(){;}
  anaRun(TString tag, Double_t sigma);
  virtual ~anaRun(){;}

  TTree* pmtTree;
  TPmtEvent* pmtEvent;
  

  //TSpectrum * spec = new TSpectrum(20000);
  std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); //,Int_t startBin, Int_t stopBin);
  hitMap findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma);
  std::vector<Double_t> SimpleFilter(std::vector<Double_t> in);
  std::vector<Double_t> SimpleLowPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> SimpleHighPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> MovingAverageFilter(std::vector<Double_t> signal,Int_t N);
  std::vector<Double_t> BubbleSort(std::vector<Double_t> A);

  void plotWave(Int_t ientry, Int_t pmtNum, hitMap pmtHits);

  //Double_t convolve(Double_t *x, Double_t *par);
  Int_t nSamples;
  /// The fft class to take the fourier transform.
  TVirtualFFT *fFFT;
  /// The fft class to take the inverse fourier transform.
  TVirtualFFT *fInverseFFT;
  //std::vector<Double_t> FFTFilter(Int_t ipmt,Int_t ievent,Double_t norm);
  //TH1D* FFTFilter(Int_t ipmt,Int_t ievent);
  //std::pair<TH1D*,std::vector<Double_t> > FFTFilter(Int_t ipmt,Int_t ievent,Double_t norm,std::vector<Double_t> signal);
  std::vector<std::complex<double> > FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal);
  std::vector<Double_t > inverseFFT(Int_t ipmt,Int_t ievent, std::vector<std::complex<double> > VectorComplex,std::vector<Double_t> sum);
  // histogram pointers
  TNtuple *ntHit;
  TNtuple *ntNHit;
  TNtuple *ntEvent;
  TNtuple *ntCal;
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hNHits[NPMT];
  TH1D* hPMTSignal[NPMT];

  TH1D* hPMTSum1;
  TH1D* hPMTSum2;
  
  TH1D* hNoise;
  TH1D* hFFTNoise;
  TH1D* hBase;
  TH1D* hChargeWindow;
  TH1D* hChargePeak;
  TH1I* hPeakNWidth;
  TH1I* hAllNWidth;
  TH1D* hBaseline;
  TH1D* hBaselineFFT;
  TH1D* hHits;
  TH1D* hDeltaV;
  TH1D* hStartTime;
  TH1D* hQMax;
  TH1D* hChargePeakOff;
};

