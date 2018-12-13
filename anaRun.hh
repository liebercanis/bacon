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
#include "TRandom3.h"

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
  enum {minLength=7,maxHalfLength=200};
  enum {baseLineHalfWindow = 200}; // even integer
  double fsigma;
  Int_t nSigma;
  Int_t aveWidth;
  Int_t windowSize;

  Double_t baseline[NPMT];
  Double_t sDev[NPMT];


 
  anaRun(){;}
  anaRun(TString tag, Int_t maxEvents);
  virtual ~anaRun(){;}

  TTree* pmtTree;
  TPmtEvent* pmtEvent;
  

  // for TSpectrum baseline
  TSpectrum * spec;
  Double_t* source;
  TRandom3 *ran;
  //
  std::vector<Int_t> findPeaks(std::vector<Double_t> v, Double_t threshold,Double_t sthreshold); //,Int_t startBin, Int_t stopBin);
  hitMap findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma);
  std::vector<Int_t>    windowPeaks(std::vector<Double_t> v, Int_t nwindow, Double_t rms); 
  std::vector<Double_t> SimpleFilter(std::vector<Double_t> in);
  std::vector<Double_t> SimpleLowPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> SimpleHighPassFilter(std::vector<Double_t> in, Double_t a);
  std::vector<Double_t> MovingAverageFilter(std::vector<Double_t> signal,Int_t N);
  std::vector<Double_t> BubbleSort(std::vector<Double_t> A);
  std::vector<Double_t> getTBaseline(TH1D* hPMTRaw, TH1D* hBaseline, Double_t& ave, Double_t& aveSigma);
  std::vector<Double_t> getBaseline(std::vector<Double_t> digi, TH1D* hBaselineFit, TH1D* hBaseline, Double_t& ave, Double_t& aveSigma);
  void getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma);
 
  void plotWave(Int_t iev, Int_t pmtNum, hitMap pmtHits);
  void sumWave(Int_t ipmt);

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
  TNtuple *ntPulse;
  TH1D* hFFT[NPMT];
  TH1D* hHitQ[NPMT];
  TH1D* hNHits[NPMT];
  TH1D* hPMTRaw[NPMT];
  TH1D* hPMTSignal[NPMT];
  TH1D* hSum[NPMT];
  TH1D* hLife[NPMT];
  
  TH1D* hWeight;
  TH1D* hWeightOne;
  TH1D* hNoise;
  TH1D* hFFTNoise;
  TH1D* hBase;
  TH1D* hChargeWindow;
  TH1D* hChargePeak;
  TH1I* hPeakNWidth;
  TH1I* hNWidthCut;
  TH1I* hAllNWidth;
  TH1I* hSlideWidth;
  TH1D* hSlideQSum;
  TH1D* hSlideHigh;
  TH1D* hSlideLow;
  TH1D* hBaseline[NPMT];
  TH1D* hBaselineFit[NPMT];
  TH1D* hHits;
  TH1D* hDeltaV;
  TH1D* hStartTime;
  TH1D* hQMax;
  TH1D* hChargePeakOff;
};

