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
#include <TH2D.h>
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
#include "TBaconEvent.hxx"
#include "compassEvent.hxx"
#include "TPmtHit.hxx"
#include "TPulse.hxx"
#include "TPmtSimulation.hxx"
#include "TPmtSimMatchStats.hxx"

typedef std::complex<double> Complex;
typedef std::map<Double_t,TPmtHit,std::less<Double_t> >  hitMap;
typedef std::map<Double_t,TPmtHit,std::less<Double_t> >::iterator  hitMapIter;

typedef std::vector<std::pair<unsigned,unsigned> >  peakType;
typedef std::vector<std::pair<unsigned,unsigned> >::iterator  peakTypeIter;

class testCompass {
  public :

    testCompass(TString tag="31_199", Int_t maxEvents=0);
    virtual ~testCompass(){;}

    TChain* fChain;

};

