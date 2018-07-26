/**
** MG, July 2018 
**/
#ifndef TPMTHIT_DEFINED
#define TPMTHIT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store pmt hit 

class TPmtHit: public TNamed {
  public:
    TPmtHit();
    ~TPmtHit(){;}

    void clear();
    // data elements
    Int_t firstBin;
    Int_t lastBin;
    Double_t startTime;
    Double_t peakWidth;
    Double_t qpeak;
    UInt_t peakt;
    Double_t peakMaxTime;
    Int_t peakBin;
    Double_t qsum;

    ClassDef(TPmtHit,3)
};
#endif

