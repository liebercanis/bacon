/**
** MG, March 2020 
**/
#ifndef TBACONEVENT_DEFINED
#define TBACONEVENT_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>
#include "TPulse.hxx"

using namespace std;

// class to store info for the event

class TBaconEvent: public TNamed {
	public:
		TBaconEvent();
    //~TBaconEvent();

		void clear();
		// data elements
    Long64_t event;
    Int_t    channel;
    Int_t    board;
    Int_t    flags;
    Int_t    npulse;
    Int_t    nspe;
    Long64_t timeStamp;
    Double_t energy;
    Double_t wsum;
    Double_t qsum;
    Double_t q900;
    std::vector<TPulse> hits;
		ClassDef(TBaconEvent,1)
};
#endif

