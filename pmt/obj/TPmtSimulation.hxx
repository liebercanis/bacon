/**
** NM, November 2018 
**/
#ifndef TPMTSIMULATION_DEFINED
#define TPMTSIMULATION_DEFINED
#include <iostream>
#include <string>
#include <TNamed.h>
#include <vector>

using namespace std;

// class to store info for the event

class TPmtSimulation: public TNamed {
	public:
		TPmtSimulation();
    //		~TPmtSimulation();

		void clear();
		// data elements
    Int_t    event;
    UShort_t  gpsYear;
    UShort_t  gpsDay;
    UInt_t    gpsSec;
    UInt_t   gpsNs;
    Double_t sigma;
    Double_t tau1;
    Double_t tau2;
    Double_t ratio12;
    Int_t Nphotons;
		std::vector<Double_t> time;		 
		std::vector<Double_t> volt;
    std::vector<Double_t> startTime;
		ClassDef(TPmtSimulation,1)
};
#endif

