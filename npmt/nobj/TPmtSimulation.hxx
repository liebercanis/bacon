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
    void resize();
		// data elements
    Bool_t   isSim;
    Int_t    event;
    Double_t sigma;
    Double_t tau1;
    Double_t tau2;
    Double_t ratio12;
    std::vector<Int_t> Nphotons;
    Int_t pmtNum;
    std::vector<std::vector<Double_t>> startTime;
    std::vector<std::vector<Double_t>> peakTime;
    std::vector<std::vector<Double_t>> q;
		ClassDef(TPmtSimulation,2)
};
#endif
