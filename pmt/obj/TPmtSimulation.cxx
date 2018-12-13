#include "TPmtSimulation.hxx"
ClassImp(TPmtSimulation)

TPmtSimulation::TPmtSimulation(): TNamed("TPmtSimulation","TPmtSimulation")
{
  clear();
}

//TPmtSimulation::~TPmtSimulation(){}

void TPmtSimulation::clear()
{
  event=0;
  gpsYear=0;
  gpsDay=0;
  gpsSec=0;
  gpsNs=0;
  sigma=0;
  tau1=0;
  tau2=0;
  ratio12=0;
  Nphotons=0;
  time.clear();	 
  volt.clear();
  startTime.clear();
 }

