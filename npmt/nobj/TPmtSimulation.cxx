#include "TPmtSimulation.hxx"
ClassImp(TPmtSimulation)

TPmtSimulation::TPmtSimulation(): TNamed("TPmtSimulation","TPmtSimulation")
{
  clear();
}

//TPmtSimulation::~TPmtSimulation(){}

void TPmtSimulation::clear()
{
  isSim=false;
  event=0;
  sigma=0;
  tau1=0;
  tau2=0;
  ratio12=0;
  Nphotons.clear();
  pmtNum=0;
  startTime.clear();
  peakTime.clear();
  q.clear();
 }
void TPmtSimulation::resize()
{
  startTime.resize(pmtNum);
  peakTime.resize(pmtNum);
  q.resize(pmtNum);
  Nphotons.resize(pmtNum);
  for(int i = 0;i < pmtNum; i++){
    startTime[i].clear();
    peakTime[i].clear();
    q[i].clear();
  }
}
