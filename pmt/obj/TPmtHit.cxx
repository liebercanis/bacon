#include "TPmtHit.hxx"
ClassImp(TPmtHit)

TPmtHit::TPmtHit(): TNamed("TPmtHit","TPmtHit")
{
  clear();
}

//TPmtHit::~TPmtHit(){}

void TPmtHit::clear()
{
  startTime=0;
  peakWidth=0;
  qpeak=0;
  peakt=0;
  peakMaxTime=0;
  peakBin=0;
  qsum=0;
}

