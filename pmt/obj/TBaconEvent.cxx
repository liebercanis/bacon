#include "TBaconEvent.hxx"
ClassImp(TBaconEvent)

TBaconEvent::TBaconEvent(): TNamed("TBaconEvent","TBaconEvent")
{
  clear();
}

//TBaconEvent::~TBaconEvent(){}

void TBaconEvent::clear()
{
  event=0;
  channel=0;
  board=0;
  flags=0;
  npulse=0;
  nspe=0;
  timeStamp=0;
  energy=0;
  wsum=0;
  qsum=0;
  q900=0;
  hits.clear();
}

