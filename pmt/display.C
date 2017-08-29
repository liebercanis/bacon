////////////////////////////////////////////////////////
void display(Int_t  ievent=0, Int_t irun=0)
{
  enum {NPMT=2};
  enum {MAXSAMPLES=10000};
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/baconRun_%i.root",irun);
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  // get pmtTree from file 
  TTree* pmtTree;
  fin->GetObject("pmtTree",pmtTree);
  Long64_t nentries = pmtTree->GetEntries();
  cout << " number of entries is " << nentries << endl;

  // set up memory for reading
  TPmtEvent* pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  
  TH1D* hSamples[NPMT];

  TString hname, htitle;
  for(UInt_t ipmt=0; ipmt<NPMT; ++ipmt) {
    hname.Form("Samples_ch%u_run%u_ev%u",ipmt,irun,ievent);
    htitle.Form("Samples  channel %u run %u event %u ",ipmt,irun,ievent);
    hSamples[ipmt] = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
    hSamples[ipmt]->SetXTitle(" sample number ");
    hSamples[ipmt]->GetYaxis()->SetRangeUser(-.5,1.0);
  }
  hSamples[0]->SetLineColor(kBlack); 
  // get the event
  pmtTree->GetEntry(ievent);
  double offset=0.5;
  double expand=100.0;
  //fill histograms 
  for(unsigned ibin=0; ibin< MAXSAMPLES; ++ibin ) { 
    hSamples[0]->SetBinContent(ibin,offset+pmtEvent->volt1[ibin]);
    hSamples[1]->SetBinContent(ibin,expand*pmtEvent->volt2[ibin]);
  }
  cout << "... event " << pmtEvent->event << endl;
  TString canName;
  canName.Form("Run%uEvent%u",irun,ievent);
  gStyle->SetOptStat(0);
  TCanvas *can = new TCanvas(canName,canName);
  hSamples[0]->Draw();
  hSamples[1]->Draw("same");
}

