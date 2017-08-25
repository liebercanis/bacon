////////////////////////////////////////////////////////
#include "anaRun.hh"

anaRun::anaRun(Int_t irun)
{
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/baconRun_%i.root",irun);
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }

  // get pmtTree from file 
  fin->GetObject("pmtTree",pmtTree);
  Long64_t nentries = pmtTree->GetEntries();
  cout << " number of entries is " << nentries << endl;

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  
  // open file for histograms
  TString outFileName; outFileName.Form("baconRunAna_%i.root",irun);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->cd();
  printf(" opening output file %s \n",outFileName.Data());

  TString hname, htitle;
  for(UInt_t ipmt=0; ipmt<NPMT; ++ipmt) {
    hname.Form("Samples_pmt%u",ipmt);
    htitle.Form("Samples  pmt%u",ipmt);
    hSamples[ipmt] = new TH1D(hname,htitle,MAXSAMPLES,0,MAXSAMPLES);
    hSamples[ipmt]->SetXTitle(" sample number ");
  }
  
  // initialize fft 
  nFFTSize = Int_t(MAXSAMPLES);
  fFFT = TVirtualFFT::FFT(1, &nFFTSize, "R2C M K");
  fInverseFFT = TVirtualFFT::FFT(1, &nFFTSize, "C2R M K");
  
  // loop over entries 
  for (Long64_t ientry=0; ientry<nentries; ientry++) {
    pmtTree->GetEntry(ientry);
    if(ientry==0) { hFFT[0]=FFTFilter(0); hFFT[1]=FFTFilter(1);}
    cout << "... event " << pmtEvent->event << endl;
    if(ientry>1) break;
  }

  outfile->Write();
}

TH1D* anaRun::FFTFilter(Int_t ipmt)
{
  if(ipmt==0) printf(" called FFTFilter pmt %i volts size %zu  \n",ipmt, pmtEvent->volt1.size()); 
  else printf(" called FFTFilter pmt %i volts size %zu  \n",ipmt, pmtEvent->volt2.size()); 
  for(int is =0; is<nFFTSize; ++is) {
    if(ipmt==0) fFFT->SetPoint(is, pmtEvent->volt1[is]);
    else fFFT->SetPoint(is, pmtEvent->volt2[is]);
  }

  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
  }

  fFFT->Transform();
  TH1D* hfft = new TH1D(Form("FFTPmt%i",ipmt),Form("FFT PMT %i",ipmt),nFFTSize/2,0,nFFTSize/2);

  // fill samples FFT histogram && elec response in time domain
  printf(" created %s %s \n",hfft->GetName(),hfft->GetTitle());
  // skip first bin which is pedestal
  for (int i = 1; i<nFFTSize/2; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im);
    hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
  } 
  return hfft;      
}


