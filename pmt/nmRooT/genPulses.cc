
#include "genPulses.hh"

// pmt signal function
// conv of gaussian and two decaying expontials
// credit to M. Gold UNM
// /
class PulseFunction {
  public:
    //copied from root TH1 website with operator() usage syntax
    double operator()(Double_t *x, Double_t *par){
      double xx = x[0]-par[4];
      double y1=1./sqrt(2)/par[0]*(xx-par[0]*par[0]/par[1]);
      double b1 = par[4]/par[1]+0.5*pow(par[0]/par[1],2.);
      double c1=0.5*exp((-log(par[1])*par[1]+b1*par[1]-x[0])/par[1])*TMath::Erfc(-y1);//
      double y2=1./sqrt(2)/par[0]*(xx-par[0]*par[0]/par[2]);
      double b2 = par[4]/par[2]+0.5*pow(par[0]/par[2],2.);
      double c2=0.5/par[2]*exp(b2)*exp((-x[0])/par[2])*TMath::Erfc(-y2);//
      return par[5]*(c1*par[3] +  c2*(1.-par[3]) );
    }
};



genPulses::genPulses(Int_t maxEvents,Int_t meanPhotons = 0,Int_t dsNum = 2)
{
  Double_t tau3 = 1.0E-6; // triplet lifetime
  Double_t eventTime = 10.0E-6;
  Int_t nEvents = maxEvents;
  if(meanPhotons == 0){
    cout<<"please set meanPhotons/event in intial excutable"<<endl;
  }
  Int_t nbins = 10000, nPMTs =0;
  startTime = 0.0,stopTime = 10e-6,shiftTime = 1e-6;
  double binWidth = 1e-9;
  bool debug = false;//false;///true;
  TDatime time;
  rand.SetSeed(time.GetTime());
  //Double_t gaussMean = 4.e-4,gaussSigma = 1.15e-3;//4.89e-5;//.0009;
  Double_t gaussMean = 0,gaussSigma = 0;
  if(dsNum == 1){
    gaussMean = 4.e-4;
    gaussSigma = 1.15e-3;
    nPMTs = 1;
  }
  else if(dsNum == 2){
    gaussMean = -3.e-4;
    gaussSigma = 3.2e-4;
    nPMTs = 2;
  }
  TString outFileName = "";
  if(dsNum == 1){
    outFileName = TString("rootData/DS1/simEvents_DS_")+to_string(dsNum)+Form("_nPhotons_%i",meanPhotons)+Form("_nEvents_%i",maxEvents)+TString(".root");
  }
  else if(dsNum == 2){
    outFileName = TString("rootData/DS2/simEvents_DS_")+to_string(dsNum)+Form("_nPhotons_%i",meanPhotons)+Form("_nEvents_%i",maxEvents)+TString(".root"); 
  }
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s npulses %i gaussMean %f gaussSigma %f\n",outFileName.Data(),meanPhotons,gaussMean,gaussSigma);

 
  // setup output tree  
  TTree * simTree = new TTree("pmtTree","pmtTree");
  pmtSimulation = new TPmtSimulation();
  simTree->Branch("pmtSimulation",&pmtSimulation);
  pmtEvent = new TPmtEvent();
  simTree->Branch("pmtEvent", &pmtEvent);
  if(dsNum == 1) pmtSimulation->pmtNum =1;
  else if(dsNum == 2) pmtSimulation->pmtNum=2;
  pmtSimulation->resize();

  // pulse function
  PulseFunction  PulseFunc;
  //The under/over shoot is tuned using t12, if t12 > 1 under shoot, 
  double s=4.e-9;double t1=2.e-9;double t2=8e-9;double t12=.8;double mean = 100e-9;double amp = 5.e-11;
  //double pulseMax = mean+10*t2;
  TF1 * fPulse = new TF1("fpulse",PulseFunc,mean-5.*s, mean+10*t2,6);
  //TF1 * fPulse = new TF1("f",fObj,0,2.e-6,5);
  fPulse->SetLineColor(kBlue);
  fPulse->SetParameters(s,t1,t2,t12,mean,amp);
  fPulse->SetParName(0,"sigma");
  fPulse->SetParName(1,"tau1");
  fPulse->SetParName(2,"tau2");
  fPulse->SetParName(3,"ratio12");
  fPulse->SetParName(4,"mean");
  fPulse->SetParName(5,"amplitude");
  double pulseMax = fPulse->GetMaximum();
  double pulseMaxTime = fPulse->GetX(pulseMax);
  //cout<<"mean "<<mean<<", pulse max "<<fPulse->GetX(fPulse->GetMaximum())<<endl;

  //TCanvas *cbackground = new TCanvas("pulse-shape","pulse-shape");
  //fPulse->Draw("lp");
  outFile->Append(fPulse);

  //TH1D* hPulse = (TH1D*) fPulse->GetHistogram();
  TFile* fSPE;
  TH1D *hTemp[nPMTs];
  if(dsNum == 1){
    fSPE = new TFile("processedData/DS1/SPE.root");
    fSPE->GetObject("SPE_0",hTemp[0]);
  }
  else if(dsNum == 2){
    fSPE = new TFile("processedData/DS2/SPE.root");
    fSPE->GetObject("SPE_0",hTemp[0]);
    fSPE->GetObject("SPE_1",hTemp[1]);
  }
  TH1D* hPulse[nPMTs];
  if(dsNum == 1)
    hPulse [0] = new TH1D("SPE","SPE",hTemp[0]->GetNbinsX(),0,hTemp[0]->GetNbinsX()*binWidth);
  if(dsNum == 2){
    hPulse [0] = new TH1D("SPE0","SPE0",hTemp[0]->GetNbinsX(),0,hTemp[0]->GetNbinsX()*binWidth);
    hPulse [1] = new TH1D("SPE1","SPE1",hTemp[1]->GetNbinsX(),0,hTemp[1]->GetNbinsX()*binWidth);
  }
  std::vector<Double_t> qnorm,peakTime;
  qnorm.resize(nPMTs);
  peakTime.resize(nPMTs);
  int oldBin =0,binCounter =0;
  for(int j = 0; j <  nPMTs;j++){
    for(int i = 0; i < hTemp[j]->GetNbinsX();i++){
      double binVal = hTemp[j]->GetBinContent(i+1);
      double binCenter = hTemp[j]->GetBinCenter(i+1);
      int newBin = hPulse[j]->FindBin(binCenter);
      if(oldBin == newBin)binCounter++;
      else binCounter =1;
      hPulse[j]->SetBinContent(newBin,(hPulse[j]->GetBinContent(newBin)+binVal)/binCounter);
      oldBin = newBin;
    }
    qnorm[j] = hPulse[j]->Integral();
    peakTime[j] = hPulse[j]->GetBinCenter(hPulse[j]->GetMaximumBin());
    cout<<peakTime[j]<<endl;
  }
  fPulse->Print();
  if(dsNum == 1)
    printf(" PulseFunction normalization is for pmt0 %f \n",qnorm[0]);
  else if( dsNum == 2)
    printf(" PulseFunction normalization is for pmt0 %f pmt1 %f \n",qnorm[0],qnorm[1]);

  //scint function
  fScintDist = new TF1("scintDist","expo(0)+expo(2)",startTime,stopTime);
  fScintDist->SetParameter(0,TMath::Log(10));
  fScintDist->SetParameter(1,-1./6e-9);
  fScintDist->SetParameter(3,-1./1e-6);
  //TCanvas *csint = new TCanvas("scint-dist","scint-dist");
  //fScintDist->Draw("lp");
  outFile->Append(fScintDist);
  

  // bounding function 
  fBound = new TF1("BoundingFunction","expo(0)",startTime,stopTime);
  fBound->SetParameter(0,TMath::Log(20));
  fBound->SetParameter(1,-1./1e-6);
  //TCanvas *cbound = new TCanvas("bound-func","bound-func");
  //fBound->Draw("lp");
  outFile->Append(fBound);

  ///*
  fBaseSag = new TF1("Baseline_Sag","landau(0)",0,eventTime);
  fBaseSag->SetParameter(0,0.01);
  fBaseSag->SetParameter(1,shiftTime+.75e-6);
  fBaseSag->SetParameter(2, 0.25*shiftTime);
  //*/
 
  outFile->Append(fBaseSag);
 
  // histograms
  outFile->cd(); 
  hSignal = new TH1D("Signal","signal",nbins,0,stopTime);
  hWave = new TH1D("Wave","wave",nbins,0,stopTime);

  hTest = new TH1D("Test","test",1000,0,stopTime*1E6);
  hTestq = new TH1D("Testq","testq",1000,0,stopTime);
  hTime = new TH1D("Time","pulse time",nbins,0,stopTime);
  hNoise = new TH1D("Noise","noise",1000,-10*gaussSigma+gaussMean,10*gaussSigma+gaussMean);

  TH1D* hCharge = new TH1D("Charge","charge",1000,-2,8);

  for(int k = 0; k < nEvents; k++){
    // loop over events
    pmtSimulation->isSim = true;
    pmtSimulation->sigma = s;
    pmtSimulation->tau1 = t1;
    pmtSimulation->tau2 = t2;
    pmtSimulation->ratio12 = t12;
    pmtSimulation->event = k;

    for(int n = 0; n < nPMTs; n++){

      if(k%100==0) printf("... event %i\n",k);

      //generate pulse start times
      std::vector<Double_t> pulseTimes = PulseStartTime(k,meanPhotons,tau3);
      std::vector<Double_t> sig,time;

      pmtSimulation->Nphotons[n] = Int_t(pulseTimes.size());
      if(debug) cout<<"making histograms"<<endl;
      TH1D* hSignal1 = (TH1D*) hSignal->Clone(Form("SignalEv%i_pmt%i",k,n));
      hSignal1->SetTitle(TString("Signal...Number of Pulses Generated") + to_string(pulseTimes.size()));
      TH1D* hWave1 = (TH1D*) hWave->Clone(Form("WaveformEv%i_pmt%i",k,n));
      hWave1->SetTitle(TString("Waveform...Number of Pulses Generated") + to_string(pulseTimes.size()));
      TH1D* hBaselineSag = new TH1D(TString("BaselineSag")+to_string(k)+to_string(n),TString("BaselineSag")+to_string(k)+to_string(n),nbins,0,stopTime);
      hBaselineSag->SetLineColor(2);
      // loop over pulses
      Double_t deltaT = 17e-9;//Peak height is shifted relative to gen time
      if(debug) cout<<"looping over pulse times"<<endl;
      for(unsigned i = 0; i < pulseTimes.size();i++){
        double time = pulseTimes[i];
        hTime->Fill(time);
        Double_t sumq=0;
        double peakSigma;
        if(dsNum == 1) peakSigma = 0.33;
        else if(dsNum == 2) peakSigma = 0.30;
        double peakJitter = 1.+peakSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
        if(debug) cout<<"making hit"<<endl;
        for(int j=1;j<=hPulse[n]->GetNbinsX();++j) {
          double qbin = hPulse[n]->GetBinContent(j);
          double qtime = hPulse[n]->GetBinLowEdge(j);
          //add in peak height fluctations to account for possion statistics off of first anode in PMT
          qbin *= peakJitter;
          Int_t jbin = hSignal1->FindBin(qtime+time-deltaT);
          //Int_t jbin = hSignal1->FindBin(qtime+time);
          Double_t binVal = hSignal1->GetBinContent(jbin);
          hSignal1->SetBinContent(jbin,binVal+qbin);
          hWave1->SetBinContent(jbin,binVal+qbin);
          sumq += qbin;
          //cout<<"bin val "<<qbin<<", qtime "<<qtime<<" jbin "<<jbin<<" hSignal "<<hSignal1->GetBinContent(jbin)<<" time "<<time<<endl;
        }
        if(debug) cout<<"Filling Start times and charge "<<pmtSimulation->startTime[n].size()<<endl;
        //pmtSimulation->startTime.push_back(pulseTimes[i]+mean);
        pmtSimulation->startTime[n].push_back(pulseTimes[i]);
        pmtSimulation->peakTime[n].push_back(peakTime[n]+time-deltaT);
        pmtSimulation->q[n].push_back(sumq);
        Int_t ibin = hTestq->FindBin(time);
        hTestq->SetBinContent( ibin , hTestq->GetBinContent(ibin)+sumq);
        hCharge->Fill(sumq/qnorm[n]);
        if(debug) cout<<"Done Filling startTimes and charge"<<endl;
      }
      //cout<<"Peak time for hSignal "<<hSignal1->GetBinCenter(hSignal1->GetMaximumBin())<<endl;;
         /*
         hSignal1->Draw();
         TLine *line0 = new TLine(pmtSimulation->startTime[n][0],0,pmtSimulation->startTime[n][0],10e-3);
         line0->SetLineColor(2);
         line0->Draw("same");
         Double_t maxBinDebug = pmtSimulation->startTime[n][0]-mean+pulseMaxTime;//hSignal1->GetBinCenter(hSignal1->GetMaximumBin());
         TLine *line1 = new TLine(maxBinDebug,0,maxBinDebug,10e-3);
         line1->SetLineColor(3);
         line1->Draw("same");
         TLine *line2 = new TLine(hSignal1->GetBinCenter(hSignal1->FindBin(maxBinDebug) ),0,hSignal1->GetBinCenter(hSignal1->FindBin(maxBinDebug) ),10e-3);
         line2->SetLineColor(4);
         line2->Draw("same");
         cout<<maxBinDebug<<" "<<pmtSimulation->startTime[n][0]<<endl;
         */
      Double_t maxVal = hWave1->GetMaximum();
      Int_t maxBin = hWave1->GetMaximumBin();
      Double_t maxTime = hWave1->GetBinCenter(maxBin);
      fBaseSag->SetParameter(0,0.5*maxVal);
      fBaseSag->SetParameter(1,1.5*maxTime);
      //fBaseSag->SetParameter(1,maxTime);
      // add noise
      if(debug) cout<<"adding noise"<<endl;
      for(int i = 0; i <hWave1->GetNbinsX();i++){
        Double_t noise = gaussMean+gaussSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
        //uncomment for just noise and no baseline sagging
        hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1));
        hNoise->Fill(noise);
        //add in baseline sag and noise
        //hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1)-fBaseSag->Eval(hWave1->GetBinCenter(i+1)));
        //Make pulse negative
        hWave1->SetBinContent(i+1,-hWave1->GetBinContent(i+1));
        hSignal1->SetBinContent(i+1,-hSignal1->GetBinContent(i+1));
        if(meanPhotons > 10 && dsNum == 1)hBaselineSag->SetBinContent(i+1,fBaseSag->Eval(hWave1->GetBinCenter(i+1)));
        // make vectors for event
        sig.push_back(hWave1->GetBinContent(i+1));
        time.push_back(hWave1->GetBinCenter(i+1));
      }
      //hSignal1->SetMinimum(-0.001);
      hSignal1->SetLineColor(2);
      if(n == 0){
        pmtEvent->volt1 = sig;
        pmtEvent->time = time;
      }
      else if(n == 1)
        pmtEvent->volt2 = sig;

      //hWave1->Draw();
      if(k>100) {
        delete hSignal1;
        delete hWave1;
        delete hBaselineSag;
      }
    }
    simTree->Fill();
    pmtSimulation->clear();
    pmtSimulation->pmtNum = nPMTs;
    pmtSimulation->resize();
  }
  outFile->Write();
  cout<<"end of genPulses "<<outFileName<<" events " << simTree->GetEntries() << endl;
}

std::vector<Double_t> genPulses::PulseStartTime(Int_t event, Double_t meanPhotons, Double_t tau3)
{
  std::vector<Double_t> pulseTimes;

  //TH1D * hTestEv = (TH1D*) hTest->Clone(Form("test-Ev%i",event));
  Int_t counter = 0,failCounter = 0;

  Int_t nPhotons;
  if(meanPhotons == 1)
    nPhotons = 1;
  else
    nPhotons= rand.Poisson(meanPhotons);
  //Int_t nPhotons = 1;
  while(counter < nPhotons){
    failCounter++;
    if(failCounter > 1e7) break;
    Double_t x = -tau3*TMath::Log(1-rand.Rndm());
    Double_t fX = fScintDist->Eval(x);
    Double_t hX = fBound->Eval(x);
    Double_t u = rand.Rndm();
    if(u*hX > fX) continue;
    if(x+shiftTime > stopTime) continue;
    x += shiftTime;
    hTest->Fill(x*1E6);
    pulseTimes.push_back(x);
    counter++;
  } 
  //hTest->Fit(fScintDist);
  //hTest->Draw();
  return pulseTimes;
}
