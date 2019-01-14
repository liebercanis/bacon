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



genPulses::genPulses(Int_t maxEvents)
{
  Double_t tau3 = 1.0E-6; // triplet lifetime
  Int_t nEvents = maxEvents;
  TDatime time;
  rand.SetSeed(time.GetTime());
  Double_t gaussMean = 0,gaussSigma = .0009;
  Int_t Nphotons = 10;
  TString outFileName = TString("rootData/simEvents_")+TString(to_string(time.GetDate()))+Form("_%i",nEvents)+TString(".root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s npulses %i gaussMean %f gaussSigma %f\n",outFileName.Data(),Nphotons,gaussMean,gaussSigma);

  
  // setup output tree  
  TTree * simTree = new TTree("pmtTree","pmtTree");
  pmtSimulation = new TPmtSimulation();
  simTree->Branch("pmtSimulation",&pmtSimulation);
  pmtEvent = new TPmtEvent();
  simTree->Branch("pmtEvent", &pmtEvent);

  // pulse function
  PulseFunction  PulseFunc;
  double s=4.e-9;double t1=2.e-9;double t2=8e-9;double t12=1.3;double mean = 100e-9;double amp = 5e-11;
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
  fPulse->Print();
  //TCanvas *cbackground = new TCanvas("pulse-shape","pulse-shape");
  //fPulse->Draw("lp");
  outFile->Append(fPulse);

  startTime = 0,stopTime = 4e-6,shiftTime = 300e-9;
 
  //scint function
  fScintDist = new TF1("scintDist","expo(0)+expo(2)",startTime,stopTime);
  fScintDist->SetParameter(0,TMath::Log(3.));
  fScintDist->SetParameter(1,-1./6e-9);
  fScintDist->SetParameter(2,TMath::Log(1));
  fScintDist->SetParameter(3,-1./1e-6);
  //TCanvas *csint = new TCanvas("scint-dist","scint-dist");
  //fScintDist->Draw("lp");
  outFile->Append(fScintDist);
  

  // bounding function 
  fBound = new TF1("BoundingFunction","expo(0)",startTime,stopTime);
  fBound->SetParameter(0,TMath::Log(2));
  fBound->SetParameter(1,-1./1e-6);
  //TCanvas *cbound = new TCanvas("bound-func","bound-func");
  //fBound->Draw("lp");
  outFile->Append(fBound);


  fBaseSag = new TF1("Baseline_Sag","[0]*TMath::Exp(-x/[1])*(1-0.5*TMath::Sin(4*(x-[2])/(2*[1])))",300e-9,4e-6);
  // base sagging
  fBaseSag->SetParameter(0,gaussSigma*10);
  fBaseSag->SetParameter(0,1);
  fBaseSag->SetParameter(1,1e-6);
  fBaseSag->SetParameter(2,0);
  //TCanvas *csag = new TCanvas("sag-func","sag-func");
  //fBaseSag->Draw("lp");
  outFile->Append(fBaseSag);
  

  hTest = new TH1D("test","test",1000,startTime,stopTime);
  Int_t nbins = 10000;//0000;
  hSignal = new TH1D("signal","signal",nbins,0,4e-6);
  hTime = new TH1D("time","pulse time",nbins,0,4e-6);
  hWave = new TH1D("wave","wave",nbins,0,4e-6);
  hNoise = new TH1D("noise","noise",1000,-10*gaussSigma,10*gaussSigma);

  

  
  // loop over events
  for(int k = 0; k < nEvents; k++){
    
    //generate pulse start times
    std::vector<Double_t> pulseTimes = PulseStartTime(k,Nphotons,tau3);
    std::vector<Double_t> sig,time;
    
    pmtSimulation->startTime = pulseTimes;
    pmtSimulation->Nphotons = Nphotons;
    pmtSimulation->sigma = s;
    pmtSimulation->tau1 = t1;
    pmtSimulation->tau2 = t2;
    pmtSimulation->ratio12 = t12;
    pmtSimulation->event = k;


    TH1D* hSignal1 = (TH1D*) hSignal->Clone(Form("SignalEv%i",k));
    TH1D* hWave1 = (TH1D*) hWave->Clone(Form("WaveformEv%i",k));
    // loop over pulses
    for(unsigned i = 0; i < pulseTimes.size();i++){
      mean = pulseTimes[i];
      hTime->Fill(mean);
      fPulse->SetParameter(4,mean);
      //Double_t norm = 0;
      Int_t startBin = hSignal1->FindBin(mean-5.*s),stopBin = hSignal1->FindBin(mean+10*t2);
      //for(int j = startBin; j < stopBin;j++) norm += fPulse->Eval(hSignal1->GetBinCenter(j+1));
      
      for(int j= startBin; j < stopBin;j++){
        //BoxMuller transform
        hSignal1->SetBinContent(j+1,hSignal1->GetBinContent(j+1)+fPulse->Eval(hSignal1->GetBinCenter(j+1)));
        hWave1->SetBinContent(j+1,hSignal1->GetBinContent(j+1)+fPulse->Eval(hSignal1->GetBinCenter(j+1)));
      }
    }

    // set baseline sagging 
    if(-hSignal1->GetMaximum() != 0) fBaseSag->SetParameter(0,-hSignal->GetMaximum()/2);
    else fBaseSag->SetParameter(0,gaussSigma*10);
    fBaseSag->SetParameter(1,1e-6);
    fBaseSag->SetParameter(2,0);

    // add noise
    for(int i = 0; i <hWave1->GetNbinsX();i++){
      Double_t noise = gaussMean+gaussSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
      hNoise->Fill(noise);
      hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1));
      if(hWave1->GetBinCenter(i+1) > 300e-9)
       hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1)+fBaseSag->Eval(hWave1->GetBinCenter(i+1)));
      else
        hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1));
      // make vectors for event
      sig.push_back(-noise-hWave1->GetBinContent(i+1));
      time.push_back(hWave1->GetBinCenter(i+1));
    }
    pmtEvent->volt1 = sig;
    pmtEvent->time = time;
    simTree->Fill();
    pmtSimulation->clear();
  }
  fPulse->SetParameter(4,100e-9);
  fBaseSag->SetParameter(0,gaussSigma*10); 
  outFile->Write();
  cout<<"end of genPulses "<<outFileName<<" events " << simTree->GetEntries() << endl;
}

std::vector<Double_t> genPulses::PulseStartTime(Int_t event, Int_t nPhotons, Double_t tau3){
  std::vector<Double_t> pulseTimes;
  /*
  Double_t max = fScintDist->GetMaximum(startTime,stopTime);
  Double_t min = fScintDist->Eval(stopTime);
  Double_t norm = fScintDist->Integral(startTime,stopTime);
  */

  TH1D * hTestEv = (TH1D*) hTest->Clone(Form("test-Ev%i",event));

  for(Int_t c = 0; c<nPhotons; ++c) {
    Double_t x = -tau3*TMath::Log(1-rand.Rndm());
    hTestEv->Fill(x);
    pulseTimes.push_back(x);
  } 
  //hTest->Fit(fScintDist);
  //hTest->Draw();
  return pulseTimes;
}

