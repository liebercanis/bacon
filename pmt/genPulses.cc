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
  bool debug=false;
  Double_t tau3 = 1.0E-6; // triplet lifetime
  Double_t eventTime = 10.0E-6;
  Int_t nEvents = maxEvents;
  Double_t meanPhotons=10.0;
  Int_t nbins = 10000;//0000;
  startTime = 1.0,stopTime = 10e-6,shiftTime = 1e-6;
  TDatime time;
  rand.SetSeed(time.GetTime());
  Double_t gaussMean = 0,gaussSigma = .0009;
  TString outFileName = TString("rootData/mySim_")+TString(to_string(time.GetDate()))+Form("_%i",nEvents)+TString(".root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s npulses %f gaussMean %f gaussSigma %f\n",outFileName.Data(),meanPhotons,gaussMean,gaussSigma);

 
  // setup output tree  
  TTree * simTree = new TTree("pmtTree","pmtTree");
  pmtSimulation = new TPmtSimulation();
  simTree->Branch("pmtSimulation",&pmtSimulation);
  pmtEvent = new TPmtEvent();
  simTree->Branch("pmtEvent", &pmtEvent);

  // pulse function
  PulseFunction  PulseFunc;
  double s=4.e-9;double t1=2.e-9;double t2=8e-9;double t12=1.3;double mean = 100e-9;double amp = 5e-11;
  double pulseMax = mean+10*t2;
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
  //TCanvas *cbackground = new TCanvas("pulse-shape","pulse-shape");
  //fPulse->Draw("lp");
  outFile->Append(fPulse);

  //BoxMuller transform
  TH1D* hPulse = (TH1D*) fPulse->GetHistogram();

  fPulse->Print();
  Double_t qnorm = hPulse->Integral();
  printf(" PulseFunction normalization is %f \n",qnorm);


  //scint function
  fScintDist = new TF1("scintDist","expo(0)+expo(2)",startTime,stopTime);
  fScintDist->SetParameter(0,TMath::Log(3.));
  fScintDist->SetParameter(1,-1./6e-9);
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


  //fBaseSag = new TF1("Baseline_Sag","[0]*TMath::Exp(-x/[1])*(1-0.5*TMath::Sin(4*(x-[2])/(2*[1])))",300e-9,4e-6);
  fBaseSag = new TF1("Baseline_Sag","[0]*TMath::Exp(-x/[1])*(1-0.5*TMath::Sin(2*TMath::Pi()*x/[2]))",0,eventTime);
  // base sagging
  fBaseSag->SetParameter(0,gaussSigma*10);
  fBaseSag->SetParameter(1,eventTime);
  fBaseSag->SetParameter(2,eventTime/2.0);
  //TCanvas *csag = new TCanvas("sag-func","sag-func");
  //fBaseSag->Draw("lp");
 
  outFile->Append(fBaseSag);
 
  // histograms
  outFile->cd(); 
  hSignal = new TH1D("Signal","signal",nbins,0,stopTime);
  hWave = new TH1D("Wave","wave",nbins,0,stopTime);

  hTest = new TH1D("Test","test",1000,0,stopTime*1E6);
  hTestq = new TH1D("Testq","testq",1000,0,stopTime);
  hTime = new TH1D("Time","pulse time",nbins,0,stopTime);
  hNoise = new TH1D("Noise","noise",1000,-10*gaussSigma,10*gaussSigma);

  TH1D* hCharge = new TH1D("Charge","charge",1000,-2,8);

  // loop over events
  for(int k = 0; k < nEvents; k++){

    if(k%100==0) printf("... event %i\n",k);
    
    //generate pulse start times
    std::vector<Double_t> pulseTimes = pulseStartTime(k,meanPhotons,tau3);
    std::vector<Double_t> sig,time;

    pmtSimulation->Nphotons = Int_t(pulseTimes.size());
    pmtSimulation->sigma = s;
    pmtSimulation->tau1 = t1;
    pmtSimulation->tau2 = t2;
    pmtSimulation->ratio12 = t12;
    pmtSimulation->event = k;
    TH1D* hSignal1 = (TH1D*) hSignal->Clone(Form("SignalEv%i",k));
    TH1D* hWave1 = (TH1D*) hWave->Clone(Form("WaveformEv%i",k));
    // loop over pulses
    for(unsigned i = 0; i < pulseTimes.size();i++){
      double time = pulseTimes[i];
      hTime->Fill(time);
      // 
      Double_t sumq=0;
      for(int j=1;j<=hPulse->GetNbinsX();++j) {
        double qbin = hPulse->GetBinContent(j);
        double qtime = hPulse->GetBinLowEdge(j);
        Int_t jbin = hSignal1->FindBin(qtime+time);
        hSignal1->SetBinContent(jbin,hSignal1->GetBinContent(jbin)+qbin);
        hWave1->SetBinContent(jbin,hWave1->GetBinContent(jbin)+qbin);
        sumq += qbin;
      }
      pmtSimulation->startTime.push_back(pulseTimes[i]+mean);
      pmtSimulation->q.push_back(sumq);
      Int_t ibin = hTestq->FindBin(time);
      hTestq->SetBinContent( ibin , hTestq->GetBinContent(ibin)+sumq);
      hCharge->Fill(sumq/qnorm);
    }

    if(debug) printf(" int %i pulses %lu integral %f \n",k,pulseTimes.size(),hSignal1->Integral()/.05);
    // set baseline sagging 
    //if(-hSignal1->GetMaximum() != 0) fBaseSag->SetParameter(0,-hSignal->GetMaximum()/2);
    //else fBaseSag->SetParameter(0,gaussSigma*10);
    //fBaseSag->SetParameter(1,1e-6);
    //fBaseSag->SetParameter(2,0);

    // add noise
    for(int i = 0; i <hWave1->GetNbinsX();i++){
      Double_t noise = gaussMean+gaussSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
      hNoise->Fill(noise);
      hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1)+fBaseSag->Eval(hWave1->GetBinCenter(i+1)));
      //hWave1->SetBinContent(i+1,noise+hWave1->GetBinContent(i+1));
      // make vectors for event
      sig.push_back(-hWave1->GetBinContent(i+1));
      time.push_back(hWave1->GetBinCenter(i+1));
    }
    pmtEvent->volt1 = sig;
    pmtEvent->time = time;
    simTree->Fill();
    pmtSimulation->clear();
    if(k>100) {
      delete hSignal1;
      delete hWave1;
    }
  }
  fPulse->SetParameter(4,100e-9);
  fBaseSag->SetParameter(0,gaussSigma*10); 
  outFile->Write();
  cout<<"end of genPulses "<<outFileName<<" events " << simTree->GetEntries() << endl;
}

std::vector<Double_t> genPulses::pulseStartTime(Int_t event, Double_t meanPhotons, Double_t tau3)
{
  std::vector<Double_t> pulseTimes;

  //TH1D * hTestEv = (TH1D*) hTest->Clone(Form("test-Ev%i",event));
  Int_t counter = 0,failCounter = 0;

  Int_t nPhotons = rand.Poisson(meanPhotons);

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

