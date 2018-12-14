#include "genPulses.hh"

// pmt signal function
// conv of gaussian and two decaying expontials
// credit to M. Gold UNM
class MyFunction {
  public:
    //copied from root TH1 website..no idea what operator() does
    double operator()(Double_t *x, Double_t *par){
      double xx = x[0]-par[4];
      //cout<<"xx = "<<xx<<endl;
      double y1=1./sqrt(2)/par[0]*(xx-par[0]*par[0]/par[1]);
      //cout<<"y1 = "<<y1<<endl;
      double b1 = par[4]/par[1]+0.5*pow(par[0]/par[1],2.);
      //cout<<"b1 = "<<b1<<endl;
      //double c1=0.5/par[1]*exp(b1)*exp((-x[0])/par[1])*TMath::Erfc(-y1);//
      //cout<<"c1 = "<<c1<<endl;
      //cout<<"\t par[1] "<<par[1]<<" exp(b1) "<<exp(b1)<<" exp((-x[0])/par[1]) "<<exp((-x[0])/par[1])<<" TMath::Erfc(-y1) "<<TMath::Erfc(-y1)<<endl;
      double c1=0.5*exp((-log(par[1])*par[1]+b1*par[1]-x[0])/par[1])*TMath::Erfc(-y1);//
      //cout<<"c1 = "<<c1<<endl;
      //cout<<"\t 0.5* "<<par[1]<<" exp((-log(par[1])*par[1]+b1*par[1]-x[0])/par[1]) "<<exp((-log(par[1])*par[1]+b1*par[1]-x[0])/par[1])<<" TMath::Erfc(-y1) "<<TMath::Erfc(-y1)<<endl;
      double y2=1./sqrt(2)/par[0]*(xx-par[0]*par[0]/par[2]);
      //cout<<"y2 = "<<y2<<endl;
      double b2 = par[4]/par[2]+0.5*pow(par[0]/par[2],2.);
      //cout<<"b2 = "<<b2<<endl;
      double c2=0.5/par[2]*exp(b2)*exp((-x[0])/par[2])*TMath::Erfc(-y2);//
      //cout<<"c2 = "<<c2<<endl;
      return par[5]*(c1*par[3] +  c2*(1.-par[3]) );
    }
};



genPulses::genPulses(Int_t maxEvents)
{
  Int_t nEvents = maxEvents;
  TDatime time;
  //time 12:36:26 133626
  ///date 24/12/1997 19971224
  Double_t gaussMean = 0,gaussSigma = .0009;
  Int_t Nphotons = 1;
  TString outFileName = TString("rootData/simEvents_")+TString(to_string(time.GetDate()))+Form("_%i",nEvents)+TString(".root");
  //TString outFileName = "simEvents.100.100photons.root";
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s npulses %i gaussMean %f gaussSigma %f\n",outFileName.Data(),Nphotons,gaussMean,gaussSigma);
  
 
  TTree * simTree = new TTree("pmtTree","pmtTree");
  pmtSimulation = new TPmtSimulation();
  simTree->Branch("pmtSimulation",&pmtSimulation);
  pmtEvent = new TPmtEvent();
  simTree->Branch("pmtEvent", &pmtEvent);

  TRandom2 rand;
  rand.SetSeed(time.GetTime());
  MyFunction  fObj;
  double s=4.e-9;double t1=2.e-9;double t2=8e-9;double t12=1.3;double mean = 100e-9;double amp = 5e-11;
  

  for(int k = 0; k < nEvents; k++){

    Int_t nbins = 10000;//0000;
    TH1F* hSignal = new TH1F(TString("response_")+to_string(k),TString("response_")+to_string(k),nbins,0,4e-6);
    
    std::vector<Double_t> pulseTimes = PulseStartTime(Nphotons,k+time.GetTime());
    std::vector<Double_t> sig,time;
    
    pmtSimulation->startTime = pulseTimes;
    pmtSimulation->Nphotons = Nphotons;
    pmtSimulation->sigma = s;
    pmtSimulation->tau1 = t1;
    pmtSimulation->tau2 = t2;
    pmtSimulation->ratio12 = t12;
    pmtSimulation->event = k;

  
    for(unsigned i = 0; i < pulseTimes.size();i++){
      mean = pulseTimes[i];
      TF1 * f1 = new TF1("fbaseline",fObj,mean-5.*s, mean+10*t2,6);
      //TF1 * f1 = new TF1("f",fObj,0,2.e-6,5);
      f1->SetLineColor(kBlue);
      f1->SetParameters(s,t1,t2,t12,mean,amp);
      f1->SetParName(0,"sigma");
      f1->SetParName(1,"tau1");
      f1->SetParName(2,"tau2");
      f1->SetParName(3,"ratio12");
      f1->SetParName(4,"mean");
      f1->SetParName(5,"amplitude");
      //cout<<f1->GetParName(0)<<" = "<<f1->GetParameter(0)<<" "<<f1->GetParName(1)<<" = "<<f1->GetParameter(1)<<" "<<f1->GetParName(2)<<" = "<<f1->GetParameter(2)<<" "<<f1->GetParName(3)<<" = "<<f1->GetParameter(3)<<" "<<f1->GetParName(4)<<" = "<<f1->GetParameter(4)<<" "<<endl;
      f1->Draw("lp");
      //cout<<f1->Eval(mean)<<endl;
      //printf("conv.C(%.2f,%.2f,%.2f,%.2f,%.2f)\n",s,t1,t2,t12,mean);
      //printf(" sigma %f tau1 %f tau2 %f fraction 1 to 2 %f mean %f integral %f \n",s,t1,t2,t12,mean,f1->Integral(mean-5.*s,mean+100*t2));
      outFile->cd();
      //for(int i = 0; i < nbins;i++) hSignal->SetBinContent(i+1,f1->Eval(hSignal->GetBinCenter(i)) );
      Double_t norm = 0;
      Int_t startBin = hSignal->FindBin(mean-5.*s),stopBin = hSignal->FindBin(mean+10*t2);
      for(int j = startBin; j < stopBin;j++) {
        //cout<<"norm "<<norm<<", time "<<hSignal->GetBinCenter(j+1)<<" f1 = "<< f1->Eval(hSignal->GetBinCenter(j+1))<<endl;
        norm += f1->Eval(hSignal->GetBinCenter(j+1));
      }
      for(int j= startBin; j < stopBin;j++){
        //BoxMuller transform
        //hSignal->SetBinContent(j+1,hSignal->GetBinContent(j+1)+f1->Eval(hSignal->GetBinCenter(j+1))/norm);
        hSignal->SetBinContent(j+1,hSignal->GetBinContent(j+1)+f1->Eval(hSignal->GetBinCenter(j+1)));
      }
    }
    TF1 * fBaseSag = new TF1("Baseline_Sag","[0]*TMath::Exp(-x/[1])*TMath::Sin((x-[2])/(2*[1]))",300e-9,4e-6);
    fBaseSag->SetParameter(0,-hSignal->GetMaximum()/2);
    fBaseSag->SetParameter(1,1e-6);
    fBaseSag->SetParameter(2,0);
    for(int i = 0; i <hSignal->GetNbinsX();i++){
      Double_t noise = gaussMean+gaussSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
      //hSignal->SetBinContent(i+1,hSignal->GetBinContent(i+1));
      //hSignal->SetBinContent(i+1,noise+hSignal->GetBinContent(i+1));
      if(hSignal->GetBinCenter(i+1) > 300e-9)
        hSignal->SetBinContent(i+1,noise+hSignal->GetBinContent(i+1)+fBaseSag->Eval(hSignal->GetBinCenter(i+1)));
      else
        hSignal->SetBinContent(i+1,noise+hSignal->GetBinContent(i+1));
      sig.push_back(-noise-hSignal->GetBinContent(i+1));
      time.push_back(hSignal->GetBinCenter(i+1));
    }
    pmtEvent->volt1 = sig;
    pmtEvent->time = time;
    simTree->Fill();
    pmtSimulation->clear();
  hSignal->Draw();
  }
  outFile->Write();
  cout<<"end of genPulses "<<outFileName<<" events " << simTree->GetEntries() << endl;
}

std::vector<Double_t> genPulses::PulseStartTime(Int_t nPhotons,Double_t randSeed){
  std::vector<Double_t> pulseTimes;
  Double_t startTime = 0,stopTime = 4e-6,shiftTime = 300e-9;
  TF1 *fScintDist = new TF1("scintDist","expo(0)+expo(2)",startTime,stopTime);
  fScintDist->SetParameter(0,TMath::Log(3.));
  fScintDist->SetParameter(1,-1./6e-9);
  fScintDist->SetParameter(2,TMath::Log(1));
  fScintDist->SetParameter(3,-1./1e-6);
  Double_t max = fScintDist->GetMaximum(startTime,stopTime);
  Double_t min = fScintDist->Eval(stopTime);
  Double_t norm = fScintDist->Integral(startTime,stopTime);
  //TF1 *fBound = new TF1("BoundingFunction","expo(0)+expo(2)",startTime,stopTime);
  TF1 *fBound = new TF1("BoundingFunction","expo(0)",startTime,stopTime);
  fBound->SetParameter(0,TMath::Log(2));
  fBound->SetParameter(1,-1./1e-6);
  //fBound->Draw();
  //fScintDist->Draw("same");
  //fScintDist->SetParameter(4,1./norm);
  TH1F * hTest = new TH1F("test","test",1000,startTime,stopTime);
  TRandom2 rand;
  rand.SetSeed(randSeed);
  Int_t counter = 0,failCounter = 0;

  while(counter < nPhotons){
    failCounter++;
    if(failCounter > 1e7) break;
    //Double_t x = rand.Rndm()*stopTime;
    Double_t x = -1e-6*TMath::Log(1-rand.Rndm());
    //max = fBound->Eval(x);
    Double_t fX = fScintDist->Eval(x);
    Double_t hX = fBound->Eval(x);
    Double_t u = rand.Rndm();
    //cout<<"x "<<x<<", fX "<<fX<<", hX "<<hX<<", u "<<u<<", u*hX "<<u*hX<<endl;
    if(u*hX > fX) continue;
    //Double_t y = rand.Rndm()*(max-min)+min;
    if(x+shiftTime > stopTime) continue;
    x += shiftTime;
    hTest->Fill(x);
    pulseTimes.push_back(x);
    counter++;
  } 
  //hTest->Fit(fScintDist);
  //hTest->Draw();
  hTest->Write();
  delete hTest;
  return pulseTimes;
}

