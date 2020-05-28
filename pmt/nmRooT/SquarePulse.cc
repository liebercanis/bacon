#include <TROOT.h>
#include <TVirtualFFT.h>
#include <TChain.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFormula.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>    // std::sort

std::vector<Double_t> LowPassFrequency(std::vector<Double_t> input, Double_t cutoff, Double_t sampleRate){

  std::vector<Double_t> output;
  output.resize(input.size());
  output[0] = input[0];
  Double_t RC = 1.0/(cutoff*2*3.1415);
  Double_t dt = sampleRate;
  Double_t alpha = dt/(RC+dt);
  for(int i = 1; i< input.size();i++){
    output[i] = output[i-1] +(alpha*(input[i]-output[i-1]));  
  }
  return output;
}

void SquarePulse(){
  Double_t gausMean = 0,gausSigma = .1;
  Int_t Npulses = 1e3;//1e5;
  Double_t deltaT = 1e-9;
  Int_t nBins =100,pulseWidth = 1, signalStart = 50, noiseStart = 10;
  TString outFileName = "SquarePulseOut.root";
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  TRandom2 rand;
  TDatime time;
  rand.SetSeed(time.GetTime());
  
  for(int k = 36;k<=36;k++){
    cout<<"Generating "<<Npulses<<" pulses of width "<<k<<" and noise sigma of "<<gausSigma<<endl;
    pulseWidth = k;
    signalStart = noiseStart+pulseWidth+10;
    if(signalStart > nBins) cout<<"signal Start is greater than nBins"<<endl;
    TH1D * hSig = new TH1D(TString("SignalArea_PulseWidth-")+to_string(k),TString("SignalArea_PulseWidth-")+to_string(k),nBins,pulseWidth*0.5,pulseWidth*1.5);
    TH1D * hNoise = new TH1D(TString("NoiseArea_IntragrationWindow-")+to_string(k),TString("NoiseArea_IntragrationWindow-")+to_string(k),nBins,-20*gausSigma,20*gausSigma);
    TH1D * hSigFilter = new TH1D(TString("SignalFilter-")+to_string(k),TString("SignalFilter-")+to_string(k),nBins,pulseWidth*0.5,pulseWidth*1.5);
    TH1D * hNoiseFilter = new TH1D(TString("NoiseFilter-")+to_string(k),TString("NoiseFilter-")+to_string(k),nBins,-20*gausSigma,20*gausSigma);

    for(int i =0; i < Npulses ;i++){
      if(i%(int (Npulses*.1)) == 0)cout<<"Processed "<<i<<" pulses"<<endl;
      std::vector<Double_t> sig;
      TH1D * hPulse = new TH1D(TString("SquarePulse")+to_string(i)+to_string(k),TString("SquarePulse")+to_string(i)+to_string(k),nBins,0,nBins*deltaT);
      TH1D * hPulseFilter = new TH1D(TString("SquarePulseFilter")+to_string(i)+to_string(k),TString("SquarePulseFilter")+to_string(i)+to_string(k),nBins,0,nBins*deltaT);
      Double_t sumSig = 0, sumNoise = 0;
      for(int j = 0; j < nBins;j++){
        Double_t noise = gausMean+gausSigma*sqrt(-2.0*log(rand.Rndm()))*cos(2*TMath::Pi()*rand.Rndm()) ;
        if(j >=signalStart && j< signalStart+pulseWidth){
          //pulse height of 1
          hPulse->SetBinContent(j,1);
        }
        hPulse->SetBinContent(j,hPulse->GetBinContent(j)+noise);
        sig.push_back(hPulse->GetBinContent(j)+noise);
        if(j >= noiseStart && j < noiseStart + pulseWidth) sumNoise += hPulse->GetBinContent(j);
        if(j >= signalStart && j < signalStart + pulseWidth) sumSig += hPulse->GetBinContent(j);
      }
      hSig->Fill(sumSig);
      hNoise->Fill(sumNoise);

      sig = LowPassFrequency(sig,50e6,deltaT);
      sumNoise = 0; sumSig = 0;
      for(int j = 0;j< sig.size();j++){
        if(j >= noiseStart && j < noiseStart + pulseWidth) sumNoise += sig[j];
        if(j >= signalStart && j < signalStart + pulseWidth) sumSig += sig[j];
        hPulseFilter->SetBinContent(j+1,sig[j]);
      }
      hSigFilter->Fill(sumSig);
      hNoiseFilter->Fill(sumNoise);

      if(i > 10){
        delete hPulse;
        delete hPulseFilter;
      }
      //hPulse->Draw();
    }
  }
outFile->Write();
}

