#include "TGraph.h"
#include "Riostream.h"
void Integrate_and_Plot() {

  Int_t nfiles = 1;
  TString fileName1 = "ChannelOne.txt";
  TString fileName2 = "ChannelTwo.txt";
  TString dir = gSystem->UnixPathName("/data1/nmcfadde/PMT/data/Photon_Counting_9/");
  //TString dir = gSystem->UnixPathName("/home/nmcfadde/RooT/BACon/PMT/data/DarkRate_4/");
  TString outputFileName = TString("PMT-")+dir.ReplaceAll("/",".")+TString("root");
  //some fancy shit to get the name from the director
  outputFileName.Remove(3,25);
  cout<<outputFileName<<endl;
  dir.ReplaceAll(".","/");
  TFile *outfile = new TFile(outputFileName,"recreate"); 
  Int_t photonBins = 400;
  TH1F* hPhotonCount = new TH1F("Charge of PMT signal - baseline","Charge of PMT signal - baseline",photonBins,1e-10,5e-13);
  TH1F* hT0 = new TH1F("T0_Pulser","Pulser T0 for 1.0 V threshold",100,-1.8e-7,-3.5e-7);
  TH1F* hSignal = new TH1F("Charge of PMT signal + baseline","Charge of PMT signal + baseline",photonBins,1e-10,5e-13);
  TH1F* hNoise = new TH1F("Noise","Integral of noise",photonBins,1e-10,5e-13);
  TH1F* hVoltsMax = new TH1F("Volts_Max","Max Voltage PMT",50,0, -100.e-3);
  TH1F* hBaselineOffSet = new TH1F("Intercept for Baseline Subtraction","Intercept for Baseline Subtraction", photonBins,1e-10,5e-13);
  TH1F* hBaselineSlope = new TH1F("Slope for Baseline Subtraction","Slope for Baseline Subtraction",photonBins,1e-10,5e-13);
  cout<<"Working directory is "<<dir<<endl;
  while(1){
    ifstream in;
    in.open(dir+fileName2 );
    if(in.is_open() ){

      TString name; name.Form("%d",nfiles);
      //10,000 bins for each histogram, labview bins data this way
      TH1F * hPMTSignal = new TH1F("PMTSignal"+name,"",1000,-1.0e-6,1.0e-6);
      TH1F * hPulserSignal = new TH1F("PulseSignal"+name,"",10000,-1.0e-6,1.0e-6);

      std::vector<Float_t> timeTrigger,voltsTrigger;
      bool pulseSwitch = false;
      Int_t nlines = 0;
      Double_t T0;
      while (1) {
        Float_t time,volts;
        in >> time >> volts;
        timeTrigger.push_back(time);
        voltsTrigger.push_back(volts);
        if(volts > 1.&& !pulseSwitch ){
          hT0->Fill(time);
          T0 = time;
          pulseSwitch = true;
        }
        Int_t tBin = hPulserSignal->FindBin(time);
        hPulserSignal->SetBinContent(tBin,hPulserSignal->GetBinContent(tBin)+volts);
        if (!in.good()) break;
        nlines++;
      }

      in.close();

      in.open(dir+fileName1);
      std::vector<Float_t> timePMT,voltsPMT;
      Float_t voltsMax = 9999.;
      while(1){
        Float_t time,volts;
        in>>time>>volts;
        timePMT.push_back(time);
        voltsPMT.push_back(volts);
        Int_t tBin = hPMTSignal->FindBin(time);
        hPMTSignal->SetBinContent(tBin,hPMTSignal->GetBinContent(tBin)+volts);
        if(hPMTSignal->GetBinContent(tBin) < voltsMax) voltsMax = hPMTSignal->GetBinContent(tBin);

        if(!in.good() ) break;
      }
      
      hVoltsMax->Fill(voltsMax);
      
      if(nfiles%500 == 0 && 0){
         TString nameTest; nameTest.Form("%d",nfiles);
        TCanvas* cPhotonCountTrigger = new TCanvas("Trigger"+nameTest,"Trigger"+nameTest);
        hPulserSignal->Draw();
        //TGraph* grTriggerEvent = new TGraph(timeTrigger.size(),&(timeTrigger[0]),&(voltsTrigger[0]));
        //grTriggerEvent->Draw();

        TCanvas* cPhotonCountPMT = new TCanvas("PMTasd"+nameTest,"PMTqwe"+nameTest);
        hPMTSignal->Draw();
        //TGraph* grPMTEvent = new TGraph(timePMT.size(),&(timePMT[0]),&(voltsPMT[0]));
        //grPMTEvent->Draw();
      }
      
      if(pulseSwitch){
        TString nameTest; nameTest.Form("%d",nfiles);
                   
        Int_t T0Bin = hPMTSignal->FindBin(T0);
        
        Double_t TWindow = 300.e-9;//ns window
        Double_t TDelaySignal = 500.e-9;//ns signal delay set by pulser
        Double_t TDelayNoise = 100.e-9;//ns noise delay choosen by Neil to avoid pulser noise
        
        Int_t TFinalBin = hPMTSignal->FindBin(T0+TDelaySignal+TWindow);
        Int_t TIntialBin = hPMTSignal->FindBin(T0+TDelaySignal);
        
        Int_t T1IntialNoiseBin =  hPMTSignal->FindBin(T0+TDelaySignal-TDelayNoise);
        Int_t T1FinalNoiseBin = TIntialBin;
        Int_t T1AverBin = (T1IntialNoiseBin+T1FinalNoiseBin)/2;
        
        Int_t T2IntialNoiseBin = TFinalBin;
        Int_t T2FinalNoiseBin = hPMTSignal->FindBin(T0+TDelaySignal+TWindow+TDelayNoise);
        Int_t T2AverBin = (T2IntialNoiseBin+T2FinalNoiseBin)/2;       
       
        Double_t averageNoise1 = 0;
        for(int i = T1IntialNoiseBin; i < T1FinalNoiseBin; i++){
          averageNoise1 += hPMTSignal->GetBinContent(i);
        }
        averageNoise1 /= (T1FinalNoiseBin-T1IntialNoiseBin);

        Double_t averageNoise2 = 0;
        for(int i = T2IntialNoiseBin; i < T2FinalNoiseBin; i++){
          averageNoise2 += hPMTSignal->GetBinContent(i);
        }
        averageNoise2 /= (T2FinalNoiseBin-T2IntialNoiseBin);
        
        //Noise slope
        Double_t BaselineSlope = (averageNoise1-averageNoise2)/(T1AverBin-T2AverBin);
        Double_t BaselineOffset = averageNoise1-BaselineSlope*T1AverBin;
        //If the slope is calculated both ways and the difference is greater than 1%
        if(std::fabs((BaselineOffset - (averageNoise2-BaselineSlope*T2AverBin))/BaselineOffset) > .01 ){
          cout<<"Something is wrong with your BaselineSlope slope incerpet thing!"<<endl;
          cout<<T1AverBin<<" "<<T2AverBin<<endl;
          cout<<BaselineOffset<<"    "<<averageNoise2-BaselineSlope*T2AverBin<<endl;
          break;
        }
        Double_t integral = 0,baseline = 0,signal = 0;;
        for(int i = TIntialBin; i < TFinalBin; i++){
          //dt*(V_pmt-(dV_noise/dt * t - V0_noise) )
          integral += hPMTSignal->GetBinWidth(i)*(hPMTSignal->GetBinContent(i) - (BaselineSlope*hPMTSignal->GetBinCenter(i)-BaselineOffset) );
          baseline += hPMTSignal->GetBinWidth(i)*(BaselineSlope*hPMTSignal->GetBinCenter(i)-BaselineOffset);
          signal   += hPMTSignal->GetBinWidth(i)*(hPMTSignal->GetBinContent(i) );
           if(nfiles == 10){
             cout<<"Noise Time Bin One start/stop"<<T1IntialNoiseBin<<"/"<<T1FinalNoiseBin<<endl;
             cout<<"Signal Time Bin start/stop"<<TIntialBin<<"/"<<TFinalBin<<endl;
             cout<<"Noise Time Bin Two start/stop"<<T2IntialNoiseBin<<"/"<<T2FinalNoiseBin<<endl;
           }

        }
        hNoise->Fill(baseline/-50.);
        hSignal->Fill(signal/-50.);
        hPhotonCount->Fill(integral/-50.0);
        hBaselineSlope->Fill(BaselineSlope);
        hBaselineOffSet->Fill(BaselineOffset);

        if(nfiles%100 == 0){ 
          //TCanvas* cPhotonCountTrigger = new TCanvas("Trigger"+nameTest,"Trigger"+nameTest);
          //hPulserSignal->Draw();
          //TCanvas* cPhotonCountPMT = new TCanvas("PMT"+nameTest,"PMT"+nameTest);
          //hPMTSignal->Draw();
        }

        
      }

      fileName1 = TString("ChannelOne_")+TString(name)+TString(".txt");
      fileName2 = TString("ChannelTwo_")+TString(name)+TString(".txt");
      if(nfiles%100 == 0) cout<<nfiles<<" Files processed"<<endl;
      nfiles++;
    }
    else if (nfiles < 1){
      cout<<"File Name/Path "<<dir+fileName1<<" is bad "<<endl;
      break;
    }
    else{
      cout<<"\t nfiles = "<<nfiles - 1 <<endl;
      break;
    }

    
  }
  
  TCanvas* cPhotonCount = new TCanvas("Photon Count","Photon Count");
  hPhotonCount->GetXaxis()->SetTitle("Charge");
  hPhotonCount->Draw();

  //TCanvas* cT0 = new TCanvas("T0","T0");
  //hT0->Draw();

  //TCanvas* cVoltsMax = new TCanvas("Max Volts PMT","Max Volts PMT");
  //cVoltsMax->SetLogx();
  //hVoltsMax->Draw();
  
  TCanvas* cSignal = new TCanvas("PMT Signal Integrated","PMT Signal Integrated");
  hSignal->Draw();

  TCanvas* cNoise = new TCanvas("Noise Signal Integrated","Noise Signal Integrated");
  hNoise->Draw();

  TCanvas* cAverageNoise = new TCanvas("AverageNoise","AverageNoise");
  hBaselineSlope->Draw();

  //TCanvas* cPhotonCountNoiseSubtraction = new TCanvas("hBaselineOffSet","hBaselineOffSet");
  //hBaselineOffSet->Draw();

  cout<<"Ratio of Inegral of Signal to Signal of Noise "<<hSignal->Integral("width")/hNoise->Integral("width")<<endl;
  
  outfile->Write();
}
