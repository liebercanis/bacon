#include "anaSPE.hh"

anaSPE::anaSPE(Int_t startRun,Int_t stopRun,Int_t dsNum = 2){

  TString outFileName = TString("/home/nmcfadde/RooT/PMT/processedData/DS")+to_string(dsNum)+TString("/SPE")+to_string(startRun)+TString(".root");
  TFile *outFile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
 
  Double_t timeBin = .8e-9,nBins = 1000;
  TH1D* hSignal[2];
  hSignal[0] = new TH1D(TString("SPE_0"),"",nBins,0,nBins*timeBin);
  hSignal[1] = new TH1D(TString("SPE_1"),"",nBins,0,nBins*timeBin);

  std::vector<Double_t> speCounter;speCounter.resize(2);
  for(int z = startRun; z<=stopRun;z++)
  {
    //time 12:36:26 133626
    ///date 24/12/1997 19971224
    Double_t globalStartTime = -999;
    TString fileDir = TString("/home/nmcfadde/RooT/PMT/rootData/DS")+to_string(dsNum)+TString("/");
    TString fileName;
    if(dsNum == 1){
      if(z >= 1000 && z <2000)
        fileName = TString("baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
      else if(z>=2000 && z < 4000)
        fileName = TString("baconRun_10kBins_1us_20mV_div_-30mV_thresh_") + to_string(z) +TString(".root");
      else if(z>=4000 && z < 6000)
        fileName = TString("baconRun_10kBins_1us_20mV_div_-6mV_thresh_") + to_string(z) +TString(".root");
      else if(z>=10000)
        fileName = TString("baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+ to_string(z) +TString(".root");
      else if (z >= 100 && z <=272)
        fileName = TString("baconRun_10kBins_1us_10mV_")+ to_string(z) +TString(".root");
      else
        fileName = "simEvents_20190508_10000.root";
    }
    else if(dsNum == 2){
      if(z >= 100 && z < 199) fileName = TString("LEDRun_")+to_string(z)+TString(".root");
    }


    TFile *fin = new TFile(fileDir+fileName,"readonly");

    if(fin->IsZombie()) {
      printf(" File is zombie %s%s\n",fileDir.Data(),fileName.Data());
      continue;
    }
    else 
      printf(" looking for file %s%s\n",fileDir.Data(),fileName.Data());

    TNamed * eventInfo;

    fin->GetObject("eventInfo",eventInfo);
    if(eventInfo != NULL)
      cout<<"eventInfo is "<<eventInfo->GetTitle()<<endl;
    else if(z != 0){
      cout<<"Skipping file "<<fileName<<" because it has an attitude problem"<<endl;
      continue;
    }

    // get pmtTree from file 
    TTree * pmtTree = new TTree();
    fin->GetObject("pmtTree",pmtTree);
    Long64_t nentries = pmtTree->GetEntries();
    cout << " number of entries is " << nentries << endl;

    // set up memory for reading
    simEvent = new TPmtSimulation();
    pmtEvent = new TPmtEvent();
    pmtTree->SetBranchAddress("pmtSimulation", &simEvent);
    cout<<"Simulation branch set"<<endl;
    pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
    cout<<"Event branch set"<<endl;

    TDatime time;

    //switch to output file

    //event info from pulse finding
    TString filename;
    TString filedir = TString("/home/nmcfadde/RooT/PMT/processedData/DS")+to_string(dsNum)+TString("/");
    if(dsNum == 1){
      if(z>= 100 && z < 1000)
        filename = TString("SimAnaResults.baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
      else if(z >= 1000 && z <2000)
        filename = TString("SimAnaResults.baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
      else if(z >= 2000 && z <4000)
        filename = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(z)+TString(".root");
      else if(z >= 4000 && z <6000)
        filename = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-6mV_thresh_")+to_string(z)+TString(".root");
      else if(z>=10000)
        filename = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+to_string(z)+TString(".root");
      else{ 
        filename = "SimAnaResults.simEvent_1000.root";//processedData/SimAnaResults.simEvents_20190308_10000.root";
      }
    }
    else if(dsNum == 2){
      if(z >= 100 && z < 199) filename = TString("SimAnaResults.LEDRun_")+to_string(z)+TString(".root");
    }
    else{
      cout<<"DSnum is not 1 or 2...you set it to "<<dsNum<<", which is not a data set number...you are so dumb. fo rel"<<endl;
    }

    TFile *f = new TFile(filedir+filename);
    cout<<"Run "<<z<<" opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
   
    f->GetObject("eventInfo",eventInfo);
    if(eventInfo == NULL){
      cout<<"Error..Warning..Fatal.."<<endl;
      cout<<"Run "<<z<<" was not properly processed"<<endl;
      cout<<"Error..Warning..Fatal.."<<endl;
    }
    string fullTitle = eventInfo->GetTitle();
     // set up memory for reading
    TTree * processedTree = new TTree();
    f->GetObject("processedTree",processedTree);
    processedEvent = new TPmtEvent();
    processedTree->SetBranchAddress("processedEvent", &processedEvent);
    cout<<"ProcessedEvent branch set"<<endl;
       ///*
    double DateYear = stod(fullTitle.substr(0,4));
    double DateMonth = stod(fullTitle.substr(5,2));
    double DateDay = stod(fullTitle.substr(8,2));
    double EndTimeHour = stod(fullTitle.substr(11,2));
    double EndTimeMinute = stod(fullTitle.substr(14,2));
    double RunTime = stod(fullTitle.substr(25) );//,fullTitle.size() - 1);
    double globalTimeStop = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.) + RunTime/(3600.*24.);
    double globalTimeStart = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.);
    double GetterStartTime = 0;
    if(z< 1000){
      //019-1-11, endTime 11:27//run 126
      GetterStartTime = 2019.*365.+monthVector[1]+11.+11./24.+10/(24.*60);
    }
    else if(z<30000 && z >= 10000){
      GetterStartTime = 2019.*365.+monthVector[2]+15.+12./24.+41/(24.*60);
    }
    else if(z >= 30000 && z <40000)
      GetterStartTime = 2019.*365.+monthVector[3]+6.+10./24.+31./(24.*60);
    else if(z>=40000 && 50000)
      GetterStartTime = 2019.*365.+monthVector[4]+4.+13./24.+32./(24.*60); 
    if(globalStartTime < 0) globalStartTime = globalTimeStart;

    Float_t Sdev,baseline,integral,deltaT,irun,tMax= 0,vMaxEvent;
    TTree *t2 = (TTree*)f->Get("ntupleEvent");
    t2->SetBranchAddress("irun",&irun);
    t2->SetBranchAddress("Sdev",&Sdev);
    t2->SetBranchAddress("baseline",&baseline);
    t2->SetBranchAddress("integral",&integral);
    t2->SetBranchAddress("deltaT",&deltaT);
    t2->SetBranchAddress("tMax",&tMax);
    t2->SetBranchAddress("vMax",&vMaxEvent);
    //intialize branch
    Int_t NEvents = (Int_t)t2->GetEntries();
    t2->GetEntry(0);

    TTree *t1 = (TTree*)f->Get("ntuplePulse");

    Float_t ientry,pmtNum,sigma,nhits,charge,startTime,peakWidth,vMax,vMaxTime,T0,sigmaP,sigmaQ; 

    t1->SetBranchAddress("ientry",&ientry);
    t1->SetBranchAddress("pmtNum",&pmtNum);
    t1->SetBranchAddress("nhits",&nhits);
    t1->SetBranchAddress("charge",&charge);
    t1->SetBranchAddress("startTime",&startTime);
    t1->SetBranchAddress("vMax",&vMax);
    t1->SetBranchAddress("vMaxTime",&vMaxTime);
    t1->SetBranchAddress("T0",&T0);
    //t1->SetBranchAddress("Sdev",&Sdev);
    t1->SetBranchAddress("peakWidth",&peakWidth);
    Int_t NEntries = (Int_t)t1->GetEntries();
    outFile->cd();
    for(int i = 0; i < NEntries; i++){
      //if(ientry == 1277)  break;
      t1->GetEntry(i);
      t2->GetEntry(2*ientry-pmtNum);
      pmtTree->GetEntry(ientry);
      processedTree->GetEntry(ientry);
      //Get waveform
      std::vector<std::vector<Double_t> > signal;signal.resize(2);
      
      /*
      if(timeBin != deltaT){
        cout<<"deltaT "<<deltaT<<" does not match timeBin "<<timeBin<<endl;
        return;
      }
      */
      /*
      signal[0] = pmtEvent->volt1;
      signal[1] = pmtEvent->volt2;
      std::vector<Double_t> time  = pmtEvent->time;
      Int_t nBins = pmtEvent->volt1.size();
      //Double_t deltaT = (pmtEvent->time[NEvents -1] - pmtEvent->time[0])/(NEvents-1);
      Double_t deltaT = (pmtEvent->time[nBins -1] - pmtEvent->time[0])/(nBins-1);
      //cout<<"PMT info...nBins "<<nBins<<", deltaT "<<deltaT<<endl;
      */
      signal[0] = processedEvent->volt1;
      signal[1] = processedEvent->volt2;
      std::vector<Double_t> time  = processedEvent->time;
      Int_t nBins = processedEvent->volt1.size();
      //Double_t deltaT = (processedEvent->time[NEvents -1] - processedEvent->time[0])/(NEvents-1);
      Double_t deltaT = (processedEvent->time[nBins -1] - processedEvent->time[0])/(nBins-1);
      //cout<<"PMT info...nBins "<<nBins<<", deltaT "<<deltaT<<endl;


      Int_t startBin = startTime/deltaT ;
      Int_t stopBin = (startTime+peakWidth)/deltaT ;
      Int_t peakBin = vMaxTime/deltaT-1;
      Double_t fudgeFactor = vMax/signal[pmtNum][peakBin];
      Double_t q = 0;
      //raw data does not WMA taken from baseline
      for(int j = startBin; j <= stopBin; j++) q += (signal[pmtNum][j]);
      fudgeFactor = 1./(q*deltaT/charge);
      fudgeFactor = 1.;
      //Volts  Channel Mean    Sigma
      //1.3kV  0       7.4e-12 2.3e-12
      //1.3kV  1       7.5e-12 2.2e-12
      //1.4kV  0       1.1e-11 3.6e-12
      //1.4kV  1       1.2e-11 3.4e-12
      //1.5kV  0       1.9e-11 6.2e-12 
      //1.5kV  1       2.0e-11 7.6e-12
      //1.6kV  0       3.2e-11 1.13e-11
      //1.6kV  1       3.6e-11 1.4e-11
      //******Gain of ~10 *****
      //1.2kV  0       3.2e-11 1.2e-11
      //1.2kV  1       3.5e-11 1.4e-11
      //1.3kV  0       6.0e-11 2.0e-11
      //1.3kV  1       6.7e-11 2.2e-11
      //1.4kV  0       1.2e-10 4.2e-11
      //1.4kV  1       1.3e-10 4.5e-11
      //1.5kV  0       2.1e-10 7.5e-11
      //1.5kV  1       2.2e-10 7.3e-11
      //1.6kV  0       3.5e-10 1.1e-10
      //1.6kV  1       3.5e-10 1.2e-10
      //10k resistor swapped for 250k
      //1.6kV  0       1.5e-11 5.48e-12
      //1.6kV  1       1.6e-11 6.0e-12
      //Grounded to Can at PMT + aluminum shield (fix to X-talk)
      //1.6kV  0       3.5e-11 1.3e-11
      //1.6kV  1       4.3e-11 1.4e-11
      //Swapped 50 -> 200, 200 -> 500
      //1.6kV  0       3.4e-11 1.2e-11
      //1.6kV  1       4.0e-11 1.5e-11
      //Double_t mean = 1.9e-11,sigma = 6.2e-12;
      std::vector<Double_t> mean,sigma;
      if(dsNum == 1){
        mean  = {4.60e-11};
        sigma = {1.5e-11};
      }
      if(dsNum == 2){
        if(startRun == 113){
          mean = {1.9e-11,2.0e-11};
          sigma = {6.2e-12,7.6e-12};
        }
        else if(startRun == 114){
          mean = {3.2e-11,3.6e-11};
          sigma = {1.13e-11,1.4e-11};
        }
        else if(startRun == 115){
          mean = {1.1e-11,1.2e-11};
          sigma = {3.6e-12,3.4e-12};
        }
        else if(startRun == 116){
          mean = {7.4e-12,7.5e-12};
          sigma = {2.3e-12,2.2e-12};
        }
        else if(startRun ==119 || startRun == 121){
          mean = {6.6e-11,7.2e-11};
          sigma = {2.5e-11,2.5e-11};
        }
        else if(startRun == 120){
          mean = {3.2e-11,3.5e-11};
          sigma = {1.2e-11,1.4e-11};
        }
        else if(startRun == 122){
          mean = {1.2e-10,1.3e-10};
          sigma = {4.2e-11,4.5e-11};
        }
        else if(startRun == 123){
          mean = {2.1e-10,2.2e-10};
          sigma = {7.5e-11,7.3e-11};
        }
        else if(startRun == 124){
          mean = {3.5e-10,3.5e-10};
          sigma = {1.2e-10,1.1e-11};
        }
        else if(startRun == 125){
          mean = {1.5e-11,1.6e-11};
          sigma = {5.5e-12,6e-12};
        }
        else if(startRun == 126){
          mean = {3.5e-11,4.3e-11};
          sigma = {1.3e-11,1.4e-11};
        }
        else if(startRun == 127){
          mean = {3.4e-11,4.0e-11};
          sigma = {1.3e-11,1.5e-11};
        }
        else if(startRun == 128){
          mean = {3.57e-11,4.0e-11};
          sigma = {1.22e-11,1.2e-11};
        }
        else if(startRun == 129){
          mean = {5.56e-11,6.23e-11};
          sigma = {1.50e-11,1.90e-11};
        }
        else if(startRun == 132){
          mean = {5.62e-11,5.38e-11};
          sigma = {1.72e-11,1.74e-11};
        }


      }
      Double_t timeCut = 0;//7e-6 for data
      if(dsNum == 1) timeCut = 7e-6;
      for(int j = startBin; j <= stopBin+20; j++){
        if(j == signal[pmtNum].size()) continue;
        if(startTime > timeCut && charge > mean[pmtNum] - sigma[pmtNum] && charge < mean[pmtNum] + sigma[pmtNum]){
          Int_t shiftBin = peakBin - 30;
          Double_t binVal = hSignal[(int)pmtNum]->GetBinContent(j+1-shiftBin);
          hSignal[(int)pmtNum]->SetBinContent(j+1-shiftBin,-signal[pmtNum][j]+binVal);
          //cout<<pmtNum<<" "<<-signal[pmtNum][j]<<" "<<hSignal[(int)pmtNum]->GetBinContent(j+1-startBin+30)<<" "<<j+1-startBin+30<<endl;
        }
      }
      if(startTime > timeCut && charge > mean[pmtNum] - sigma[pmtNum] && charge < mean[pmtNum] +sigma[pmtNum]) speCounter[pmtNum]++;
      //cout<<"start Time "<<startTime<<" endTime "<<startTime+peakWidth<<", peakValue "<<-vMax<<endl;
      //cout<<"Found start "<<processedEvent->time[startBin]<<", foundStop "<<processedEvent->time[stopBin]<<", peak value "<<signal[pmtNum][peakBin]
      //<<" peakBin "<<peakBin<<", vMaxTime "<<vMaxTime<<endl;
      if(i%1000 == 0) cout<<ientry<<endl;
    }
  }
  for(int j =0; j < 2; j++){
    for(int i = 0; i < hSignal[j]->GetNbinsX(); i++){
      Double_t binVal = hSignal[j]->GetBinContent(i+1);
      hSignal[j]->SetBinContent(i+1,binVal/speCounter[j]);
      //cout<<"Channel "<<j<<" signal "<<binVal<<" speCounter "<<speCounter[j]<<endl;
    }
  }
  cout<<speCounter[0]<<" single photons were found channel 0"<<endl;
  cout<<speCounter[1]<<" single photons were found channel 1"<<endl;
  hSignal[0]->Draw();
  if(hSignal[1]->GetMaximum() > hSignal[0]->GetMaximum()){
    hSignal[0]->SetMaximum(hSignal[1]->GetMaximum()*1.10);
    //cout<<"Channel 1 max "<<hSignal[1]->GetMaximum()<<endl;
  }
  hSignal[1]->SetLineColor(2);
  hSignal[1]->Draw("same");
  outFile->Write();
}
