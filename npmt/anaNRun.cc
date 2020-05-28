#include "anaNRun.hh"

anaNRun::anaNRun(Int_t dsNum,Int_t runStart,Int_t runStop)
{


  TFile * fBackground = new TFile("runAna9999_9999_DS_2.root");
  TH1D *hBackground;
  fBackground->GetObject("Charge_Run9999",hBackground);
  hBackground->Rebin(5);
  //fBackground->Close();
  //delete fBackground;


  std::vector<double> meanSPE; 
  std::vector<double> sigmaSPE; 
  if (dsNum == 1){
    meanSPE  = {4.60e-11};
    sigmaSPE = {1.5e-11};
  }
  else if(dsNum == 2 && runStart <3000 && runStart > 9999){
    meanSPE = {5.62e-11,5.38e-11};
    sigmaSPE = {1.72e-11,1.74e-11};
  }
  else if(dsNum == 2 && runStart >=3000 && runStart < 99999){
    meanSPE = {1.03800e-10,1.11016e-10};
    sigmaSPE = {3.26488e-11,3.44283e-11};
  }
  bool isSim = false;
  bool debug = false;
  //debug = true;
  
 
  //Int_t irunStop = irunStart;
  TString outFileName2 = TString("TBAcon_")+to_string(runStart)+TString("_")+to_string(runStop)+TString("_DS_")+to_string(dsNum)+TString(".root");
  TFile *outfile2 = new TFile(outFileName2,"recreate");
  outfile2->cd();
  TBacon =  new TTree("TBacon"," bacon data ");
  baconEvent = new TBaconEvent();
  TBacon->Branch("bevent",&baconEvent); 
  
  
  TString outFileName = TString("runAna")+to_string(runStart)+TString("_")+to_string(runStop)+TString("_DS_")+to_string(dsNum)+TString(".root");
  TFile *outfile = new TFile(outFileName ,"recreate");
  outfile->cd();

  TH1D * hTripletAverage[3];
  hTripletAverage[0] = new TH1D("TLong","TLong",5000,0.1e-6,4.5e-6);
  hTripletAverage[1] = new TH1D("TTransfer","TTransfer",5000,0.1e-6,4.5e-6);
  hTripletAverage[2] = new TH1D("Tslow","Tslow",5000,0.1e-6,4.5e-6);
  TH1D *hTotalChargeVrTime;
  std::vector<Double_t> chargePerRun;

  TH1D * hTripletTotalAverage  = new TH1D("TripletTotalAverage","TripletTotalAverage",400,0,10e-6);
  TH1D * hTotalChargeAverage = new TH1D("TotalChargeAverage","TotalChargeAverage",1000,10,10000);
  TH1D * hTotalChargeMuon = new TH1D("TotalChargeMuon","TotalChargeMuon",1000,10,10000);
  TH1D * hTotalChargeRand = new TH1D("TotalChargeRand","TotalChargeRand",1000,10,10000);

  TH1D *hMuVmax = new TH1D("MuVmax","muon vmax",1000,0,1);

  for(int z = runStart; z <= runStop; z++){
    if(z >= 5183 && z <= 5217) continue;
    if(z>5300 && z%2 == 0 && z<10000) continue;
    //event info from pulse finding
    TString fileName;
    TString fileDir = TString("processedData/DS")+to_string(dsNum)+TString("/");

    if(dsNum == 1){
      if(z>= 100 && z < 1000)
        fileName = TString("SimAnaResults.baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
      else if(z >= 1000 && z <2000)
        fileName = TString("SimAnaResults.baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
      else if(z >= 2000 && z <4000)
        fileName = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(z)+TString(".root");
      else if(z >= 4000 && z <6000)
        fileName = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-6mV_thresh_")+to_string(z)+TString(".root");
      else if(z>=10000)
        fileName = TString("SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+to_string(z)+TString(".root");
      else if(z == 0){
        fileName = TString("SimAnaResults.simEvents_DS_1_nPhotons_1_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 1){
        fileName = TString("SimAnaResults.simEvents_DS_1_nPhotons_10_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 2){
        fileName = TString("SimAnaResults.simEvents_DS_1_nPhotons_100_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 3){
        fileName = TString("SimAnaResults.simEvents_DS_1_nPhotons_1000_nEvents_10000.root");
        isSim = true;
      }

    }
    else if(dsNum == 2){
      if(z >= 100 && z < 199) fileName = TString("SimAnaResults.LEDRun_")+to_string(z)+TString(".root");
      else if(z == 0){
        fileName = TString("SimAnaResults.simEvents_DS_2_nPhotons_1_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 1){
        fileName = TString("SimAnaResults.simEvents_DS_2_nPhotons_10_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 2){
        fileName = TString("SimAnaResults.simEvents_DS_2_nPhotons_100_nEvents_10000.root");
        isSim = true;
      }
      else if(z == 3){
        fileName = TString("SimAnaResults.simEvents_DS_2_nPhotons_1000_nEvents_10000.root");
        isSim = true;
      }
      else if(z >= 1000 && z < 2000){
        fileName = TString("SimAnaResults.FillRun_")+to_string(z)+TString(".root");
      }
      else if(z >= 2000 && z < 9999){
        fileName = TString("SimAnaResults.RecirculateRun_")+to_string(z) + TString(".root");
      }
      else if(z == 9999){
        fileName = TString("SimAnaResults.PanelNoiseRun.root");
      }
      else if(z >= 10000 && z<19000){
        fileName = TString("SimAnaResults.MuonRun_")+to_string(z) + TString(".root");
      }
      else if(z >= 19999 && z<30000){
        fileName = TString("SimAnaResults.XenonDoping1ppm_")+to_string(z) + TString(".root");
      }
      else if(z >= 30000 && z<40000){
        fileName = TString("SimAnaResults.XenonDoping10ppm_1ppmN2_")+to_string(z) + TString(".root");
      }
      else if(z >= 40000){
        fileName = TString("SimAnaResults.XenonDoping10ppmGasTest_")+to_string(z) + TString(".root");
      }
    }
    else{
      cout<<"DSnum is not 1 or 2...you set it to "<<dsNum<<", which is not a data set number...you are so dumb. fo rel"<<endl;
    }

    //if(runStart == runStop && !isSim) debug = true;
    
    if(isSim) cout<<"This is a simulations"<<endl;
    else cout<<"This is real data.. Event number "<<z<<endl;
    TFile * infile = new TFile(fileDir+fileName);
    if(infile->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
    //---------Date and Time of run!----------//
    TNamed * eventInfo;

    infile->GetObject("eventInfo",eventInfo);
    if(eventInfo == NULL){
      cout<<"Error..Warning..Fatal.."<<endl;
      cout<<"Run "<<z<<" was not properly processed"<<endl;
      cout<<"Error..Warning..Fatal.."<<endl;
      if(!isSim) continue;
      else eventInfo = new TNamed("title","2019-02-22 13:57 0.03999996185");
    }
    TH1D * hSummedWaveform[2];
    infile->GetObject("SumWaveform0",hSummedWaveform[0]);
    hSummedWaveform[0]->SetTitle(Form("SumWaveform0_%i",z));
    infile->GetObject("SumWaveform1",hSummedWaveform[1]);
    hSummedWaveform[1]->SetTitle(Form("SumWaveform1_%i",z));
    TF1 * fTripFitSummed = new TF1("LongExpoPlusConstantSummed","([0])*exp(-x/[1]) + [2]",2.5e-6,9.5e-6);
    Double_t maxVal = hSummedWaveform[0]->GetMaximum();
    Double_t minVal = hSummedWaveform[0]->GetMinimum();
    fTripFitSummed->SetParameter(0,maxVal);
    fTripFitSummed->SetParLimits(0,0.1*maxVal,2*maxVal);;
    fTripFitSummed->SetParameter(1,1e-6);
    fTripFitSummed->SetParLimits(1,0.5e-6,2.5e-6);
    fTripFitSummed->SetParameter(2,minVal);
    if(minVal < 0) minVal = 0;
    fTripFitSummed->SetParLimits(2,0,10);
    for(int i = 0; i < 2;i++){
      if(i == 0){
        hSummedWaveform[i]->Fit("LongExpoPlusConstantSummed","Q0","",3e-6,9e-6);
        tripletLifetimePMTSum0.push_back(fTripFitSummed->GetParameter(1));
        tripletLifetimeErrorPMTSum0.push_back(fTripFitSummed->GetParError(1));
        cout<<"Summed Waveform Triplet fit PMT 0 "<<fTripFitSummed->GetParameter(1)<<" ";
      }
      else if(i == 1 && 0){
        hSummedWaveform[i]->Fit("LongExpoPlusConstantSummed","Q0","",3e-6,9e-6);
        tripletLifetimePMTSum1.push_back(fTripFitSummed->GetParameter(1));
        tripletLifetimeErrorPMTSum1.push_back(fTripFitSummed->GetParError(1));
        cout<<"...PMT 1 "<<fTripFitSummed->GetParameter(1)<<endl;
      }
      //hSummedWaveform[i]->Draw();
    }
    string fullTitle = eventInfo->GetTitle();
    double DateYear = stod(fullTitle.substr(0,4));
    double DateMonth = stod(fullTitle.substr(5,2));
    double DateDay = stod(fullTitle.substr(8,2));
    double EndTimeHour = stod(fullTitle.substr(11,2));
    double EndTimeMinute = stod(fullTitle.substr(14,2));
    double RunTime = stod(fullTitle.substr(25) );//,fullTitle.size() - 1);
    double globalTimeStop = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.) + RunTime/(3600.*24.);
    double globalTimeStart = DateYear*365.+ monthVector[DateMonth] + DateDay + EndTimeHour/24. +EndTimeMinute/(24.*60.);
    double GetterStartTime = 0;
    //Year,month,day,hour,minute -- all in  units of days
    //GetterStartTime = 2019.*365.+monthVector[11]+4.+14./24.+10/(24.*60);
    if(z >=5000 && z < 10000) 
      GetterStartTime = 2019.*365.+monthVector[11]+19.+11./24.+21/(24.*60);
    else if(z >= 19999 && z< 30000 ) 
      GetterStartTime = 2019.*365.+monthVector[12]+9.+11./24.+16/(24.*60);
    else if(z >=30000 ){
      GetterStartTime = 2020.*365.+monthVector[1]+4.+9./24.+30/(24.*60);
    }
    //cout<<"GetterStartTime "<<GetterStartTime<<", globalTimeStart "<<globalTimeStart<<", globalTimeStop "<<globalTimeStop;
    if(z == 19999){
      runTime.push_back(globalTimeStart-GetterStartTime);
      runTime.push_back(globalTimeStart-GetterStartTime+0.11999500/24);
      cout<<"Length of run "<<(0.11999500)<<"(hours)"<<", time since start/stop "<<globalTimeStart-GetterStartTime
        <<"/"<<globalTimeStart-GetterStartTime+0.11999500/24 <<endl;
    }else if(z == 20000){
      runTime.push_back(globalTimeStart-GetterStartTime+0.11999500/24);
      runTime.push_back(globalTimeStop-GetterStartTime );
      cout<<"Length of run "<<(globalTimeStop-globalTimeStart)*24.<<"(hours)"<<", time since start/stop "<<globalTimeStart-GetterStartTime+0.11999500/24
        <<"/"<<globalTimeStop-GetterStartTime<<endl;
    }
    else{
      runTime.push_back(globalTimeStart-GetterStartTime);
      runTime.push_back(globalTimeStop-GetterStartTime);
      cout<<"Length of run "<<(globalTimeStop-globalTimeStart)*24.<<"(hours)"<<", time since start/stop "<<globalTimeStart-GetterStartTime
        <<"/"<<globalTimeStop-GetterStartTime<<endl;
    }



    if(debug) cout<<"Loading ProcessedTree"<<endl;
    TTree * processedTree = new TTree();
    infile->GetObject("processedTree",processedTree);

    if(processedTree == NULL){
      cout<<"ProcessedTree is NULL meaning the file needs to be rerun through anaSim or the run number does not exist!"<<endl;
      continue;
    }
    if(debug) cout<<"Loading pmtRun"<<endl;
    pmtRun = new TPmtRun();
    processedTree->SetBranchAddress("pmtRun",&pmtRun);
    if(debug) cout<<"Loading processedEvents"<<endl;
    pmtEvent = new TPmtEvent();
    processedTree->SetBranchAddress("processedEvent",&pmtEvent);

    if(debug) cout<<"Creating Histograms"<<endl;
    TH1D* hTripletChargeWeightedSummed = new TH1D(Form("TripletChargeWeightedSummed_Run%i",z),Form("TripletChargeWeightedSummed_Run%i",z),400,0,10e-6);
    TH1D* hTripletChargeWeightedPMT[2];
    hTripletChargeWeightedPMT[0] = new TH1D(Form("TripletChargeWeightedPMT0_Run%i",z),Form("TripletChargeWeightedPMT0_Run%i",z),400,0,10e-6);
    hTripletChargeWeightedPMT[1] = new TH1D(Form("TripletChargeWeightedPMT1_Run%i",z),Form("TripletChargeWeightedPMT1_Run%i",z),4000,0,4e-6);
    TH1D* hCharge           = new TH1D(Form("Charge_Run%i",z),Form("Charge_Run%i",z),1000,10,10000);
    TH1D* hF40Summed                   = new TH1D(Form("F40Summed_Run%i",z),Form("F40Summed_Run%i",z),100,0,1);
    
    TH1D* hSPE[2];
    hSPE[0] = new TH1D(Form("SPE_PMT0_Run%i",z),Form("SPESummed_Run%i",z),100,0,3e-9);
    hSPE[1] = new TH1D(Form("SPE_PMT1_Run%i",z),Form("SPESummed_Run%i",z),100,0,3e-9);

    Int_t NEvents = (Int_t)processedTree->GetEntries();
    Int_t eventCounter = NEvents;
    if(debug) cout<<"Processing Events"<<endl;
    outfile->cd();

    
    
    Double_t tMaxCut = 5e-6,tMinCut = 0.00e-6,vMaxEventCut = 10,vMinCut = 1e-3,peakWidthCut = 0,nHits = 10;//25e-9;
    for(int ientry = 0; ientry < NEvents; ientry++){
      if(z == 19999 && ientry > 500) continue;
      if(z == 20000 && ientry < 500) continue;
      //Process pulse finding information
      processedTree->GetEntry(ientry);
      Int_t nPMTs = pmtRun->nPMTs;
      Double_t deltaT = pmtRun->deltaT;

      Int_t triggerTime = 0;
      //Process waveform information
      std::vector<Double_t> volts;
      std::vector<Double_t> time = pmtEvent->time;

      for(int j = 0; j < nPMTs;j++){
        baconEvent->clear();
        if(debug) cout<<"Looping over PMTs"<<endl;
        if(j == 0) volts = pmtEvent->volt1;
        else if(j == 1) volts = pmtEvent->volt2;

        Int_t nPulses = pmtRun->charge[j].size();
        Double_t T0 = pmtRun->T0[j];
        Double_t C0 = 0;
        if(nPulses > 0) C0 = pmtRun->charge[j][0];
        Double_t totalCharge = pmtRun->totalCharge[j];
        Double_t tMax = pmtRun->tMax[j];
        Double_t vMax = pmtRun->vMax[j];
        Double_t muonVmax = 9999;
        Double_t startCharge = 0;
        for(int i = 0; i < nPulses;i++){
          Double_t charge     = pmtRun->charge[j][i];
          Double_t startTime  = pmtRun->startTimes[j][i];
          Double_t peakWidth  = pmtRun->peakWidths[j][i];
          if(startTime > 9e-7) break;
          startCharge += charge;
        }
        if(startCharge/totalCharge > 0.05)continue;
        if(pmtRun->vMax.size() == 2){
          muonVmax = pmtRun->vMax[1];
        }
        hMuVmax->Fill(muonVmax);
        if(muonVmax > 500e-3 && muonVmax != 9999 && j == 0){
          hTotalChargeMuon->Fill(totalCharge/meanSPE[j]);
        }
        else if(muonVmax != 9999 && j == 0){
          hTotalChargeRand->Fill(totalCharge/meanSPE[j]);
        }
        triggerTime += T0;
        if(j == 0) hCharge->Fill(totalCharge/meanSPE[j]);

        //if(z < 19999 || z > 20010)
        //if(z < 20020 || z > 20026)
        //if(z < 20040 || z > 20050)
        if( (z < 20060 || z > 20070) && j == 0){
          hTotalChargeAverage->Fill(totalCharge/meanSPE[j]);
        }
        if(nPulses < nHits && j == 0){
          if(debug) cout<<"nHits cut...nPulses "<<nPulses<<", nHits cut "<<nHits<<", PMT "<<j<<endl;
          eventCounter--;
          continue;
        }
        //avoid events that start late
        if(tMax > tMaxCut && j == 0){
          if(debug) cout<<"tMaxCut "<<tMaxCut<<", pulse value "<<tMax<<", PMT "<<j<<endl;
          eventCounter--;
          continue;
        }

        baconEvent->event=ientry;
        baconEvent->run=z;
        baconEvent->npulse=nPulses;
        baconEvent->npmt=j;
        baconEvent->totQ = pmtRun->totalCharge[j];
        baconEvent->spe = meanSPE[j];
        baconEvent->muVmax=muonVmax;


        baconEvent->T0=pmtRun->T0[j];
        baconEvent->totalCharge=pmtRun->totalCharge[j];
        baconEvent->tMax=       pmtRun->tMax[j];
        baconEvent->vMax=       pmtRun->vMax[j];
        baconEvent->cMax=       pmtRun->cMax[j];
        baconEvent->baseline=   pmtRun->baseline[j];
        baconEvent->sDev=       pmtRun->sDev[j]; 


        if(ientry%1000 == 0 || ientry == NEvents - 1) printf(" ev %i  pmt %i totQ %E \n ",ientry, j,baconEvent->totQ);

        for(int i = 0; i < nPulses;i++){
          if(debug && i == 0) cout<<"Looping over Pulses"<<endl;
          Double_t charge     = pmtRun->charge[j][i];
          Double_t startTime  = pmtRun->startTimes[j][i];
          Double_t peakWidth  = pmtRun->peakWidths[j][i];
          Double_t peakHeight = pmtRun->peakHeights[j][i];
          Double_t peakTime   = pmtRun->peakTimes[j][i];
          TPulse thePulse;
          thePulse.time=startTime;
          thePulse.tpeak=peakTime;
          thePulse.pwidth=peakWidth;
          thePulse.peak=peakHeight;
          thePulse.q=charge;
          //thePulse.qerr=phitQErr;


          //SPE Fill
          if(peakTime > 7e-6) hSPE[j]->Fill(charge);

          //cut on small peaks
          debug = false;
          if(peakHeight < vMinCut && j == 0){
            if(debug) cout<<"vMinCut "<<vMinCut<<", pulse value "<<peakHeight<<", PMT "<<j<<endl;
            continue;
          }
          if(peakWidth < peakWidthCut && j == 0) {
            if(debug) cout<<"peakWidthCut "<<peakWidthCut<<", pulse value "<<peakWidth<<", PMT "<<j<<endl;
            continue;
          }

          //avoid events that start early
          if(startTime < tMinCut && j == 0){
            if(debug) cout<<"tMinCut "<<tMinCut<<", pulse value "<<startTime<<", PMT "<<j<<endl;
            continue;
          }
          //avoid events with large maximum
          //i.e. events that go out of range
          if(vMax >= vMaxEventCut && j == 0){
            if(debug) cout<<"vMaxEventCut "<<vMaxEventCut<<", pulse value "<<vMax<<", PMT "<<j<<endl;
            continue;
          }
          //*/
          //Fill Triplet
          //if(peakTime > 6e-6 && peakTime < 7e-6) cout<<"charge "<<charge<<endl;
          Int_t peakTimeBin = hTripletChargeWeightedSummed->FindBin(peakTime);
          Int_t nPhotons = charge/meanSPE[j];
          //triplet fit for both channels
          hTripletChargeWeightedSummed->Fill(peakTime,charge/meanSPE[j]);
          //hTripletChargeWeightedSummed->Fill(peakTime,charge);
          //hTripletChargeWeightedPMT[j]->Fill(peakTime,charge);
          if(j == 0){
            hTripletChargeWeightedPMT[0]->Fill(peakTime,charge/meanSPE[j]);
            hTripletChargeWeightedPMT[1]->Fill(peakTime,charge/meanSPE[j]);
          }
          baconEvent->hits.push_back(thePulse);
        }
        //if(debug) cout<<"F40 cut "<<endl;
        //pulse finding loop
        //F40 cut time window
        Double_t F40 = 0;
        Double_t F40Window = 40e-9;
        Double_t F40Start = 1e-6,F40Stop = F40Start+ F40Window;
        for(int i = F40Start/deltaT; i <= F40Stop/deltaT; i++){
          F40 += -volts[i]*deltaT;  
        }
        //if(debug) cout<<" Filling final histogram"<<endl;
        hF40Summed->Fill(F40/totalCharge);
        baconEvent->F40=F40;
        TBacon->Fill();
      }

      if(ientry%1000 == 0 || ientry == NEvents - 1) cout<<"\t processed "<<ientry<<" events "<< TBacon->GetEntries() 
        << " totQ 0 " << pmtRun->totalCharge[0]
          << " totQ 1 " << pmtRun->totalCharge[1] << endl;
      //if((ientry+1)%(int(NEvents*0.33)) == 0 || ientry == NEvents - 1 ) cout<<"\tprocessed "<<ientry+1<<" events"<<endl;
    }

    Double_t startFit = 2.5e-6;
    Double_t stopFit  = 10e-6;
    TF1 *fSingleExpoPlusConstant = new TF1("SingleExpoPlusConstant","([0])*exp(-x/[1]) + [2]",1.3e-6,stopFit);

    Double_t tripMaxVal = hTripletChargeWeightedPMT[0]->GetMaximum();
    Double_t tripMinVal = hTripletChargeWeightedPMT[0]->GetMinimum();
    fSingleExpoPlusConstant->SetParameter(0,tripMaxVal);
    fSingleExpoPlusConstant->SetParLimits(0,0.1*tripMaxVal,2*tripMaxVal);;
    fSingleExpoPlusConstant->SetParameter(1,.1e-6);
    fSingleExpoPlusConstant->SetParLimits(1,0.1e-6,2.5e-6);
    fSingleExpoPlusConstant->SetParameter(2,tripMinVal);
    fSingleExpoPlusConstant->SetParLimits(2,0.1*tripMinVal,10*tripMinVal);
   
    TF1 *fSingleExpo = new TF1("SingleExpo","([0])*exp(-x/[1])",startFit,stopFit);
    
    fSingleExpo->SetParameter(0,tripMaxVal);
    fSingleExpo->SetParLimits(0,0.1*tripMaxVal,2*tripMaxVal);
    fSingleExpo->SetParameter(1,1e-6);
    fSingleExpo->SetParLimits(1,4e-7,2e-6);

    TF1 *fDoubleExpoPlusConstant = new TF1("DoubleExpoPlus","([0])*exp(-x/[1]) +[2] - ([3])*exp(-x/[4])+([5])*exp(-x/[6])",1e-6,stopFit);
    fDoubleExpoPlusConstant->SetParameter(0,tripMaxVal);
    fDoubleExpoPlusConstant->SetParLimits(0,1e-2*tripMaxVal,tripMaxVal);
    fDoubleExpoPlusConstant->SetParameter(1,0.1e-6);
    fDoubleExpoPlusConstant->SetParLimits(1,1e-6,4e-6);
    fDoubleExpoPlusConstant->SetParameter(2,tripMinVal);
    fDoubleExpoPlusConstant->SetParLimits(2,1e-2*tripMinVal,1e2*tripMinVal); 
    fDoubleExpoPlusConstant->SetParameter(3,tripMaxVal);
    fDoubleExpoPlusConstant->SetParLimits(3,.01*tripMaxVal,1e3*tripMaxVal);
    fDoubleExpoPlusConstant->SetParameter(4,3.5e-7);
    fDoubleExpoPlusConstant->SetParLimits(4,1e-7,1e-6);
    fDoubleExpoPlusConstant->SetParameter(5,tripMaxVal);
    fDoubleExpoPlusConstant->SetParLimits(5,.1*tripMaxVal,1e4*tripMaxVal);
    fDoubleExpoPlusConstant->SetParameter(6,3.5e-7);
    fDoubleExpoPlusConstant->SetParLimits(6,1e-7,1e-6);

    TF1 *fDoubleExpoPlusConstantMod = new TF1("DoubleExpoPlusMod","([0])*exp(-x/[1]) +[2] - ([3])*exp(-x/[4])",1e-6,stopFit);
    fDoubleExpoPlusConstantMod->SetParameter(0,tripMaxVal);
    fDoubleExpoPlusConstantMod->SetParLimits(0,1e-2*tripMaxVal,tripMaxVal);
    fDoubleExpoPlusConstantMod->SetParameter(1,0.1e-6);
    fDoubleExpoPlusConstantMod->SetParLimits(1,1e-6,4e-6);
    fDoubleExpoPlusConstantMod->SetParameter(2,tripMinVal);
    fDoubleExpoPlusConstantMod->SetParLimits(2,1e-2*tripMinVal,1e2*tripMinVal); 
    fDoubleExpoPlusConstantMod->SetParameter(3,tripMaxVal);
    fDoubleExpoPlusConstantMod->SetParLimits(3,tripMaxVal,1e3*tripMaxVal);
    fDoubleExpoPlusConstantMod->SetParameter(4,3.5e-7);
    fDoubleExpoPlusConstantMod->SetParLimits(4,1e-7,1e-6);
    //debug = true;
    //TF1 *fAfterPulsing = new TF1("AfterPulsing","([0])*exp(-x/[1]) + [2] + [3]*TMath::Landau(x,[4],[5])",2e-6,10e-6);
    //TF1 *fAfterPulsing = new TF1("AfterPulsing","([0])*exp(-x/[1]) + [2] + [3]*TMath::Gaus(x,[4],[5])",1.5e-6,10e-6);
    TF1 *fAfterPulsing = new TF1("AfterPulsing","([0])*exp(-x/[1]) + [2] + [3]*TMath::Gaus(x,[4],[5]) + [6]*TMath::Landau(x,[7],[8])",1.25e-6,10e-6);
    fAfterPulsing->SetParameter(0,tripMaxVal);
    fAfterPulsing->SetParLimits(0,0.1*tripMaxVal,2*tripMaxVal);
    fAfterPulsing->SetParameter(1,1e-6);
    fAfterPulsing->SetParLimits(1,4e-7,2.5e-6);
    fAfterPulsing->SetParameter(2,tripMinVal);
    fAfterPulsing->SetParLimits(2,0.1*tripMinVal,10*tripMinVal);
    //Gaus const, mean, sigma
    fAfterPulsing->SetParameter(3,0.01*tripMaxVal);
    fAfterPulsing->SetParLimits(3,0.001*tripMaxVal,0.1*tripMaxVal);
    fAfterPulsing->SetParameter(4,2.e-6);
    fAfterPulsing->SetParLimits(4,1.8e-6,2.2e-6);
    fAfterPulsing->SetParameter(5,4e-8);
    fAfterPulsing->SetParLimits(5,1e-8,10e-8);
    //Landau const, MPV, sigma
    fAfterPulsing->SetParameter(6,0.1*tripMaxVal);
    fAfterPulsing->SetParLimits(6,0.001*tripMaxVal,10.*tripMaxVal);
    fAfterPulsing->SetParameter(7,1.5e-6);
    fAfterPulsing->SetParLimits(7,1.25e-6,1.75e-6);
    fAfterPulsing->SetParameter(8,5e-8);
    fAfterPulsing->SetParLimits(8,1e-9,10e-8);    
      
    TCanvas * cTripletFit = NULL;
    if(debug) cTripletFit = new TCanvas("Triplet Fit","Triplet Fit");
    else delete cTripletFit;

    Double_t fitTriplet0 = 0, fitTripletError0 = 0;
    Double_t fitTriplet1 = 0, fitTripletError1 = 0, const1 = 0, offset = 0;
    Double_t tripletMean = 0;//fitTriplet1 = fDoubleExpoPlusConstant->GetParameter(6);
    Double_t tripletSigma = 0;//fitTripletError1 = fDoubleExpoPlusConstant->GetParError(6);


    TString fitStr = "";
    if(debug) fitStr ="";
    else fitStr = "Q";
    if(isSim) 
      hTripletChargeWeightedSummed->Fit("SingleExpoPlusConstant","","",startFit,stopFit);
    else{
      hTripletChargeWeightedSummed->Fit("SingleExpoPlusConstant",fitStr,"",startFit,stopFit);
      //hTripletChargeWeightedPMT[0]->Fit("SingleExpoPlusConstant",fitStr,"",startFit,stopFit);
      //fitTriplet0 = fSingleExpoPlusConstant->GetParameter(1);
      //fitTripletError0 = fSingleExpoPlusConstant->GetParError(1);
      hTripletChargeWeightedPMT[0]->Fit("AfterPulsing",fitStr,"",1.5e-6,stopFit);
      //refit if something goes wrong
      if(5 < fAfterPulsing->GetChisquare()/fAfterPulsing->GetNDF()){
        cout<<"refitting because ChiSquare/NDF is "<<fAfterPulsing->GetChisquare()/fAfterPulsing->GetNDF()
          <<" Moving fit to start after first after pulsing "<<endl;

        hTripletChargeWeightedPMT[0]->Fit("AfterPulsing",fitStr,"",1.75e-6,stopFit);      
      }
      //hTripletChargeWeightedPMT[0]->Fit("SingleExpoPlusConstant",fitStr,"",startFit,stopFit);
      //fitTriplet1 = fSingleExpoPlusConstant->GetParameter(1);
      //fitTripletError1 = fSingleExpoPlusConstant->GetParError(1);
      if(z>20019){
        hTripletChargeWeightedPMT[0]->Fit("DoubleExpoPlus",fitStr,"",1.05e-6,stopFit);
        fitTriplet1 = fDoubleExpoPlusConstant->GetParameter(4);
        fitTripletError1 = fDoubleExpoPlusConstant->GetParError(4);
        fitTriplet0 = fDoubleExpoPlusConstant->GetParameter(1);//fAfterPulsing->GetParameter(1);
        fitTripletError0 = fDoubleExpoPlusConstant->GetParError(1);//fAfterPulsing->GetParError(1);
        tripletMean = fDoubleExpoPlusConstant->GetParameter(6);
        tripletSigma = fDoubleExpoPlusConstant->GetParError(6);
      }
      else if(z>=20000){
        hTripletChargeWeightedPMT[0]->Fit("DoubleExpoPlusMod",fitStr,"",1.2e-6,stopFit);
        fitTriplet1 = fDoubleExpoPlusConstantMod->GetParameter(4);
        fitTripletError1 = fDoubleExpoPlusConstantMod->GetParError(4);
        tripletMean = fDoubleExpoPlusConstantMod->GetParameter(1);
        tripletSigma = fDoubleExpoPlusConstantMod->GetParError(1);

      }
      /*
      else if(z >= 5000 && z < 9999){
        fitTriplet0 = 
        fitTriplet1 = 
        tripletMean = 
      }
      */
      //const1 = fSingleExpoPlusConstant->GetParameter(0);
      //offset = fSingleExpoPlusConstant->GetParameter(2);

    }

    //double x1 = 1.3e-6,x2 = 10e-6;
    //double val = offset*x2 - (const1*fitTriplet1)*exp(-x2/fitTriplet1) - (offset*x1 - (const1*fitTriplet1)*exp(-x1/fitTriplet1));
    //cout<<"Trip "<<fitTriplet1<<", const "<<const1<<", offset "<<offset<<" integral "<<val<<endl;

    Int_t bin1 = hTripletChargeWeightedPMT[0]->FindBin(900e-9);
    Int_t bin2 = hTripletChargeWeightedPMT[0]->FindBin(1300e-9);
    Int_t binFinal = hTripletChargeWeightedPMT[0]->FindBin(10e-6);
    cout<<"Short time integral "<<hTripletChargeWeightedPMT[0]->Integral(bin1,bin2,"width");
    cout<<", Long  time integral "<<fSingleExpoPlusConstant->Integral(1.3e-6,10e-6)<<endl;
    
    Double_t total = hTripletChargeWeightedPMT[0]->Integral(bin1,bin2,"width") + fSingleExpoPlusConstant->Integral(1.3e-6,10e-6);
    cout<<"Total charge "<<binFinal*total/((10e-6-900e-9)*eventCounter)<<
      ", short "<<binFinal*hTripletChargeWeightedPMT[0]->Integral(bin1,bin2,"width")/((10e-6-900e-9)*eventCounter) <<
      ", long "<<binFinal*fSingleExpoPlusConstant->Integral(1.3e-6,10e-6)/((10e-6-900e-9)*eventCounter) <<endl;
    
    //if(z < 19999 || z > 20010)
    //if(z < 20020 || z > 20026)
    //if(z < 20040 || z > 20050)
    if(z < 20060 || z > 20070){

      hTripletTotalAverage->Add(hTripletChargeWeightedPMT[0]);
      cout<<"adding a histogem with Z = "<<z<<endl;
    }

    hTripletAverage[0]->Fill(fitTriplet0);
    tripletLifetimePMT0.push_back(fitTriplet0);
    tripletLifetimeErrorPMT0.push_back(fitTripletError0);
    hTripletAverage[1]->Fill(fitTriplet1);
    tripletLifetimePMT1.push_back(fitTriplet1);
    tripletLifetimeErrorPMT1.push_back(fitTripletError1);
    //Double_t tripletMean = (fitTriplet0+fitTriplet1)/2.;
    //Double_t tripletSigma = sqrt(pow(tripletMean-fitTriplet0,2)+pow(tripletMean-fitTriplet1,2));

    //Double_t tripletMean = fitTriplet1;
    //Double_t tripletSigma = fitTripletError1;
    tripletLifetime.push_back(tripletMean);
    hTripletAverage[2]->Fill(tripletMean);
    tripletLifetimeError.push_back(tripletSigma);
    
    cout<<"T_long "<<fitTriplet0<<"+/-"<<fitTripletError0<<", T_D "<<fitTriplet1<<"+/-"<<fitTripletError1<<", T_slow "<<tripletMean<<"+/-"<<tripletSigma<<endl;

    hCharge->Rebin(5);
    hBackground->Scale(hCharge->GetMaximum()/hBackground->GetMaximum());
    hCharge->Add(hBackground,-1);
    hCharge->GetXaxis()->SetRangeUser(10,3000);
    //outfile->cd();
    TF1 *fLandau = new TF1("myLandau","[0]*TMath::Landau(x,[1],[2])",10,3000);
    fLandau->SetParLimits(0,1,1e5);
    fLandau->SetParameter(0,1000);
    fLandau->SetParLimits(1,200,700);
    fLandau->SetParameter(1,400);
    fLandau->SetParLimits(2,10,500);
    fLandau->SetParameter(2,200);
    hCharge->Fit("myLandau","","",150,3000);
    chargePerRun.push_back(fLandau->GetParameter(1));
    //cout<<"Triplet Fit values for AfterPulsing function :"<<fitTriplet0<<", Single Exponetial w/ constant "<<
      //fitTriplet1<<", averaged "<<tripletMean<<"+/-"<<tripletSigma<<endl;
       //Write Histogram for each run
    hTripletChargeWeightedSummed->Write();
    hCharge->Write();
    hTripletChargeWeightedPMT[0]->Write();
    hTripletChargeWeightedPMT[1]->Write();
    hSummedWaveform[0]->Write();
    hSummedWaveform[1]->Write();

    hSPE[0]->Write();
    if(dsNum == 2)
      hSPE[1]->Write();
    hF40Summed->Write();

    if(debug){
      cTripletFit->cd();
      //What is a logy?
      cTripletFit->SetLogy();
      //hTripletChargeWeightedSummed->Draw();
      hTripletChargeWeightedPMT[0]->Draw("");
      //hTripletChargeWeightedPMT[0]->SetLineColor(2);
      //hTripletChargeWeightedPMT[1]->Draw("same");

      TCanvas * cTotalCharge = new TCanvas("TotalCharge","TotalCharge");
      cTotalCharge->cd();
      hCharge->Draw();

      TCanvas * cSPE = new TCanvas("SPE","SPE");
      cSPE->cd();
      hSPE[0]->Draw();
      hSPE[1]->SetLineColor(2);
      hSPE[1]->Draw("same");

      TCanvas * cF40 = new TCanvas("F40","F40");
      cF40->cd();
      hF40Summed->Draw();
    }


  //Final Loop over all runs
  }
  hTripletvTime = new TH1D("Triplet_versus_Time","Triplet_versus_Time",runTime.size()-1,&(runTime[0]));
  
  hTripletvTimeSum0 = new TH1D("Triplet_versus_Time_PMTSUM0","Triplet_versus_Time_PMTSUM0",runTime.size()-1,&(runTime[0]));
  hTripletvTimeSum1 = new TH1D("Triplet_versus_Time_PMTSUM1","Triplet_versus_Time_PMTSUM1",runTime.size()-1,&(runTime[0])); 

  hTotalChargeVrTime = new TH1D("TotalCharge_Versus_Time","TotalCharge_Versus_Time",runTime.size()-1,&(runTime[0]));
  for(unsigned int i = 0; i < runTime.size(); i++){
    if((i-1)%2 == 0){
      hTripletvTime->SetBinContent(i,tripletLifetimePMT0[(i-1)/2.]);
      hTripletvTime->SetBinError(i,tripletLifetimeErrorPMT0[(i-1)/2.]);

      hTripletvTimeSum0->SetBinContent(i,tripletLifetimePMT1[(i-1)/2.]);
      hTripletvTimeSum0->SetBinError(i,tripletLifetimeErrorPMT1[(i-1)/2.]);

      hTripletvTimeSum1->SetBinContent(i,tripletLifetime[(i-1)/2.]);
      hTripletvTimeSum1->SetBinError(i,tripletLifetimeError[(i-1)/2.]);

      hTotalChargeVrTime->SetBinContent(i,chargePerRun[(i-1)/2.]);
      hTotalChargeVrTime->SetBinError(i,sqrt(chargePerRun[(i-1)/2.]));
    }
  }

  Double_t tripMin = hTripletvTime->GetMinimum(1e-9);
  Double_t tripMax = hTripletvTime->GetMaximum();
  
  hTripletvTime->SetMaximum(tripMax*1.05);
  hTripletvTime->SetMinimum(tripMin*0.95);
  hTotalChargeVrTime->GetYaxis()->SetRangeUser(200,1.1*hTotalChargeVrTime->GetMaximum());
  hTotalChargeVrTime->Draw();

  TF1 *fFilter = new TF1("Filter Rate","[0]/(1+[1]*[0]*[3]*exp(-(x)/[2]))",0,30);
  fFilter->SetParameter(0,1.6e-6);
  fFilter->SetParLimits(0,1.e-6,2e-6);
  fFilter->FixParameter(0, 1.487e-6);
  fFilter->SetParName(0,"tau'");
  //fFilter->SetParameter(1,40);
  //fFilter->SetParLimits(1,10,1e15);
  fFilter->FixParameter(1,1.225e5);
  fFilter->SetParName(1,"Birk's Constant");
  
  fFilter->SetParameter(2,3);
  fFilter->SetParLimits(2,0,100);
  fFilter->SetParName(2,"Filter time scale");
  
  fFilter->SetParLimits(3,0.1,10);
  fFilter->SetParameter(3,1);
  fFilter->SetParName(3,"Concentration (ppm)");

  hTripletvTime->Write();
  hTripletvTimeSum0->Write();
  hTripletvTimeSum1->Write();
  Double_t scaleInt = hTripletTotalAverage->Integral(0,400);
  //hTripletTotalAverage->Scale(1/scaleInt);
  hTripletTotalAverage->Write();
  hTotalChargeMuon->Write();
  hTotalChargeRand->Write();
  hTotalChargeVrTime->Write();

  //Write Histograms for Data set
  hTotalChargeAverage->Write();
  //hTripletAverage->Write();
  outfile->Write();
  outfile2->Write();
}
