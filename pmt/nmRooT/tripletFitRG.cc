//#include "anaEvents.hh"

//anaEvents::ananEvents( Int_t irun){
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
#include <string>

void tripletFitRG(){
 
  //enter file run number and code will loop over data sets
  Int_t irunStart = 30100;
  Int_t irunStop  = 30120;
  bool simFlag = false;

  TString outFileName = TString("1usTripletFit");
  if( irunStart == 100 && irunStop == 272) outFileName = outFileName + TString("_100ppmN2.root");
  else if( irunStart == 3000 && irunStop < 4000) outFileName = outFileName + TString("_30ShortRuns.root");
  else if( irunStart == 4000 && irunStop < 5000) outFileName = outFileName + TString("_PulseTrigger.root");
  else if (irunStart == 10000 && irunStop < 11000) outFileName = outFileName + TString("_20ppmN2.root");
  else outFileName = outFileName + TString("Temp.root");
  TFile *outfile = new TFile(outFileName ,"recreate");
  std::vector<Double_t> monthVector = {0,31,59,90,120,151,181,212,243,273,304,334,365};
  std::vector<Double_t> time,totalCharge,totalChargeErr,promptTime,promptTimeError,intermediateTime,intermediateTimeError,longTime,longTimeError;
  std::vector<Double_t> triplet,tripletErr;
  std::vector<Double_t> runNumber;
  Double_t globalStartTime = -999; 
  TH1D * hSinglePEAverage = new TH1D("SinglePEAverage","SinglePEAverage",200,0,2e-10); 
  TH1D * hTotalChargeAverage = new TH1D("TotalChargeAverage","TotalChargeAverage",150,0,450);
  
  bool cutDebug = false;

  for(int z =irunStart;z<=irunStop;z++){
    if(z>=20000 && z<30000) continue;
    if( z == 10117) continue;
    if( z == 10274) continue;
    if( z == 10330) continue;
    if( z == 10353) continue;
    if( z == 10358) continue;
    if( z == 5002) continue;
    if( z == 2046) continue;
    if( z == 185) continue;
    if( z == 159) continue;
    //these runs are for cobalt
    if( z >= 262 && z <= 271) continue;
    if( z == 107) continue;
    if( z == 108) continue;
    if( z == 109) continue;
    //complex way of reading run data
    TString filename;
    if(z>= 100 && z < 1000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
    else if(z >= 1000 && z <2000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
    else if(z >= 2000 && z <4000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(z)+TString(".root");
    else if(z >= 4000 && z <6000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-6mV_thresh_")+to_string(z)+TString(".root");
    else if(z>=10000)
      filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+to_string(z)+TString(".root");
    else{ 
      //filename = "processedData/SimAnaResults.simEvents_20190308_10000.root";
      filename = "processedData/SimAnaResults.simEvents_20190308_10000.root";
      simFlag = true;
      z = irunStop;
    }
    TFile *f = new TFile(filename);
    cout<<"Run "<<z<<" opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
   
    TNamed * eventInfo;

    f->GetObject("eventInfo",eventInfo);
    if(eventInfo == NULL){
      cout<<"Error......"<<endl;
      cout<<"Run "<<z<<" was not properly processed or does not exist"<<endl;
      if(!simFlag) continue;
      else eventInfo = new TNamed("title","2019-02-22 13:57 0.03999996185");

    }
    runNumber.push_back(z);
    string fullTitle = eventInfo->GetTitle();
    //cout<<"eventInfo is "<<eventInfo->GetTitle()<<endl;
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
    else if(z >= 10000){
      GetterStartTime = 2019.*365.+monthVector[2]+15.+12./24.+41/(24.*60);
    }
    if(globalStartTime < 0) globalStartTime = globalTimeStart;
    time.push_back(globalTimeStart-GetterStartTime);
    time.push_back(globalTimeStop-GetterStartTime);
    cout<<"Time Since run zero "<<time[time.size() -2]<<", date "<<DateYear<<"-"<<DateMonth<<"-"<<DateDay<<", endTime "<<EndTimeHour<<":"<<EndTimeMinute<<", runTime "<<RunTime<<endl;

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
    //t1->SetBranchAddress("baseline",&baseline);

    Int_t NEntries = (Int_t)t1->GetEntries();
    TString runStr = to_string(z);
    TH1F * hTotalCharge;
    TH1F * hPromptTime;
    TH1F * hLongTime; 
    if(z <10000){
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,100,700);
      hPromptTime = new TH1F("PromptTime"+runStr,"PromptTime"+runStr,100,20,200);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,100,700);
    }
    else{
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,0,450);
      hPromptTime = new TH1F("PromptTime"+runStr,"PromptTime"+runStr,100,0,150);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,0,300);
    }
    TH1F * hArrivalTime = new TH1F("ArrivalTime weighted by charge"+runStr,"ArrivalTime weighted by charge"+runStr,400,0,10e-6);
    hArrivalTime->Sumw2();
    TH1F * hArrivalTimeUnweighted = new TH1F("ArrivalTime Unweighted"+runStr,"ArrivalTime Unweighted"+runStr,200,0,10e-6);
    TH1F * hSinglePhoton = new TH1F("SinglePhoton"+runStr,"SinglePhoton"+runStr,200,0,2e-10);
   
    hArrivalTime->Reset();
    TF1 *fSinglePhoton = new TF1("SinglePE","landau(0)",0,2e-9);
    outfile->cd();
    Int_t currentID = 0,pastID = -1;
    std::vector<Double_t> eventSkip;eventSkip.resize(NEvents);
    std::vector<Double_t> maxPeakTime;maxPeakTime.resize(NEvents);
    Double_t peakMax = -9999, peakTime = -9999,maxTime = -9999,Integral = 0,promptCharge = 0,intermediateCharge = 0,longCharge = 0;
    Double_t IntegralRun = 0,promptChargeRun = 0,longChargeRun = 0;
    Double_t tMaxCut = 1.25e-6,tMinCut = 0.6e-6,vMaxEventCut = 6e-3,vMinCut = 2e-3,peakWidthCut = 0;//25e-9;
    //taken from low triplet runs at late times
    Double_t singlePE = 5.37e-11;//fSinglePhoton->GetParameter(1);
    Double_t singlePEerr = 1.67e-11;//fSinglePhoton->GetParameter(2);
    //TF1 *hGaus = new TF1("Gaus","gaus(0)",0,200);

    //Double_t meanNhits = hGaus->GetParameter(1);
    //Double_t sigmaNhits = hGaus->GetParameter(2);
    std::vector<Double_t> deltaQ,Q,dQNorm;
    deltaQ.resize(hArrivalTime->GetNbinsX());Q.resize(hArrivalTime->GetNbinsX());dQNorm.resize(hArrivalTime->GetNbinsX());
    //cout<<hGaus->GetChisquare()/hGaus->GetNDF()<<endl;
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      t2->GetEntry(ientry);
      ///*
      //pulse cuts go here
      if(vMax < vMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", vMax < 0 cut, vMax = "<<vMax<<endl;
        continue;
      }
      //one sigma cut on charge
      if(charge < singlePE-2*singlePEerr){
        if(cutDebug)cout<<"ientry "<<ientry<<", one sigma cut on charge... "<<charge<<", less than "<<singlePE-singlePEerr<<endl;
        continue;
      }
      if(peakWidth < peakWidthCut) {
        if(cutDebug)cout<<"ientry "<<ientry<<", Negative Peak Width"<<endl;
        continue;
      }
      //avoid events that start late
      if(tMax > tMaxCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on tMax "<<tMax<<endl;
        continue;
      }
      //avoid events that start late
      if(T0 < tMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on tMax "<<tMax<<endl;
        continue;
      }
      if(vMaxEvent < vMaxEventCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on vMaxEvent "<<vMaxEvent<<endl;
        continue;
      }
      //*/

      currentID = ientry;
      //End of event Fill,
      if(pastID != currentID && currentID != 0){
        pastID = currentID;
        if(cutDebug)cout<<"ientry "<<ientry<<", tMax "<<tMax<<", vMaxEvent "<<vMaxEvent<<endl;
        if(singlePE && promptCharge > 0 && longCharge > 0){
          Integral = promptCharge+intermediateCharge+longCharge;
          hTotalCharge->Fill(Integral/singlePE);
          hTotalChargeAverage->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
        }
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
      //*/
      if(startTime > maxTime) maxTime = startTime;
      Int_t hitBin = hArrivalTime->FindBin(startTime-tMax);
      if(Sdev == 0)Sdev = 1.e-3;
      hArrivalTimeUnweighted->Fill(startTime,charge);
      hArrivalTime->Fill(startTime-tMax,charge);
      Integral+=charge;
      Double_t shortTime = 900e-9,shortLeftWindow = 150e-9,shortRightWindow = 150e-9;
      if(startTime > shortTime - shortLeftWindow && startTime < shortTime + shortRightWindow){
        promptCharge += charge;
        promptChargeRun += charge;
      }
      else if( startTime > 900e-9 +100e-9 && startTime < 900e-9 +200e-9){
        intermediateCharge += charge;
      }
      else if( startTime > shortTime + shortRightWindow){
        longCharge   += charge;
        longChargeRun += charge;
      }
      if(i == NEntries -1){
        if(singlePE){
          hTotalCharge->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
        }
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
      //*/
    }
    
    Double_t bestChi = 9999999,bestTrip = 0,bestErrTrip = 0,bestStart,bestSinglet = 0,bestErrSinglet = 0;
    TF1 *f1;

    Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin());
    if(z>1000)
      f1 = new TF1("myfit","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])",startFitTime,maxTime);
    else
      f1 = new TF1("myfit","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])",startFitTime,maxTime);
    f1->SetParLimits(0,0,1e-10);
    f1->SetParameter(1,1e-6);
    f1->SetParLimits(1,4e-7,2e-6);
    f1->SetParLimits(2,0,1e-10);
    f1->SetParameter(3,1e-7);
    f1->SetParLimits(3,5e-9,5e-7);
    f1->SetParLimits(4,0,1);

    //hArrivalTime->Fit("myfit","Q","",0e-9,maxTime-1000e-9);
    cout<<"startFitTime "<<startFitTime<<", stopTime "<<maxTime-1000e-9<<endl;
    hArrivalTime->Fit("myfit","Q","",startFitTime,maxTime-1000e-9);
    int i =1;
    do{
      cout<<"Refitting...old Triplet is "<<f1->GetParameter(1);
      hArrivalTime->Fit("myfit","Q","",startFitTime + 10e-9*i,maxTime-1000e-9*(i+1));
      cout<<"...new Triplet is "<<f1->GetParameter(1)<<", Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<", err "<<f1->GetParError(1)/f1->GetParameter(1)<<endl;
      i++;
      if(i> 10){
        cout<<"warning...fit did not get below chiSquare/ndf = 2"<<endl;
        break;
      }
    }
    while(f1->GetChisquare()/f1->GetNDF() > 2 || f1->GetParError(1)/f1->GetParameter(1) > 0.1);
    if(z >= 110 && z <= 140) hArrivalTime->Fit("myfit","Q","",50e-9,2000e-9);
    bestTrip = f1->GetParameter(1);
    bestErrTrip = f1->GetParError(1);
    bestSinglet = f1->GetParameter(3);
    bestErrSinglet = f1->GetParError(3);
    gStyle->SetOptFit(1111); 
    triplet.push_back(bestTrip);
    //tripletErr.push_back(bestErrTrip+bestTrip*0.0272821);
    tripletErr.push_back(bestErrTrip);
    cout<<"Triplet = "<<bestTrip<<" +/- "<<bestErrTrip<<", Singlet = "<<bestSinglet<<"+/-"<<bestErrSinglet<<" Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<" ratio = "<<f1->GetParameter(0)/f1->GetParameter(2)<<endl;

    hSinglePhoton->Write();
    Double_t sigmaCharge = 0;
    Double_t sigmaErr = 0;
    Double_t promptErr = 0;
    Double_t intermediateTimeErr = 0;
    Double_t longChargeErr = 0;
    sigmaCharge  = hArrivalTime->IntegralAndError(1,hArrivalTime->GetNbinsX(),sigmaErr,"");
    sigmaCharge /= NEvents*singlePE;
    //totalCharge.push_back(sigmaCharge);
    sigmaErr/= NEvents*singlePE;
    totalChargeErr.push_back(sigmaErr);
    promptCharge = hArrivalTime->IntegralAndError(1,hArrivalTime->FindBin(50e-9),promptErr,"");
    promptCharge/= NEvents*singlePE;
    //promptTime.push_back(promptCharge);
    promptErr/= NEvents*singlePE;
    promptTimeError.push_back(promptErr);
    longCharge   = hArrivalTime->IntegralAndError(hArrivalTime->FindBin(50e-9),hArrivalTime->GetNbinsX(),longChargeErr,"");
    longCharge/= NEvents*singlePE;
    //longTime.push_back(longCharge);
    longChargeErr/= NEvents*singlePE;
    longTimeError.push_back(longChargeErr);
    
    totalCharge.push_back(IntegralRun/(NEvents*singlePE));
    promptTime.push_back(promptChargeRun/(NEvents*singlePE));
    longTime.push_back(longChargeRun/(NEvents*singlePE));



    cout<<"Total Charge "<<IntegralRun/(NEvents*singlePE)<<"+/-"<<sigmaErr<<", Singlet "<<promptChargeRun/(NEvents*singlePE)<<"+/-"<<promptErr<<", Triplet "<<longChargeRun/(NEvents*singlePE)<<"+/-"<<longChargeErr<<endl;
    
    hTotalCharge->Write();
    hPromptTime->Write();
    hLongTime->Write();
    hArrivalTimeUnweighted->Write();
    hArrivalTime->Write();
  }
  //end of loop over run data
  std::sort(time.begin(),time.end());
  TH1D * hTripletvTime = new TH1D("Triplet_versus_Time","Triplet_versus_Time",time.size()-1,&(time[0]));
  TH1D * hChargevTime = new TH1D("Charge_versus_Time","Total_Photo_Electron_versus_Time",time.size()-1,&(time[0]));
  TH1D * hPromtvTime = new TH1D("Prompt_versus_Time","Promt_versus_Time",time.size()-1,&(time[0]));
  TH1D * hLongvTime = new TH1D("Long_versus_Time","Long_versus_Time",time.size()-1,&(time[0]));
  Double_t tripMean = 0, tripErr = 0;
  cout<<"Filling time histograms"<<endl;
  for(int i = 0; i < time.size(); i++){
    if((i-1)%2 == 0){
      tripMean += triplet[(i-1)/2.];
      hTripletvTime->SetBinContent(i,triplet[(i-1)/2.]);
      hTripletvTime->SetBinError(i,tripletErr[(i-1)/2.]);
      hChargevTime->SetBinContent(i,totalCharge[(i-1)/2]);
      hChargevTime->SetBinError(i,totalChargeErr[(i-1)/2]);
      hPromtvTime->SetBinContent(i,promptTime[(i-1)/2]);
      hPromtvTime->SetBinError(i,promptTimeError[(i-1)/2]);
      //hIntermediatevTime->SetBinContent(i,intermediateTime[(i-1)/2]);
      //hIntermediatevTime->SetBinError(i,intermediateTimeError[(i-1)/2]);
      hLongvTime->SetBinContent(i,longTime[(i-1)/2]);
      hLongvTime->SetBinError(i,longTimeError[(i-1)/2]);
    }
  }
  tripMean /= triplet.size();
  std::vector<Double_t> timeMod, tripletMod, tripletModErr;
  Double_t startTime = 0,stopTime = 0;
  
  hTripletvTime->SetMarkerStyle(4);
  hTripletvTime->SetMarkerColor(4);
  hTripletvTime->Draw();
  
  TCanvas * c2 = new TCanvas("c2","c2");
  c2->cd();
  hChargevTime->Draw();
  hPromtvTime->Draw("same");
  hPromtvTime->SetLineColor(2);
  hLongvTime->Draw("same");
  hLongvTime->SetLineColor(3);
  
  TCanvas * cA = new TCanvas("cA","cA");
  cA->cd();
  if(irunStart == 100)
  cout<<"Filling Run Number histograms"<<endl;
  TH1D * hTripletvCharge = new TH1D("Tripet_versus_Photo_Electron_Yield","Triplet_v_P.E.",1000,400e-9,1600e-9);
  TH1D * hTripletvRunNumber = new TH1D("Tripet_versus_RunNumber","Triplet_v_RunNumber.",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTotalChargevRunNumber = new TH1D("Total_Charge_versus_RunNumber","Total_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTripletChargevRunNumber = new TH1D("Triplet_Charge_versus_RunNumber","Triplet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hSingletChargevRunNumber = new TH1D("Singlet_Charge_versus_RunNumber","Singlet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));

  for(int i = 0; i < triplet.size(); i++){
    hTripletvRunNumber->SetBinContent(i,triplet[i]);
    hTripletvRunNumber->SetBinError(i,tripletErr[i]);
    hTotalChargevRunNumber->SetBinContent(i,totalCharge[i]);
    hTripletChargevRunNumber->SetBinContent(i,longTime[i]);
    hSingletChargevRunNumber->SetBinContent(i,promptTime[i]);
    tripErr += (tripMean-triplet[i])*(tripMean-triplet[i]);
  }
  tripErr = std::sqrt(tripErr/triplet.size());
  cout<<"mean triplet "<<tripMean<<", error "<<tripErr<<", fraction "<<tripErr/tripMean<<endl;
  TCanvas * c3 = new TCanvas("c3","c3");
  c3->cd();
  hTotalChargevRunNumber->Draw();
  hTripletChargevRunNumber->Draw("same");
  hSingletChargevRunNumber->Draw("same");


  TCanvas * c4 = new TCanvas("c4","c4");
  c4->cd();
  hTripletvRunNumber->Draw();
  //hSinglePEAverage->Write();
  
outfile->Write();
}
