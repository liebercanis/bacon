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

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}



TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }


   //his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot
   his->Fit(FunName,"Q");   // fit within specified range, use ParLimits, plot, do not report
   
   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}

void tripletFit1us(Int_t dsNum = 1){
 
  Int_t irunStart = 100;//30000;//40000;//100;//10000;
  Int_t irunStop  = 105;//30230;//40176;//272;//10784;
  bool simFlag = false;
  TString outFileName = TString("1usTripletFit");
  if( irunStart == 100 && irunStop == 271) outFileName = outFileName + TString("_100ppmN2.root");
  else if( irunStart == 3000 && irunStop < 4000) outFileName = outFileName + TString("_30ShortRuns.root");
  else if( irunStart == 4000 && irunStop < 5000) outFileName = outFileName + TString("_PulseTrigger.root");
  else if (irunStart == 10000 && irunStop < 11000) outFileName = outFileName + TString("_20ppmN2.root");
  else if (irunStart == 1000 && irunStop == 1007) outFileName = outFileName + TString("_Muon.root");
  else outFileName = outFileName + TString("Temp.root");
  TFile *outfile = new TFile(outFileName ,"recreate");
  //Double_t time [] = {0,0.885,1.181,3.013,4.267,5.267,6.893/*,7.881*/,7.989,8.840,9.840,10.254,10.99,11.80,12.766,13.24,14.01,14.246,15.22,16.11,16.83,17.25,17.31,18.05,18.52,19.03,19.55,19.99,20.5,21.13,21.55,33.03,33.5};
  std::vector<Double_t> monthVector = {0,31,59,90,120,151,181,212,243,273,304,334,365};
  std::vector<Double_t> time,totalCharge,totalChargeErr,promptTime,promptTimeError,intermediateTime,intermediateTimeError,longTime,longTimeError,singleCharge,singleChargeError,A_s,A_l;
  //time.push_back(0);
  std::vector<Double_t> triplet,tripletErr;
  std::vector<Double_t> tripletUnweighted,tripletErrUnweighted;
  std::vector<Double_t> runNumber;
  Double_t globalStartTime = -999; 
  TH1D * hSinglePEAverage = new TH1D("SinglePEAverage","SinglePEAverage",200,0,2e-10);
  TH1D * hTripletAverage = new TH1D("AverageTriplet","AverageTriplet",100,1.3e-6,1.6e-6);
  TH1D * hTotalChargeAverage = new TH1D("TotalChargeAverage","TotalChargeAverage",150,0,4000);
  TH1D * hEventsPerRun = new TH1D("EventsPerRun","Events per Run",irunStop-irunStart+1,irunStart,irunStop);
  TH1D * hStartTime = new TH1D("startTimes","startTime",200,0e-9,2e-6);
  TH1D * hTmaxTimes = new TH1D("tMaxTimes","tMaxTimes ",400,0,4e-6);
    //TString filename; filename.Form("PulseAnalysis_%i.root",irunStart);
    //TString filename = "PulseAnalysisLAr_1_1.root";
    //TString filename = "Out_run_3.root";
    //TString filename = "Out_run_4.root";
    //TString filename = "Out_run_5.root";
    //TString filename = "Out_run_6.root";
  bool cutDebug = false;//true;//false;//false;//false;
  bool plotDebug = false;//false;
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
    /*
    if( z == 107) continue;
    if( z == 108) continue; 
    if( z == 109) continue; 
    */
    //TString filename = TString("dataSet1/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet2/Out_run_")+to_string(z)+TString(".root");
    //TString filename = TString("dataSet3/Out_run_")+to_string(z)+TString(".root");
    //TString filename = "SimAnaResults.124103-20181206.root;
    //TString filename = TString("Out_run_")+to_string(z)+TString(".root");
    TString filename;
    if(dsNum == 1){
      if(z>= 100 && z < 1000)
        filename = TString("processedData/DS1/SimAnaResults.baconRun_10kBins_1us_10mV_")+to_string(z)+TString(".root");
      else if(z >= 1000 && z <2000)
        filename = TString("processedData/DS1/SimAnaResults.baconRun_10kBins_1us_20mV_muon_")+to_string(z)+TString(".root");
      else if(z >= 2000 && z <4000)
        filename = TString("processedData/DS1/SimAnaResults.baconRun_10kBins_1us_20mV_div_-30mV_thresh_")+to_string(z)+TString(".root");
      else if(z >= 4000 && z <6000)
        filename = TString("processedData/DS1/SimAnaResults.baconRun_10kBins_1us_20mV_div_-6mV_thresh_")+to_string(z)+TString(".root");
      else if(z>=10000)
        filename = TString("processedData/DS1/SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_")+to_string(z)+TString(".root");
      else{ 
        //filename = "processedData/SimAnaResults.simEvents_20190308_10000.root";
        filename = "processedData/DS1/SimAnaResults.simEvent_1000.root";//processedData/SimAnaResults.simEvents_20190308_10000.root";
        simFlag = true;
        z = irunStop;
        plotDebug = true;
        //filename = TString("processedData/SimAnaResults.baconRun_10kBins_1us_20mV_div_-7.2mV_thresh_20ppmN2_HADD_")+to_string(z)+TString(".root");
      }
    }
    else if (dsNum == 2){
      if(z<100){
        filename = TString("processedData/DS2/SimAnaResults.simEvents_20190630_10000.root");
        simFlag = true;
      }
      //filename = TString("processedData/DS2/SimAnaResults");
    }
    TFile *f = new TFile(filename);
    cout<<"Run "<<z<<" opening file "<<filename<<endl;
    if(f->IsZombie()){
      cout<<"cannot open file"<<endl;
    }
   
    TNamed * eventInfo;

    f->GetObject("eventInfo",eventInfo);
    if(eventInfo == NULL){
      cout<<"Error..Warning..Fatal.."<<endl;
      cout<<"Run "<<z<<" was not properly processed"<<endl;
      cout<<"Error..Warning..Fatal.."<<endl;
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
    else if(z<30000 && z >= 10000){
      GetterStartTime = 2019.*365.+monthVector[2]+15.+12./24.+41/(24.*60);
    }
    else if(z >= 30000 && z <40000)
      GetterStartTime = 2019.*365.+monthVector[3]+6.+10./24.+31./(24.*60);
    else if(z>=40000 && 50000)
      GetterStartTime = 2019.*365.+monthVector[4]+4.+13./24.+32./(24.*60); 
    if(globalStartTime < 0) globalStartTime = globalTimeStart;
    time.push_back(globalTimeStart-GetterStartTime);
    time.push_back(globalTimeStop-GetterStartTime);
    cout<<"Time Since run zero "<<time[time.size() -2]<<", date "<<DateYear<<"-"<<DateMonth<<"-"<<DateDay<<", endTime "<<EndTimeHour<<":"<<EndTimeMinute<<", runTime "<<RunTime<<endl;

    Float_t Sdev,baseline,integral,deltaT,irun,tMax= 0,vMaxEvent,thatTotalChargeDude;
    TTree *t2 = (TTree*)f->Get("ntupleEvent");
    t2->SetBranchAddress("irun",&irun);
    t2->SetBranchAddress("Sdev",&Sdev);
    t2->SetBranchAddress("baseline",&baseline);
    t2->SetBranchAddress("integral",&integral);
    t2->SetBranchAddress("deltaT",&deltaT);
    t2->SetBranchAddress("tMax",&tMax);
    t2->SetBranchAddress("vMax",&vMaxEvent);
    t2->SetBranchAddress("totalCharge",&thatTotalChargeDude);
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
    TH1F * hIntermediateTime;
    TH1F * hLongTime; 
    if(z <10000){
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,0,300);
      hIntermediateTime = new TH1F("IntermediateTime"+runStr,"IntermediateTime"+runStr,150,0,500);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,100,700);
    }
    else{
      hTotalCharge = new TH1F("TotalCharge"+runStr,"Total_Charge_per_event_"+runStr ,150,0,450);
      hIntermediateTime = new TH1F("IntermediateTime"+runStr,"IntermediateTime"+runStr,150,0,450);
      hLongTime = new TH1F("LongTime"+runStr,"LongTime"+runStr,150,0,300);
    }
    hPromptTime = new TH1F("PromptTime"+runStr,"PromptTime"+runStr,50,0,50);
    TH1F * hArrivalTime = new TH1F("ArrivalTime weighted by charge"+runStr,"ArrivalTime weighted by charge"+runStr,400,0,10e-6);
    hArrivalTime->Sumw2();
    TH1F * hArrivalTimeUnweighted = new TH1F("ArrivalTime Unweighted"+runStr,"ArrivalTime Unweighted"+runStr,200,0,10e-6);
    TH1F * hSinglePhoton = new TH1F("SinglePhoton"+runStr,"SinglePhoton"+runStr,200,0,2e-10);
   
    Double_t singlePEThreshold = 0.99, sumSPE = 0,startTimeSPE = 7e-6;

    hArrivalTime->Reset();
    TF1 *fSinglePhoton;// = new TF1("SinglePE","landau(0)",0,2e-9);
    if(dsNum == 1)
      fSinglePhoton = new TF1("SinglePE","landau(0)",0,2e-9);
    else if(dsNum == 2)
      fSinglePhoton = new TF1("SinglePE","gaus(0)",1e-11,10e-11);
    outfile->cd();
    Int_t currentID = 0,pastID = -1,firstEntry = -9999;
    Double_t peakMax = -9999, peakTime = -9999,maxTime = -9999,Integral = 0,promptCharge = 0,intermediateCharge = 0,longCharge = 0;
    Double_t IntegralRun = 0,promptChargeRun = 0,longChargeRun = 0,eventCounter = 0;
    Double_t tMaxCut = 10e-6,tMinCut = 0.6e-6,vMaxEventCut = 2e-3,vMinCut = 2e-3,peakWidthCut = 0;//25e-9;
    Double_t shortTime = 900e-9,shortLeftWindow = 150e-9,shortRightWindow = 50e-9;
    //taken from low triplet runs at late times
    Double_t singlePE = 4.65e-11;//5.37e-11;
    Double_t singlePEerr = 1.52e-11;//1.67e-11;
    if(dsNum == 2){
      singlePE = 5.56e-11;
      singlePEerr = 1.70e-11;
    }

    currentID = pastID = eventCounter = 0;
    for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);
      t2->GetEntry(eventCounter);
      if(firstEntry < 0) firstEntry = ientry;

      currentID = ientry;
      //End of previous event Fill,
      if(pastID != currentID && firstEntry != currentID){
        if(cutDebug)cout<<"ientry "<<pastID<<", tMax "<<tMax<<", vMaxEvent "<<vMaxEvent<<", promptCharge "<<promptCharge<<", TotalCharge "<<Integral<<endl;
        
        if(singlePE && Integral > 0){
          Integral = promptCharge+intermediateCharge+longCharge;
          hTotalCharge->Fill(Integral/singlePE);
          hTotalChargeAverage->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hIntermediateTime->Fill(intermediateCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
          if(cutDebug) cout<<"ientry "<<pastID<<", Filling Charge histogram"<<endl;
        }
        hStartTime->Fill((T0));
        hTmaxTimes->Fill(tMax);
        pastID = currentID;
        eventCounter++;
        //Update the event data with the current event
        t2->GetEntry(eventCounter);
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
      
      //Fill SPE histogram before any cuts are considered
      //True SPE is taken from a sum of many runs and inserted after first analysis
      if(startTime > startTimeSPE){
        hSinglePhoton->Fill(charge);
        hSinglePEAverage->Fill(charge);
      }
      //pulse cuts go here
      //
      //cut on small peaks
      if(vMax < vMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", vMax < 0 cut, vMax = "<<vMax<<endl;
        continue;
      }
      //two sigma cut on charge
      if(charge < singlePE-2.*singlePEerr){
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
      //avoid events that start early
      if(T0 < tMinCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on tMax "<<tMax<<endl;
        continue;
      }
      //avoid events with small maximum
      if(vMaxEvent < vMaxEventCut){
        if(cutDebug)cout<<"ientry "<<ientry<<", cut on vMaxEvent "<<vMaxEvent<<endl;
        continue;
      }
      
      if(startTime > maxTime) maxTime = startTime;
      //cout<<startTime - (tMax-T0)/2.<<endl;
      //Int_t hitBin = hArrivalTime->FindBin(startTime-(tMax-T0)/2. );
      Int_t hitBin = hArrivalTime->FindBin(vMaxTime);//startTime);
      //Int_t hitBin = hArrivalTime->FindBin(startTime);
      //cout<<hitBin<<" "<<startTime - (tMax-T0)/2.<<endl;
      if(Sdev == 0)Sdev = 1.e-3;
      //err like sqrt(N)^2 + uncertainty of pulse width ^2
      //err = sqrt(charge*singlePE + std::pow(Sdev*peakWidth,2) );
      
      hArrivalTimeUnweighted->Fill(startTime);
      hArrivalTime->SetBinContent(hitBin,hArrivalTime->GetBinContent(hitBin)+charge);
      Double_t chargeErr = std::sqrt(std::pow(Sdev*peakWidth,2) + charge*singlePE);
      hArrivalTime->SetBinError(hitBin,sqrt(pow(hArrivalTime->GetBinError(hitBin),2)+chargeErr*chargeErr));
      //hArrivalTime->Fill(startTime-tMax,charge);
      hEventsPerRun->Fill(z);//,1/RunTime);
      Integral+=charge;
      IntegralRun+=charge;
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

      //Fill on very last entry
      if(i == NEntries -1){
        if(cutDebug)cout<<"ientry "<<pastID<<", tMax "<<tMax<<", vMaxEvent "<<vMaxEvent<<", promptCharge "<<promptCharge<<", TotalCharge "<<Integral<<endl;
        if(singlePE && Integral > 0){
          hTotalCharge->Fill(Integral/singlePE);
          hPromptTime->Fill(promptCharge/singlePE);
          hIntermediateTime->Fill(intermediateCharge/singlePE);
          hLongTime->Fill(longCharge/singlePE);
          if(cutDebug) cout<<"ientry "<<pastID<<", Filling Charge histogram"<<endl;
        }
        Integral = 0;
        promptCharge = 0;
        intermediateCharge = 0;
        longCharge = 0;
      }
    }
    hSinglePhoton->Fit(fSinglePhoton,"QN","",0.02e-9,.1e-9);
    singleCharge.push_back(fSinglePhoton->GetParameter(1));
    singleChargeError.push_back(fSinglePhoton->GetParameter(2));
    cout<<"SinglePE= "<<fSinglePhoton->GetParameter(1)<<"+/-"<<fSinglePhoton->GetParameter(2)<<" chi/ndf "<<fSinglePhoton->GetChisquare()/fSinglePhoton->GetNDF()<<endl;

    //Rebin is stats are to low for triplet fit
    Double_t nRunBins = hArrivalTime->GetEntries();
    Double_t tripNBins = hArrivalTime->GetNbinsX();
    if(1./sqrt(nRunBins/tripNBins) > .1){
      do{
        hArrivalTime->Rebin(2);
        tripNBins = hArrivalTime->GetNbinsX();
        if(tripNBins < 10) break;
      }
      while(1./sqrt(nRunBins/tripNBins) > .1);
    }


    Double_t  bestTrip = 0,bestErrTrip = 0,bestStart,bestSinglet = 0,bestErrSinglet = 0;
    TF1 *f1,*f2;
    
    TString fitString;
    if(plotDebug)
      fitString = TString("Q");
    else 
      fitString = TString("Q0");
    Double_t startFitTime = hArrivalTime->GetBinCenter(hArrivalTime->GetMaximumBin());
    startFitTime = 3e-6;
    
    //simulation does not have a flat component
    if(z<100)
      f1 = new TF1("myfit","([0])*exp(-x/[1])",startFitTime,maxTime);
    else
      f1 = new TF1("myfit","([0])*exp(-x/[1]) + [2]",startFitTime,maxTime);
      //f1 = new TF1("myfit","([0])*exp(-x/[1]) + ([2])*exp(-x/[3])+[4]",startFitTime,maxTime);
    
    f1->SetParLimits(0,0,1);
    f1->SetParameter(1,1e-6);
    f1->SetParLimits(1,4e-7,2e-6);
    if(z>= 100){
      f1->SetParameter(2,1e-10);
      f1->SetParLimits(2,1e-12,1e-8);
      /*
      f1->SetParLimits(2,0,1);
      f1->SetParameter(3,1e-7);
      f1->SetParLimits(3,5e-9,5e-7);
      f1->SetParameter(4,1e-10);
      f1->SetParLimits(4,1e-12,1e-8);
      */
    }
    Double_t startFitTimeUnweighted = hArrivalTimeUnweighted->GetBinCenter(hArrivalTimeUnweighted->GetMaximumBin());
    //f2 = new TF1("myfitUnweighted","([0]/[1])*exp(-x/[1]) + ([2]/[3])*exp(-x/[3])+[4]",startFitTimeUnweighted,maxTime);
    f2 = new TF1("myfitUnweighted","([0])*exp(-x/[1]) + [2]",startFitTimeUnweighted,maxTime);
    
    f2->SetParLimits(0,10,1e10);
    f2->SetParameter(1,1e-6);
    f2->SetParLimits(1,4e-7,2e-6);
    f2->SetParLimits(2,0,1e10);

    /*
    f2->SetParLimits(0,10,1e10);
    f2->SetParameter(1,1e-6);
    f2->SetParLimits(1,4e-7,2e-6);
    f2->SetParLimits(2,10,1e10);
    f2->SetParameter(3,1e-7);
    f2->SetParLimits(3,5e-9,5e-7);
    f2->SetParameter(4,1e-10);
    f2->SetParLimits(4,1e-12,1e-8);
    */
    
    hArrivalTimeUnweighted->Fit("myfitUnweighted",fitString,"",3.5e-6,maxTime-1000e-9);
    if(z>=106 && z<=130) startFitTime = 2e-6;
    hArrivalTime->Fit("myfit",fitString,"",startFitTime,maxTime-1000e-9);
    /*
    int i =1;
    do{
      cout<<"Refitting...old Triplet is "<<f1->GetParameter(1);
      if(plotDebug){
        hArrivalTime->Fit("myfit","Q","",startFitTime + 10e-9*i,maxTime-1000e-9*(i+1));
        hArrivalTimeUnweighted->Fit("myfit","Q","",startFitTimeUnweighted + 10e-9*i,maxTime-1000e-9*(i+1));
      }
      else{
        hArrivalTime->Fit("myfit",fitString,"",startFitTime + 10e-9*i,maxTime-1000e-9*(i+1));
        hArrivalTimeUnweighted->Fit("myfit",fitString,"",startFitTimeUnweighted + 10e-9*i,maxTime-1000e-9*(i+1));
      }
      cout<<"...new Triplet is "<<f1->GetParameter(1)<<", Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<", err "<<f1->GetParError(1)/f1->GetParameter(1)<<endl;
      i++;
      if(i> 10){
        cout<<"warning...fit did not get below chiSquare/ndf = 2"<<endl;
        break;
      }
    }
    while(f1->GetChisquare()/f1->GetNDF() > 3 || f1->GetParError(1)/f1->GetParameter(1) > 0.1);
    */
    bestTrip = f1->GetParameter(1);
    bestErrTrip = f1->GetParError(1);
    //bestSinglet = f1->GetParameter(3);
    //bestErrSinglet = f1->GetParError(3);
    gStyle->SetOptFit(1111); 
    triplet.push_back(bestTrip);
    tripletErr.push_back(bestErrTrip);
    //cout<<"Triplet = "<<bestTrip<<" +/- "<<bestErrTrip<<", Singlet = "<<bestSinglet<<"+/-"<<bestErrSinglet<<" Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<" ratio = "<<f1->GetParameter(0)/f1->GetParameter(2)<<endl;
    cout<<"Triplet = "<<bestTrip<<" +/- "<<bestErrTrip<<" Chi/ndf = "<<f1->GetChisquare()/f1->GetNDF()<<" ratio = "<<f1->GetParameter(0)/f1->GetParameter(2)<<endl;
    hTripletAverage->Fill(bestTrip);
    tripletUnweighted.push_back(f2->GetParameter(1));
    tripletErrUnweighted.push_back(f2->GetParError(1));


    if(z< 100) hArrivalTime->Draw();
    hSinglePhoton->Write();

    Double_t sigmaCharge = 0;
    Double_t sigmaErr = 0;
    Double_t promptErr = 0;
    Double_t intermediateTimeErr = 0;
    Double_t longChargeErr = 0;
    
    sigmaCharge  = hArrivalTime->IntegralAndError(1,hArrivalTime->GetNbinsX(),sigmaErr,"");
    sigmaCharge /= eventCounter*singlePE;
    sigmaErr/= eventCounter*singlePE;
    totalChargeErr.push_back(sigmaErr);
    
    promptCharge = hArrivalTime->IntegralAndError(shortLeftWindow-hArrivalTime->GetMaximumBin(),shortRightWindow+hArrivalTime->GetMaximumBin(),promptErr,"");
    promptCharge/= eventCounter*singlePE;
    promptErr/= eventCounter*singlePE;
    promptTimeError.push_back(promptErr);
    
    longCharge   = hArrivalTime->IntegralAndError(shortRightWindow+hArrivalTime->GetMaximumBin(),hArrivalTime->GetNbinsX(),longChargeErr,"");
    longCharge/= eventCounter*singlePE;
    longChargeErr/= eventCounter*singlePE;
    longTimeError.push_back(longChargeErr);

    totalCharge.push_back(IntegralRun/(eventCounter*singlePE));
    promptTime.push_back(promptChargeRun/(eventCounter*singlePE));
    longTime.push_back(longChargeRun/(eventCounter*singlePE));

    cout<<"Total Charge "<<IntegralRun/(eventCounter*singlePE)<<"+/-"<<sigmaErr<<", Singlet "<<promptChargeRun/(eventCounter*singlePE)<<"+/-"<<promptErr<<", Triplet "<<longChargeRun/(eventCounter*singlePE)<<"+/-"<<longChargeErr<<" NEvents"<<eventCounter<<endl;
    hTotalCharge->Write();
    hPromptTime->Write();
    hIntermediateTime->Write();
    hLongTime->Write();
    hArrivalTimeUnweighted->Write();
    hArrivalTime->Write();
  }

  std::sort(time.begin(),time.end());
  TH1D * hTripletvTime = new TH1D("Triplet_versus_Time","Triplet_versus_Time",time.size()-1,&(time[0]));
  TH1D * hTripletUnweightedvTime = new TH1D("TripletUnweighted_versus_Time","TripletUnweighted_versus_Time",time.size()-1,&(time[0]));
  TH1D * hChargevTime = new TH1D("Charge_versus_Time","Total_Photo_Electron_versus_Time",time.size()-1,&(time[0]));
  TH1D * hPromtvTime = new TH1D("Prompt_versus_Time","Promt_versus_Time",time.size()-1,&(time[0]));
  //TH1D * hIntermediatevTime = new TH1D("Intermediate_versus_Time","Intermediate_versus_Time",time.size()-1,&(time[0]));
  TH1D * hLongvTime = new TH1D("Long_versus_Time","Long_versus_Time",time.size()-1,&(time[0]));
  TH1D * hSinglePEvTime = new TH1D("SinglePE_versus_Time","SinglePE_versus_Time",time.size()-1,&(time[0]));
  Double_t tripMean = 0, tripErr = 0;
  Double_t SPEMean = 0, SPEErr = 0;
  for(int i = 0; i < time.size(); i++){
    if((i-1)%2 == 0){
      tripMean += triplet[(i-1)/2.];
      SPEMean += singleCharge[(i-1)/2.];
      SPEErr +=  singleChargeError[(i-1)/2];
      hTripletvTime->SetBinContent(i,triplet[(i-1)/2.]);
      hTripletvTime->SetBinError(i,tripletErr[(i-1)/2.]);
      hTripletUnweightedvTime->SetBinContent(i,tripletUnweighted[(i-1)/2.]);
      hTripletUnweightedvTime->SetBinError(i,tripletErrUnweighted[(i-1)/2.]);
      hChargevTime->SetBinContent(i,totalCharge[(i-1)/2]);
      hChargevTime->SetBinError(i,totalChargeErr[(i-1)/2]);
      hPromtvTime->SetBinContent(i,promptTime[(i-1)/2]);
      hPromtvTime->SetBinError(i,promptTimeError[(i-1)/2]);
      //hIntermediatevTime->SetBinContent(i,intermediateTime[(i-1)/2]);
      //hIntermediatevTime->SetBinError(i,intermediateTimeError[(i-1)/2]);
      hLongvTime->SetBinContent(i,longTime[(i-1)/2]);
      hLongvTime->SetBinError(i,longTimeError[(i-1)/2]);
      hSinglePEvTime->SetBinContent(i,singleCharge[(i-1)/2]);
      hSinglePEvTime->SetBinError(i,singleChargeError[(i-1)/2]);
    }
  }
  tripMean /= triplet.size();
  SPEMean /= singleCharge.size();
  SPEErr /= singleChargeError.size();
  std::vector<Double_t> timeMod, tripletMod, tripletModErr;
  Double_t startTime = 0,stopTime = 0;
  for(int i = 0; i < runNumber.size(); i++){
    if(runNumber[i] == 161) startTime = time[2*i];
    if(runNumber[i] == 166) stopTime = time[2*i];
    if(runNumber[i] < 161){ 
      timeMod.push_back(time[2*i]);
      timeMod.push_back(time[2*i + 1]);
      tripletMod.push_back(triplet[i]);
      tripletModErr.push_back(tripletErr[i]);
      //cout<<"Erasing "<<runNumber[i]<<" time "<<time[2*i]<<"..."<<time[2*i+1]<<" triplet "<<triplet[i]<<endl;
    }
    else if (runNumber[i] > 165){
      timeMod.push_back(time[2*i] - (stopTime-startTime));
      timeMod.push_back(time[2*i + 1] - (stopTime-startTime));
      //cout<<"time "<<time[2*i]<<", time + 1 "<<time[2*i+1]<<", delta time "<<(stopTime-startTime)<<endl;
      tripletMod.push_back(triplet[i]);
      tripletModErr.push_back(tripletErr[i]);
    }
  }
  TH1D * hTripletvTimeModified = new TH1D("Triplet_versus_TimeModified","Triplet_versus_TimeModified",timeMod.size()-1,&(timeMod[0]));
  for(int i = 0; i < timeMod.size(); i++){
    if((i-1)%2 == 0){
      //cout<<timeMod[i-1]<<"..."<<timeMod[i]<<" tripelt "<<tripletMod[(i-1)/2]<<endl;
      hTripletvTimeModified->SetBinContent(i,tripletMod[(i-1)/2.]);
      hTripletvTimeModified->SetBinError(i,tripletModErr[(i-1)/2.]);
    }
  }

  TF1 *fFilter = new TF1("Filter Rate","[0]/(1+[1]*[0]*[3]*exp(-(x)/[2]))",0,30);
  fFilter->SetParameter(0,1.6e-6);
  fFilter->SetParLimits(0,1.e-6,2e-6);
  fFilter->FixParameter(0, 1.419e-6);
  fFilter->SetParName(0,"tau'");
  //fFilter->SetParameter(1,40);
  //fFilter->SetParLimits(1,10,1e15);
  fFilter->FixParameter(1,1.11e5);
  fFilter->SetParName(1,"Birk's Constant");
  fFilter->SetParameter(2,3);
  fFilter->SetParLimits(2,0,100);
  fFilter->SetParName(2,"Filter time scale");
  if(irunStart == 100)
    fFilter->SetParameter(3,7.4);
  else if(irunStart == 10000)
    fFilter->SetParameter(3,1.4);
  else{
    fFilter->SetParameter(3,1);
    fFilter->SetParLimits(3,1,100);
  }
  fFilter->SetParName(3,"Concentration (ppm)");
  
  //Saturated Filter Rate beta = 530 ppm from WArP
  TF1 *fFilterSaturate = new TF1("Birk's Law Saturated","[0]/(1+[1]*[0]*[3]*[4]*(1-exp(-exp(-(x)/[2])/[4]) ) )",0,30);
  fFilterSaturate->SetParameter(0,1.6e-6);
  fFilterSaturate->SetParLimits(0,1.e-6,2e-6);
  fFilterSaturate->SetParName(0,"tau'");
  //fFilterSaturate->SetParameter(1,40);
  //fFilterSaturate->SetParLimits(1,10,1e15);
  fFilterSaturate->FixParameter(1,1.11e5);
  fFilterSaturate->SetParName(1,"Birk's Constant");
  fFilterSaturate->SetParameter(2,3);
  fFilterSaturate->SetParLimits(2,0,100);
  fFilterSaturate->SetParName(2,"Filter time scale");
  if(irunStart == 100)
    fFilterSaturate->SetParameter(3,7.4);
  else if(irunStart == 10000)
    fFilterSaturate->SetParameter(3,1.4);
  else{
    fFilterSaturate->SetParameter(3,1);
    fFilterSaturate->SetParLimits(3,1,100);
  }
  fFilterSaturate->FixParameter(4,530);

  Double_t tripMin = hTripletvTime->GetMinimum(1e-9);
  Double_t tripMax = hTripletvTime->GetMaximum();
  
  hTripletvTime->SetMaximum(tripMax*1.05);
  hTripletvTime->SetMinimum(tripMin*0.95);
  hTripletvTime->SetMarkerStyle(4);
  hTripletvTime->SetMarkerColor(4);
  hTripletvTime->Fit(fFilter,"Q","",0,17);
  hTripletvTime->Draw();

  TCanvas * c5 = new TCanvas("c5","c5");
  c5->cd();
  tripMin = hTripletUnweightedvTime->GetMinimum(1e-9);
  tripMax = hTripletUnweightedvTime->GetMaximum();
  hTripletAverage->Draw(); 
  hTripletUnweightedvTime->SetMaximum(tripMax*1.05);
  hTripletUnweightedvTime->SetMinimum(tripMin*0.95);
  //hTripletUnweightedvTime->Draw();
  
  if(irunStart == 100){
    //hTripletvTimeModified->Draw("same");
    hTripletvTimeModified->Fit(fFilter,"Q0","",0,hTripletvTime->GetBinContent(hTripletvTime->GetMaximumBin()));
  }
  
  TCanvas * c2 = new TCanvas("c2","c2");
  c2->cd();
  hChargevTime->SetMinimum(1);
  hChargevTime->Draw();
  hPromtvTime->Draw("same");
  hPromtvTime->SetLineColor(2);
  hLongvTime->Draw("same");
  hLongvTime->SetLineColor(3);
  //hIntermediatevTime->Draw("same");
  //hIntermediatevTime->SetLineColor(4);
  
  TH1D * hTripletvCharge = new TH1D("Tripet_versus_Photo_Electron_Yield","Triplet_v_P.E.",1000,400e-9,1600e-9);
  TH1D * hTripletvRunNumber = new TH1D("Tripet_versus_RunNumber","Triplet_v_RunNumber.",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTotalChargevRunNumber = new TH1D("Total_Charge_versus_RunNumber","Total_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hTripletChargevRunNumber = new TH1D("Triplet_Charge_versus_RunNumber","Triplet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  TH1D * hSingletChargevRunNumber = new TH1D("Singlet_Charge_versus_RunNumber","Singlet_Charge_versus_RunNumber",runNumber.size()-1,&(runNumber[0]));
  Double_t survial = 0;
  for(int i = 0; i < triplet.size(); i++){
    Int_t tripletBin = hTripletvCharge->FindBin(triplet[i]);
    hTripletvCharge->SetBinContent(tripletBin,totalCharge[i]);
    hTripletvCharge->SetBinError(tripletBin,totalChargeErr[i]);
    hTripletvRunNumber->SetBinContent(i,triplet[i]);
    hTripletvRunNumber->SetBinError(i,tripletErr[i]);
    hTotalChargevRunNumber->SetBinContent(i,totalCharge[i]);
    hTripletChargevRunNumber->SetBinContent(i,longTime[i]);
    hSingletChargevRunNumber->SetBinContent(i,promptTime[i]);
    tripErr += (tripMean-triplet[i])*(tripMean-triplet[i]);
  }
  tripErr = std::sqrt(tripErr/triplet.size());
  cout<<"mean triplet "<<tripMean<<", error "<<tripErr<<", fraction "<<tripErr/tripMean<<endl;
  cout<<"mean amplitude "<<survial/triplet.size()<<endl;
  cout<<"mean SPE "<<SPEMean<<", error "<<SPEErr<<endl;
  TCanvas * c3 = new TCanvas("c3","c3");
  c3->cd();
  hEventsPerRun->Draw();
  //hTripletvCharge->Draw();


  TCanvas * c4 = new TCanvas("c4","c4");
  c4->cd();
  hTripletvRunNumber->Draw();
  TF1 *fTripleGauss = new TF1("tripleGauss","gaus(0)+gaus(3)+gaus(6)",0,2e-10);
  fTripleGauss->SetParameters(3.70215e+03,9.93090e-12,6.57871e-12,5.06168e+03,4.60472e-11,1.52139e-11,1.98257e+03,7.03236e-11,2.76522e-11);
  fTripleGauss->SetParLimits(0,10,1e6);
  fTripleGauss->SetParLimits(1,.1e-12,2e-11);
  fTripleGauss->SetParLimits(2,1e-12,1e-10);
  fTripleGauss->SetParLimits(3,10,1e6);
  fTripleGauss->SetParLimits(4,4e-11,6e-11);
  fTripleGauss->SetParLimits(5,1e-11,1e-10);
  fTripleGauss->SetParLimits(6,10,1e6);
  fTripleGauss->SetParLimits(7,7e-11,10e-11);
  fTripleGauss->SetParLimits(8,1e-11,1e-10);

  TCanvas * c7 = new TCanvas("SPE","SPE");
  c7->cd();

  hSinglePEAverage->Fit(fTripleGauss,"","",0,1.2e-10);
  hSinglePEAverage->Draw();
  //hSinglePEAverage->Write();
  //
  TCanvas * c6 = new TCanvas("c6","c6");
  c6->cd();
  hTmaxTimes->Draw();
  hTmaxTimes->SetLineColor(2);
  hStartTime->Draw("same");
  
outfile->Write();
}
