/* MG revised */
#include "anaCompass.hh"
#include <TFitResult.h>

struct hit_t {
  int nev;
  int nhits;
  int order;
  int istart;
  int nwidth;
  int good;
  int kind;
  float time;
  float thit;
  float q;
  float qerr;
  float peak;
  float qsum;
  float q900;
};
struct hit_t theHit;

anaCompass::anaCompass(TString tag, Int_t maxEvents)
{
  isSimulation=false;
  bool debug=false;
  printf(" starting anaCompass tag %s \n",tag.Data());
  
  TString fileName; fileName.Form("rootData/compass/compass_run_%s.root",tag.Data());
 // TString fileName; fileName.Form("rootData/DS4/run_31/UNFILTERED/compass_run_%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fChain = new TChain("Data");
  fChain->Add(fileName.Data());
  fChain->GetListOfFiles()->ls();
  fCurrent = -1;
  Init();
  fChain->GetListOfBranches()->ls();
  Long64_t nentries = fChain->GetEntriesFast();

  int printInterval=100;
  int nHists=0;
  unsigned noPeakEventCount=0;
  nMaxHistEvents=10;
  Int_t nMaxPulseHist=100;
  Int_t pulseHistCount=0;
  firstChargeCut=0.1;
  lifeChargeCut=0.025;
  microSec=1.0E6;
  ran = new TRandom3();
  Double_t simHitMatchTime=0.01E-6;
  Double_t sigma=0;
  fsigma=sigma;
  if(sigma==0) fsigma=5;
  derivativeSigma=2.5;
  windowSize=7;
  //windowSize=2;
  nSigma=5;
  aveWidth=20;
  spec = new TSpectrum();
  //Int_t irunStop = irunStart;
  

  //Int_t irunStop = irunStart;
  TString outFileName2 ; outFileName2.Form("TBacon-anaCompass-%s.root",tag.Data(),maxEvents,derivativeSigma);
  TFile *outfile2 = new TFile(outFileName2,"recreate");
  outfile2->cd();
  TBacon =  new TTree("TBacon"," bacon data ");
  baconEvent = new TBaconEvent();
  TBacon->Branch("bevent",&baconEvent);


  TString outFileName ; outFileName.Form("anaCompass-%s.root",tag.Data(),maxEvents,derivativeSigma);
  TFile *outfile = new TFile(outFileName,"recreate");
  outfile->mkdir("pulses");
  hLatePulse = new TH1D("LatePulse","late pulse",2*PWIDTH,0,2*PWIDTH);
  hPromptPulse = new TH1D("PromptPulse","prompt pulse",2*PWIDTH,0,2*PWIDTH);
  hLateNoisePulse = new TH1D("LateNoisePulse","late noise pulse",2*PWIDTH,0,2*PWIDTH);
  outfile->cd();

  printf(" opening output file %s and %s \n",outFileName.Data(),outFileName2.Data());

   
  //ntBase = new TNtuple("ntBase","base","iw:w:b:bnon:bneil:width");
  TNtuple *ntMatchTime = new TNtuple("ntMatchTime"," match time ","ihit:dt");
  ntCal =  new TNtuple("ntCal","ntuple Cal","iev:ipmt:base:sigma:dbase:dsigma");
  // now a TTree
  treeHit = new TTree("THit","Tree with hits");
  treeHit->Branch("hits",&theHit,"nev/I:nhits/I:order/I:istart/I:nwidth/I:good/I:kind/I:time/F:thit/F:q/F:qerr/F:peak/F:qsum/F:q900/F");


  //
  ntNHit = new TNtuple("ntNHit","negative ntuple Hit","npmt:nhits:order:istart:time:q:nwidth:peak:good:kind");
  ntDer =  new TNtuple("ntDer"," deriviative ","t:sigma:d0:kover:type");
  ntEvent= new TNtuple("ntEvent","ntuple Event","entry:n0:nspe0:wsum0:wsum1:qsum:qsum1:qsum3:qspe0:q900");
 // ntPulse= new TNtuple("ntPulse"," pulse ","sum:shigh:slow:nsamp:kover:qlow:qhigggh:klow:khigh");
  ntWave = new TNtuple("ntWave"," wave ","event:v:d");

  hHitLength = new TH1I("HitLength"," hit length",100,0,100);
  hPeakNWidth = new TH1I("PeakNWidth","PeakNWidth",100,0,100);
  hAllNWidth = new TH1I("AllNWidth","AllNWidth",100,0,100);
  hNWidthCut = new TH1I("NWidthCut","NWidth cut ",100,0,100);
  hDerWidth = new TH1I("DerWidth","sliding window pulse number Width",2*maxHalfLength,0,2*maxHalfLength);
  hDerAfter = new TH1D("DerAfter"," after cut derivative",1000,-.004,.004);
  hSlideHigh = new TH1D("SlideHigh"," sliding window high threshold ",1000,0,1);
  hSlideLow  = new TH1D("SlideLow"," sliding window low threshold",1000,0,.1);
  hSlideQSum = new TH1D("SlideQSum"," sliding window Q sum",100,0,100);

  hQFirst = new TH2D("QFirst"," q versus first time ",100,0,4,1000,0,10);
  hQFirst->GetXaxis()->SetTitle(" micro-seconds ");
  hQFirst->GetYaxis()->SetTitle(" hit charge ");

  hQStart = new TH2D("QStart"," q versus time  ",200,-4,4,1000,0,10);
  hQStart->GetXaxis()->SetTitle(" micro-seconds from first hit ");
  hQStart->GetYaxis()->SetTitle(" hit charge ");

  hNegQStart = new TH2D("NegQStart"," negative pulse q versus time  ",200,-4,4,1000,0,10);
  hNegQStart->GetXaxis()->SetTitle(" micro-seconds from first hit ");
  hNegQStart->GetYaxis()->SetTitle(" hit charge ");



  // define pmt signal histograms
  if(isSimulation) { printf(" \t\t this is simulation \n");
    simMatchStats = new TPmtSimMatchStats();
  }
  else printf(" \t\t this is real data \n");


  //switch to output file
  outfile->cd();

  //Event Vectors
  Int_t movingN1 = 0;
  Int_t movingN2 = 0;
  Double_t runningNoise1 = 0;
  Double_t runningNoise2 = 0;
  std::vector<Double_t> runningNoise(2);
  std::vector<Int_t> movingN(2);
  //cout << " XXXXXXXX movingN size " << movingN.size() << endl;
  //printf(" peak finding parameters: \n");
  //printf(" \t sigma %.2f minLength %i  \n",fsigma, minLength); 

  double maxSample[2];
  Float_t qpmt[NPMT][MAXHIT];
  Float_t tpmt[NPMT][MAXHIT];
  Int_t   npmtHit[NPMT];
  Int_t   nspe[NPMT];
  Double_t qspe[NPMT];
  Double_t q900[NPMT];
  Double_t wsum[NPMT];
  Double_t qsum[NPMT];
  Double_t qsum1[NPMT];
  Double_t qsum3[NPMT];
  Double_t qped[NPMT];
  Int_t nped[NPMT];
  Int_t totalHits[NPMT]={0,0};


  // loop over entries
  printf(" STARTING RUN %s with  events  %lld of %lld derivative smoothing = %i \n\n",tag.Data(),maxEvents,nentries,windowSize);
  if(maxEvents>0) nentries=maxEvents;
  for(Long64_t entry =0; entry< nentries ; ++ entry) {
    baconEvent->clear();
    Long64_t nb=getEvent(entry);
    nSamples = Samples->GetSize();

    baconEvent->event     = entry;
    baconEvent->channel   = Int_t(Channel);
    baconEvent->board     = Int_t(Board);
    baconEvent->flags     = Int_t(Flags);
    baconEvent->timeStamp = Long64_t(Timestamp);
    baconEvent->energy    = Double_t(Energy);



    if(entry%1000==0) printf(" .... %lld %lld %i \n",entry,nb,nSamples);
    for(int ipmt=0; ipmt<NPMT ; ++ ipmt) {
      npmtHit[ipmt]=0;
      qped[ipmt]=0;
      qsum[ipmt]=0;
      wsum[ipmt]=0;
      qsum1[ipmt]=0;
      qsum3[ipmt]=0;
      nped[ipmt]=0;
      nspe[ipmt]=0;
      qspe[ipmt]=0;
      q900[ipmt]=0;
      for(int ihit=0; ihit<MAXHIT ; ++ihit) {
        qpmt[ipmt][ihit]=0;
        tpmt[ipmt][ihit]=0;
      }
    }

    // check how many pmts we have
    int gotPMT = 1;
    
    TString name; name.Form("%s_Ev%i",tag.Data(),entry);
    // define pmt signal histograms
    if(entry==0) {
      source = new Double_t[nSamples];
      Double_t maxLife = 10.*microSec;
      timeUnit=2;
      printf(" \n\n ***** setting time unit %E maxLife %f # digis %lu  \n",timeUnit,maxLife,Samples->GetSize());
      if(isSimulation) printf(" \n\n ***** setting matching time  %E \n",simHitMatchTime);
      Double_t pmtXLow= microSec;
      Double_t pmtXHigh= 100000*microSec;
      if(isSimulation) {
        hIsSimHitMatchTime = new TH1D("IsSimHitMatch","",1000,-10*simHitMatchTime,10*simHitMatchTime);
        hAllSimHitMatchTime = new TH1D("AllSimHitMatch","",1000,-10*simHitMatchTime,10*simHitMatchTime);
        hSimHitMatched = new TH1D("SimHitMatched","",nSamples,pmtXLow,pmtXHigh);
        hSimHitNotMatched = new TH1D("SimHitNotMatched","",nSamples,pmtXLow,pmtXHigh);
        hSimHitMissed = new TH1D("SimHitMissed","",nSamples,pmtXLow,pmtXHigh);
        ntSimMatch = new TNtuple("ntSimMatch"," hit matches ","pmt:nsim:nhit:match:nnot:nmiss");
      }
      for(int ipmt=0; ipmt<gotPMT; ++ipmt) {

        hSPE[ipmt] = new TH1D(Form("SPE%i",ipmt),Form(" SPE charge PMT %i ",ipmt),500,0,.5);
        hSPE[ipmt]->GetXaxis()->SetTitle(" charge*bin (ADC-ns) ");


        hQSum[ipmt] = new TH1D(Form("QSum%i",ipmt),Form(" sum charge PMT %i ",ipmt),1500,0,10);
        hQSum[ipmt]->GetXaxis()->SetTitle(" charge ");

        hQEarly[ipmt] = new TH1D(Form("QEarly%i",ipmt),Form(" sum charge PMT %i ",ipmt),1500,0,10);
        hQEarly[ipmt]->GetXaxis()->SetTitle(" charge ");

        hQLate[ipmt] = new TH1D(Form("QLate%i",ipmt),Form(" sum charge PMT %i ",ipmt),1500,0,10);
        hQLate[ipmt]->GetXaxis()->SetTitle(" charge ");

        hNoise[ipmt] = new TH1D(Form("Noise%i_%s",ipmt,tag.Data()),"",1000,0,.1);

        hPMTRaw[ipmt] = new TH1D(Form("PMTRaw%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hPMTSignal[ipmt] = new TH1D(Form("PMTSignal%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hPMTDerivative[ipmt] = new TH1D(Form("PMTDeriv%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hDerivativeNoise[ipmt] = new TH1D(Form("DerivNoise%i_%s",ipmt,tag.Data()),"",1000,-.02,.02);
        hSum[ipmt] = new TH1D(Form("SumPmt%i",ipmt),"",nSamples,pmtXLow,pmtXHigh);
        hSum[ipmt]->GetXaxis()->SetTitle(" seconds ");
        hBaseline[ipmt]  = new TH1D(Form("Baseline%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hBaselineFit[ipmt]  = new TH1D(Form("BaselineFit%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hBaselineWMA[ipmt]  = new TH1D(Form("BaselineWMA%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hWeight[ipmt]  = new TH1D(Form("BaselineWeight%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        hNHits[ipmt] = new TH1D(Form("NHits%i",ipmt),Form(" number of hits PMT %i ",ipmt),10,0,50);
        hNegNHits[ipmt] = new TH1D(Form("NegNHits%i",ipmt),Form(" number of neg hits PMT %i ",ipmt),10,0,50);
        hLife[ipmt] = new TH1D(Form("Life%i",ipmt),Form(" lifetime PMT %i ",ipmt),500,0,maxLife);
        hLife[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
        hNLife[ipmt] = new TH1D(Form("NLife%i",ipmt),Form(" negative pulse lifetime PMT %i ",ipmt),500,0,maxLife);
        hNLife[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
        //
        hLifeCut[ipmt] = new TH1D(Form("LifeCut%i",ipmt),Form(" Q lifetime PMT %i ",ipmt),500,0,maxLife);
        hLifeCut[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");

        hLifeCount[ipmt] = new TH1D(Form("LifeCount%i",ipmt),Form(" Start lifetime PMT %i ",ipmt),500,0,maxLife);
        hLifeCount[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");


        hNLifeCut[ipmt] = new TH1D(Form("NLifeCut%i",ipmt),Form(" negative pulse lifetime PMT %i ",ipmt),500,0,maxLife);
        hNLifeCut[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
        hLife[ipmt]->Sumw2();
        hLifeCut[ipmt]->Sumw2();
        hNLife[ipmt]->Sumw2();

        //hLife[ipmt]->Sumw2();
        if(isSimulation) {
          hLifeSim[ipmt] = new TH1D(Form("LifeSim%i",ipmt),Form(" real pulse lifetime PMT %i ",ipmt),500,0,maxLife);
          hLifeSim[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
          hLifeSim[ipmt]->Sumw2();

          hLifeTrue[ipmt] = new TH1D(Form("LifeTrue%i",ipmt),Form(" truth lifetime PMT %i ",ipmt),500,0,maxLife);
          hLifeTrue[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
          hLifeTrue[ipmt]->Sumw2();
          
          hLifeNoise[ipmt] = new TH1D(Form("LifeNoise%i",ipmt),Form(" noise pulse lifetime PMT %i ",ipmt),500,0,maxLife);
          hLifeNoise[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
          hLifeNoise[ipmt]->Sumw2();

          hPMTSim[ipmt] = new TH1D(Form("PMTSim%i_%s",ipmt,tag.Data()),"",nSamples,pmtXLow,pmtXHigh);
        }      
      }
      // initialize fft 
      fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
      fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
    }


    //outfile->ls();

    std::vector<std::vector<Double_t> > ddigi; ddigi.resize(NPMT);
    std::vector<std::vector<Double_t> > ndigi; ndigi.resize(NPMT);
    std::vector<std::vector<Double_t> > deriv;  deriv.resize(NPMT);
    std::vector<std::vector<Double_t> > nderiv; nderiv.resize(NPMT);
    std::vector<std::vector<Double_t> > sddigi; sddigi.resize(NPMT); // baseline subtracted
    std::vector<std::vector<Double_t> > sndigi; sndigi.resize(NPMT);

    
    std::vector<std::vector<Double_t> > baselineDigi; baselineDigi.resize(NPMT);

    Double_t wbase[2];
    wbase[0]=0; wbase[1]=0;

    unsigned uEndSample=Samples->GetSize()/10;
    if(entry==0) printf(" XXXXX simple baseline starts  %u at Samplesle %u \n",uEndSample,unsigned(Samples->GetSize())-uEndSample);
    // calculate simple baseline

    for(unsigned is = 0; is < Samples->GetSize(); is++) {
      if(is> Samples->GetSize() - uEndSample) wbase[0]+= double(Samples->GetAt(is));
      hPMTRaw[0]->SetBinContent(is+1,double(Samples->GetAt(is)));
    }
    wbase[0] /= Double_t(uEndSample); 


    // loop over PMT 1
    for(unsigned is = 0; is < Samples->GetSize(); is++){
      Double_t digi0 = -( double(Samples->GetAt(is))- wbase[0] );
      //printf(" %lld %is %f \n",entry,is,digi0);
      ddigi[0].push_back(digi0);
      ndigi[0].push_back(-digi0);
    }
      
    //FFT
    if(entry<1) {
       for(int ipmt=0; ipmt<2; ++ipmt) if(ddigi[ipmt].size()>0) FFT(ipmt,entry,ddigi[ipmt]);
    }

  
    // derivative PMT
    Double_t derAve[NPMT], derSigma[NPMT];
    Double_t rawAve[NPMT], rawSigma[NPMT];
    for(int j = 0; j < gotPMT; j++){
      getAverage(ddigi[j],rawAve[j],rawSigma[j]);
      deriv[j] = differentiate( ddigi[j],windowSize);
      nderiv[j] = differentiate( ndigi[j],windowSize);
      getAverage(deriv[j],derAve[j],derSigma[j]);
      for(unsigned isample = 0; isample < deriv[j].size(); isample++){
        hPMTDerivative[j]->SetBinContent(isample+1,deriv[j][isample]);
        ntWave->Fill(entry,ddigi[j][isample],deriv[j][isample]);
      }
      ntCal->Fill(entry,j,rawAve[j],rawSigma[j],derAve[j],derSigma[j]);
      for(unsigned isample = 0; isample < deriv[j].size(); isample++) hDerivativeNoise[j]->Fill(deriv[j][isample]-derAve[j]);
    }

    
    for(int pmtNum = 0 ; pmtNum < gotPMT; pmtNum++){
      /* 
      ** peak finding 
      */
      std::vector<Int_t> peakKind;
      peakType peakList = derivativePeaks(deriv[pmtNum],windowSize,derSigma[pmtNum],peakKind);
      if(peakList.size()<1) { 
        ++noPeakEventCount; continue; 
      } 

      // get baseline
      Int_t maxwidth=0;
      std::vector<Double_t> weight = getBaselineWeights(ddigi[pmtNum].size(),peakList,maxwidth);
      std::vector<Int_t> npeakKind;
      peakType npeakList = derivativePeaks(nderiv[pmtNum],windowSize,derSigma[pmtNum],npeakKind);
      baselineDigi[pmtNum] = getBaselineWMARecursive(rawAve[pmtNum],ddigi[pmtNum],weight,unsigned(2*maxHalfLength));

     // subtract baseline
      maxSample[pmtNum]=0;
      hPMTSignal[pmtNum]->Reset();
      hBaselineWMA[pmtNum]->Reset();
      if(isSimulation) hPMTSim[pmtNum]->Reset();
      for(int i = 0; i < int(ddigi[pmtNum].size()); i++){
        if(ddigi[pmtNum][i]>maxSample[pmtNum])  maxSample[pmtNum] = ddigi[pmtNum][i];
        hBaselineWMA[pmtNum]->SetBinContent(i,baselineDigi[pmtNum][i]);
        sddigi[pmtNum].push_back(ddigi[pmtNum][i]-baselineDigi[pmtNum][i]);
        sndigi[pmtNum].push_back(ddigi[pmtNum][i]-baselineDigi[pmtNum][i]);
        hPMTSignal[pmtNum]->SetBinContent(i,(ddigi[pmtNum][i]-baselineDigi[pmtNum][i]));
      }
      getAverage(sddigi[pmtNum],sAve[pmtNum],sDev[pmtNum]);
      hNoise[pmtNum]->Fill(sDev[pmtNum]);
      

      // trim peaks
      trimPeaks(peakList,sddigi[pmtNum]);
      trimPeaks(npeakList,sndigi[pmtNum]);

      Double_t minDev = 0*sDev[pmtNum];
      Double_t maxDev = fsigma*sDev[pmtNum];
      Double_t triggerTime, firstCharge;
      hitMap  pmtHits = makeHits(peakList,peakKind,sddigi[pmtNum],sDev[pmtNum],triggerTime,firstCharge);
      // negative pulses
      Double_t ntriggerTime, nfirstCharge;
      hitMap  npmtHits = makeHits(npeakList,npeakKind, ndigi[pmtNum],sDev[pmtNum],ntriggerTime,nfirstCharge);

      for(unsigned iw=0; iw< weight.size(); ++iw) {
        //ntBase->Fill(float(iw),weight[iw],baserec[iw],basenon[iw],baselineDigi[pmtNum][iw],float(maxwidth));
        hWeight[pmtNum]->SetBinContent(iw,weight[iw]);
      }

      if(entry%printInterval==0) printf(" \t  %i ave %.4E sdev0 %.4E dave %.4E dsigma %.4E \n",
          entry,baseline[pmtNum],sDev[pmtNum],derAve[pmtNum],derSigma[pmtNum]);

      unsigned nhits = pmtHits.size();
      hNHits[pmtNum]->Fill(nhits);
      unsigned negnhits = npmtHits.size();
      hNegNHits[pmtNum]->Fill(negnhits);

 
      if(entry%printInterval==0) printf(" \t pmt  %i peaktime %lu nhits %u \n ",pmtNum,peakList.size(),nhits);

      totalHits[pmtNum] += nhits;

      if(nhits<1) {
        // pedestal pulses 
        for(int i = 0; i < int(ddigi[pmtNum].size()); i++){
          qped[pmtNum] += ddigi[pmtNum][i];
          ++nped[pmtNum];
        }
        if(nped[pmtNum]>0) qped[pmtNum] /= 20.0*Double_t(nped[pmtNum]);
        continue;
      }

     
      npmtHit[pmtNum]=nhits;
      hitMapIter hitIter;
      hitIter=pmtHits.begin();
      TPmtHit phit0 = hitIter->second;
      TPmtHit phit1;
      // next hit if it is there
      if(pmtHits.size()>1) {
        ++hitIter;
        phit1 = hitIter->second;
      }
      // pulse width cut
      Int_t nwidth0 = phit0.lastBin - phit0.firstBin +1;

      // peak 0 width cut 
      hNWidthCut->Fill(nwidth0);
      //if(nwidth0>20) continue;


      hQFirst->Fill(triggerTime,firstCharge);
      //if(firstCharge<firstChargeCut&&phit0.qsum>firstChargeCut) 
      if(triggerTime==1E9&&phit0.qsum>firstChargeCut) 
        printf(" WARNING NO FIRST PULSE event %i  pmt %i pulses %i qhit %f time %f first %f charge %f \n",
            entry,pmtNum,nhits,phit0.qsum,phit0.startTime*microSec,triggerTime,firstCharge);

   
      /* ******compare hits and simulation*******
      */
      if(isSimulation) {
        unsigned nnot=0;
        std::vector<Double_t> startTime=pmtSimulation->startTime;
        std::vector<Int_t> hitMatchNumber(pmtHits.size(),-1); 
        if(debug) printf(" event %u photons %lu hits %lu \n",entry,startTime.size(),pmtHits.size());
        // loop over hits
        unsigned ihit=0;
        for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
          TPmtHit phiti = hitIter->second;
          // loop over sim 
          for(unsigned isim = 0 ; isim < startTime.size(); ++ isim ) {
            Double_t tdiff =  phiti.peakt*timeUnit-startTime[isim];
            hAllSimHitMatchTime->Fill(tdiff);
            ntMatchTime->Fill(ihit,tdiff);
            if(ihit==0) { // only fill once per event
              hPMTSim[pmtNum]->Fill(startTime[isim]*microSec);
              hLifeTrue[pmtNum]->Fill(startTime[isim]*microSec);
            }
            if(abs(tdiff)<simHitMatchTime) {
              ++hitIter->second.good;
              hitMatchNumber[ihit]=isim;
              hIsSimHitMatchTime->Fill(tdiff);
            }
          }
          // count not matched (noise hits)
          if( hitIter->second.good==0) {
            ++nnot;
            hSimHitNotMatched->Fill(phiti.peakt*timeUnit*microSec);
          }
          ++ihit;
        }
        // count matches
        unsigned nmatch=0;
        unsigned nmiss=0;
        for(unsigned isim=0; isim<startTime.size(); ++isim) {
          bool isMatched=false;
          for(unsigned ihit=0; ihit<hitMatchNumber.size(); ++ihit) if( hitMatchNumber[ihit]==isim ) isMatched=true;
          if(isMatched){
            ++nmatch;
            hSimHitMatched->Fill(startTime[isim]*microSec);
          } else {
            ++nmiss;
            hSimHitMissed->Fill(startTime[isim]*microSec);
          }
        }
       
        ntSimMatch->Fill(float(pmtNum),float(startTime.size()),float(pmtHits.size()),float(nmatch),float(nnot),float(nmiss));
        if(pmtNum==0) simMatchStats->fill(startTime.size(),pmtHits.size(),nmatch,nnot,nmiss);
        if(entry%printInterval==0) simMatchStats->print();
        //simMatchStats->print();
      } // end is sim

      if(nHists<nMaxHistEvents) {
        plotWave(entry,pmtNum,pmtHits );
        ++nHists;
        // print hits
        if(debug) printf(" event %u hits %lu : \n",entry,pmtHits.size());
        unsigned hitNumber=0;
         for (hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
            TPmtHit phiti = hitIter->second;
            if(debug) printf("\t %3u (%5i,%5i) kind=%i q=%f \n",hitNumber++,phiti.firstBin,phiti.lastBin,phiti.kind,phiti.qsum);
        }
        
      }

      /* sum hits */
      for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
        TPmtHit phiti = hitIter->second;
        Double_t phitTime =  phiti.startTime*microSec-triggerTime;
        Double_t phitQ = phiti.qsum*timeUnit*1E9;
        hQStart->Fill(phitTime,phitQ);
        qsum[pmtNum] += phitQ;
        if(phitTime<0.9) q900[pmtNum] += phitQ;
        if(phitTime>speTimeCut) {
          hSPE[pmtNum]->Fill(phitQ);
          ++nspe[pmtNum];
          qspe[pmtNum] += phitQ;
          qsum3[pmtNum] += phitQ;
        } else {
          qsum1[pmtNum] += phitQ;
        }
      }

      // look at hits
      int hitCount=0;
      if(debug) printf(" event %u \n",entry);
      for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
        TPmtHit phiti = hitIter->second;
        Double_t phitTime =  phiti.startTime*microSec-triggerTime;
        Double_t phitQ = phiti.qsum*timeUnit*1E9;
        Double_t phitQErr = phiti.qerr*timeUnit*1E9;
        if(debug) printf(" hit %i time %f q %f good %i \n  ",hitCount, phiti.startTime*microSec, phitQ,phiti.good);
        hQStart->Fill(phitTime,phitQ);
        Int_t nwidth = phiti.lastBin - phiti.firstBin +1;

        Int_t istartBin =  hLife[pmtNum]->FindBin(phitTime); 
        //fill the hit
        theHit.nev=int(entry);
        theHit.nhits=nhits;
        theHit.order=hitCount;
        theHit.istart=istartBin;
        theHit.nwidth=nwidth;
        theHit.good=phiti.good;
        theHit.kind=phiti.kind;
        theHit.time=phiti.startTime*microSec;
        theHit.thit=phitTime;
        theHit.q=phitQ;
        theHit.qerr=phitQErr;
        theHit.peak=phiti.qpeak;
        theHit.qsum=qsum[0];
        theHit.q900=q900[0];
        if(pmtNum==0) {
          TPulse thePulse;
          thePulse.istart=istartBin;
          thePulse.nwidth=nwidth;
          thePulse.good=phiti.good;
          thePulse.kind=phiti.kind;
          thePulse.time=phiti.startTime*microSec;;
          thePulse.thit=phitTime;
          thePulse.q=phitQ;
          thePulse.qerr=phitQErr;
          if(debug) printf(" pulse %i %i %i %i %f %f %f %f  \n", istartBin, nwidth, phiti.good, phiti.kind, thePulse.time, phitTime,phitQ, phitQErr);
          baconEvent->hits.push_back(thePulse);
        }
        treeHit->Fill();

        hLife[pmtNum]->SetBinContent( istartBin, hLife[pmtNum]->GetBinContent(istartBin)+phitQ);
        hLife[pmtNum]->SetBinError(istartBin, sqrt( pow(hLife[pmtNum]->GetBinError(istartBin),2)+pow(phiti.qerr,2) ));
        //if(phitQ>lifeChargeCut) {
        if(phiti.kind==0) {
          hLifeCut[pmtNum]->SetBinContent( istartBin, hLife[pmtNum]->GetBinContent(istartBin)+phitQ);
          hLifeCut[pmtNum]->SetBinError(istartBin, sqrt( pow(hLifeCut[pmtNum]->GetBinError(istartBin),2)+pow(phiti.qerr,2) ));
          hLifeCount[pmtNum]->Fill( phitTime );
        }

        if(isSimulation) {
          if( phiti.good>0 ) {
            hLifeSim[pmtNum]->SetBinContent( istartBin, hLifeSim[pmtNum]->GetBinContent(istartBin)+phitQ);
            hLifeSim[pmtNum]->SetBinError(istartBin, sqrt( pow(hLifeSim[pmtNum]->GetBinError(istartBin),2)+pow(phiti.qerr,2) ));

          }
          else {
            hLifeNoise[pmtNum]->SetBinContent( istartBin, hLifeNoise[pmtNum]->GetBinContent(istartBin)+phitQ);
            hLifeNoise[pmtNum]->SetBinError(istartBin, sqrt( pow(hLifeNoise[pmtNum]->GetBinError(istartBin),2)+pow(phiti.qerr,2) ));

          }
        }

        // plot  pulses
        std::vector<Double_t> pulsei = phiti.getPulse(PWIDTH,sddigi[pmtNum]);
        if( pulseHistCount++ < nMaxPulseHist) {
          outfile->cd("pulses");
          TH1D* phist = new TH1D(Form("Pulse%i_Ev%i_PMT%i",hitCount,entry,pmtNum),Form("Pulse%i_Ev%i_PMT%i tine %.2f ",hitCount,entry,pmtNum, phitTime), pulsei.size(),0,pulsei.size());
         if(debug) printf(" plot puleses   %.2f size %lu \n ",phitTime,pulsei.size() );
          for(unsigned ibin=0; ibin<pulsei.size(); ++ibin ) phist->SetBinContent(ibin, pulsei[ibin]);
          outfile->cd();
        }
        if(phitTime>2.0) {
          for(unsigned ibin=0; ibin<pulsei.size(); ++ibin ) {
            if(isSimulation) {
              if( phiti.good>0 ) hLatePulse->SetBinContent( ibin, hLatePulse->GetBinContent(ibin)+pulsei[ibin]) ;
              else hLateNoisePulse->SetBinContent( ibin, hLateNoisePulse->GetBinContent(ibin)+pulsei[ibin]) ;
            } else 
              hLatePulse->SetBinContent( ibin, hLatePulse->GetBinContent(ibin)+pulsei[ibin]) ;
          }
        } else 
          for(unsigned ibin=0; ibin<pulsei.size(); ++ibin ) hPromptPulse->SetBinContent( ibin, hPromptPulse->GetBinContent(ibin)+pulsei[ibin]) ;

        ++hitCount;
      } 
      if(debug) printf(" ////// \n\n");

      // summed wave forms only if passes first charge cut
      if(firstCharge>firstChargeCut) sumWave(pmtNum);

      // negative pulse ntuple
      // hitMapIter hitIter;
      hitMapIter nhitIter;
      nhitIter=npmtHits.begin();
      TPmtHit nphit0 = nhitIter->second;

      // look at negative hits
      hitCount=0;
      for (hitMapIter nhitIter=npmtHits.begin(); nhitIter!=npmtHits.end(); ++nhitIter) {
        TPmtHit phiti = nhitIter->second;
        Double_t phitQ = phiti.qsum*timeUnit*1E9;
        //Double_t phitQErr = phiti.qerr*timeUnit*1E9;
        Int_t nwidth = phiti.lastBin - phiti.firstBin +1;
        Double_t phitTime =  phiti.startTime*microSec-ntriggerTime;
        hNegQStart->Fill(phitTime,phitQ);
        Int_t istartBin =  hNLife[pmtNum]->FindBin(phitTime);
        ntNHit->Fill(pmtNum,nhits,hitCount,istartBin, phiti.startTime*microSec,phitQ,nwidth,phiti.qpeak,phiti.good,phiti.kind);
        hNLife[pmtNum]->SetBinContent( istartBin, hNLife[pmtNum]->GetBinContent(istartBin)+phitQ);
        if(phitQ>lifeChargeCut) hNLifeCut[pmtNum]->SetBinContent( istartBin, hNLife[pmtNum]->GetBinContent(istartBin)+phitQ);
      }
      if(ntriggerTime==1E9&&nphit0.qsum>firstChargeCut) 
        printf(" WARNING NO NEGATIVE FIRST PULSE event %i  pmt %i pulses %i qhit %f time %f first %f charge %f \n",
            entry,pmtNum,negnhits,nphit0.qsum,nphit0.startTime*microSec,ntriggerTime,nfirstCharge);

      // save for npmtHit ntuple
      qpmt[pmtNum][0]=phit0.qsum;
      tpmt[pmtNum][0]=phit0.startTime*microSec;    // milli-sceconds 
      qpmt[pmtNum][1]=phit1.qsum;
      tpmt[pmtNum][1]=phit1.startTime*microSec;

      // pedestal pulses 
      for(int i = 0; i < int(ddigi[pmtNum].size()); i++){
        // remove pulses 
        bool skip = false;
        for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
          TPmtHit phiti = hitIter->second;
          if( i >= phiti.firstBin && i <= phiti.lastBin ) skip=true;
        }
        if(!skip) {
          qped[pmtNum] += ddigi[pmtNum][i];
          ++nped[pmtNum];
        }
      }
      if(nped[pmtNum]>0) qped[pmtNum] /= Double_t(aveWidth*nped[pmtNum]);

      }

      //cout << ddigi[0].size() << "   " <<  nped[0]  << "  " << nped[1] << " " << qped[0] << "  " << qped[1] <<  endl;

      hQSum[0]->Fill(qsum[0]);
      hQEarly[0]->Fill(qsum1[0]);
      hQLate[0]->Fill(qsum3[0]);
      //if(npmtHit[0]>0||npmtHit[1]>0) 
      ntEvent->Fill(entry,npmtHit[0],nspe[0],wsum[0],wsum[1],qsum[0],qsum1[0],qsum3[0],qspe[0],q900[0]);
      baconEvent->event=entry;
      baconEvent->npulse=baconEvent->hits.size();
      baconEvent->nspe=nspe[0];
      baconEvent->wsum=wsum[0];
      baconEvent->qsum=qsum[0];
      baconEvent->q900=q900[0];
      TBacon->Fill();
      if(entry%printInterval==0) printf(" total hits PMT1 %i PMT2 %i out of %u events events with no peaks %u \n",totalHits[0],totalHits[1],entry,noPeakEventCount);
    }
    printf(" total hits PMT1 %i PMT2 %i out of %lld events events with no peaks %u \n",totalHits[0],totalHits[1],nentries,noPeakEventCount);

    printf(" derivative sigma is %.2f \n",derivativeSigma);
  if(isSimulation) simMatchStats->print();

  TH1D* htimeCountNorm = (TH1D*) hLifeCount[0]->Clone("timeCountNorm");
  htimeCountNorm->SetTitle(" average photon time (ns) (PMT 0) ");
  if(nentries>0) for(int ibin=0; ibin<htimeCountNorm->GetNbinsX()+2; ++ibin) htimeCountNorm->SetBinContent(ibin, hLifeCount[0]->GetBinContent(ibin)/float(nentries));

  TH1D* hQSumNorm = (TH1D*) hQSum[0]->Clone("QSumNorm");
  hQSumNorm->SetTitle(" average Q  (PMT 0) ");
  if(nentries>0) for(int ibin=0; ibin<hQSumNorm->GetNbinsX()+2; ++ibin) hQSumNorm->SetBinContent(ibin, hQSum[0]->GetBinContent(ibin)/float(nentries));

  TH1D* hQEarlyNorm = (TH1D*) hQEarly[0]->Clone("QEarlyNorm");
  hQEarlyNorm->SetTitle(" average prompt Q  (PMT 0) ");
  if(nentries>0) for(int ibin=0; ibin<hQEarlyNorm->GetNbinsX()+2; ++ibin) hQEarlyNorm->SetBinContent(ibin, hQEarly[0]->GetBinContent(ibin)/float(nentries));
 
  TH1D* hQLateNorm = (TH1D*) hQLate[0]->Clone("QLateNorm");
  hQLateNorm->SetTitle(" average late Q  (PMT 0) ");
  if(nentries>0) for(int ibin=0; ibin<hQLateNorm->GetNbinsX()+2; ++ibin) hQLateNorm->SetBinContent(ibin, hQLate[0]->GetBinContent(ibin)/float(nentries));

  // find SPE
  TF1 *fLan = new TF1("flan","landau",speStartTimeCut,1);
  hSPE[0]->Fit("flan","R");
  Double_t SPEval= fLan->GetParameter(1);

  //
  hSPE[0]->SetTitle(Form("%s SPE peak %.3f",tag.Data(),SPEval));
  //
  Double_t meanQ=hQSumNorm->GetMean();
  Double_t meanPE=0;
  if(SPEval>0) meanPE = meanQ/SPEval;
  hQSumNorm->SetTitle(Form(" mean Q = %.3f SPE %.3f PMT 0",meanQ , meanPE));
  //
  Double_t meanEarlyQ=hQEarlyNorm->GetMean();
  Double_t meanEarlyPE=0;
  if(SPEval>0) meanEarlyPE = meanEarlyQ/SPEval;
  hQEarlyNorm->SetTitle(Form(" mean Q = %.3f SPE %.3f PMT 0",meanEarlyQ , meanEarlyPE));
  //
  Double_t meanLateQ=hQLateNorm->GetMean();
  Double_t meanLatePE=0;
  if(SPEval>0)  meanLateQ/SPEval;
  hQLateNorm->SetTitle(Form(" mean Q = %.3f SPE %.3f PMT 0",meanLateQ , meanLatePE));
  //
   printf("REPORT %s # events %lld SPE peak is %.3f  <Q> %.3f prompt %.3f late %.3f Q-SPE/event %.3f prompt %.3f late %.3f  \n",
      tag.Data(),nentries,SPEval,meanQ,meanEarlyQ,meanLateQ,meanPE,meanEarlyPE,meanLatePE);
  outfile->Purge();
  outfile->Write();
  outfile2->Write();
  return;
  }


  hitMap anaCompass::makeHits(peakType peakList, std::vector<Int_t> peakKind, std::vector<Double_t> ddigi,Double_t sigma, Double_t& triggerTime, Double_t& firstCharge) 
  {

    triggerTime=1E9;
    firstCharge=0;
  hitMap pmtHits;
  if(peakList.size()<1) return pmtHits;
  Double_t qmax=0;
 
  unsigned minLength=5;
  for(unsigned ip=0; ip<peakList.size(); ++ip) {
    unsigned klow  = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);
    //printf(" hit  %u (%u,%u) kind %i length %u \n",ip,klow,khigh,peakKind[ip],khigh-klow+1);
    hHitLength->Fill(khigh-klow+1);
    if( khigh-klow+1<minLength) {
      continue;
    }
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    Double_t qsum = 0;
    for(unsigned k=klow; k<khigh; ++k) {
      qsum += ddigi[k];
      if(ddigi[k]>qpeak) {
        peakt=k;
        qpeak = ddigi[k];
      }
    }

    TPmtHit phit;
    phit.peakBin=Int_t(peakt);
    phit.qsum=qsum;
    phit.qpeak=qpeak;
    phit.firstBin = klow;
    phit.lastBin = khigh;
    phit.peakMaxTime=peakt;
    phit.peakt=peakt;
    phit.startTime=klow;
    phit.peakWidth=khigh - klow;
    // this is N= q/qnorm and delta q = root(n)*qnorm;
    phit.qerr = sqrt(  pow(sigma*Double_t(phit.peakWidth),2)+ qnorm*qsum );
    phit.kind = peakKind[ip];

    // just use the biggest pulse 
    if(qsum>qmax) {
      qmax=qsum;
      triggerTime=phit.startTime*microSec;
      firstCharge = qsum;
    }
    Double_t hitTime = phit.startTime*microSec;
    //printf(" insert  %lu (%u,%u) %f \n",pmtHits.size(),phit.firstBin,phit.lastBin,phit.startTime*microSec);
    pmtHits.insert ( std::pair<Double_t,TPmtHit>(hitTime,phit) );
    hPeakNWidth->Fill(phit.lastBin-phit.firstBin+1);
  }

  // first time, charge from map
  /*
  hitMapIter hitIter;
  hitIter=pmtHits.begin();
  TPmtHit phit0 = hitIter->second;
  triggerTime = phit0.startTime*microSec;
  firstCharge = phit0.qsum;
  */
  triggerTime=0.0;
  if(firstCharge<firstChargeCut&&pmtHits.size()>0&&qmax>firstChargeCut) printf("\t WARNING XXXXX NO FIRST PULSE pulses %i max %f \n",int(pmtHits.size()),qmax);

  return  pmtHits;
}

std::vector<std::complex<double> > anaCompass::FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal)
{
  std::vector<std::complex<double> > VectorComplex;
  for(int is =0; is<nSamples; ++is) {
    if(ipmt==0) fFFT->SetPoint(is, signal[is]);
    else fFFT->SetPoint(is, signal[is]);
  }
  fFFT->Transform();
  TH1D* hfft = new TH1D(Form("FFTPMT%i_Ev%i",ipmt,ievent),Form("FFT Channel %i Ev %i",ipmt,ievent),nSamples/2,0,nSamples/2);
  cout << " FFT histogram " << hfft->GetName() << endl;

  std::vector<Double_t> realVec,imVec;
  for (int i = 0; i<nSamples; ++i) {
    double rl, im;
    fFFT->GetPointComplex(i,rl,im);
    std::complex<double>  c(rl,im); //.real or .imag accessors
    VectorComplex.push_back(c);
    // skip first bin which is pedestal
    //if( i> 0) hfft->SetBinContent(i+1,hfft->GetBinContent(i+1)+std::abs(c));
    if(i<nSamples/2) hfft->SetBinContent(i,hfft->GetBinContent(i)+std::abs(c));
    realVec.push_back( VectorComplex[i].real());
    imVec.push_back(VectorComplex[i].imag() );
  }
  /*
     realVec = MovingAverageFilter(realVec,5);
     imVec   = MovingAverageFilter(imVec,5);
     for(int i = 0; i < nSamples; i++){
     std::complex<double> c(realVec[i],imVec[i] );
     VectorComplex[i] = c;
     }
     */

  return VectorComplex;
}

std::vector< Double_t > anaCompass::inverseFFT(Int_t ipmt,Int_t ievent,std::vector<std::complex<double> > VectorComplex,std::vector<Double_t> sum)
{
  std::vector<Double_t > Signal;
  for(int is =0; is<nSamples; ++is) {
    if(ipmt==0) fInverseFFT->SetPoint(is, VectorComplex[is].real(),VectorComplex[is].imag() );
    else fInverseFFT->SetPoint( is, VectorComplex[is].real(),VectorComplex[is].imag() );
  }
  fInverseFFT->Transform();
  
  TH1D* hIfft = new TH1D(Form("inverseFFTPMT%i_Ev%i",ipmt,ievent),Form("Inverse FFT Channel %i Ev %i",ipmt,ievent),nSamples/2,0,nSamples/2);
  cout << " FFT histogram " << hIfft->GetName() << endl;

  Double_t norm = 0;
  for (int i = 0; i<nSamples; ++i) {
    double rl = fInverseFFT->GetPointReal(i);
    norm += rl;
    Signal.push_back(rl);
  }
  for(int i = 0; i < nSamples; i++){
    Signal[i] = Signal[i]*(sum[ipmt]/norm);
    hIfft->SetBinContent(i,hIfft->GetBinContent(i) + Signal[i]);
  }

  return Signal;
}

std::vector<Double_t> anaCompass::BubbleSort(std::vector<Double_t> A){
  int i, j, N = A.size();

  for (i = 0; i < N; i++){
    for (j = i + 1; j < N; j++){
      if (A[j] < A[i]){
        Double_t buff;
        buff = A[i];
        A[i] = A[j];
        A[j] = buff;
      }
    }
  }
  return A;
}

std::vector<Double_t> anaCompass::SimpleLowPassFilter(std::vector<Double_t> signal, Double_t alpha){
  std::vector<Double_t> FilteredSignal;
  //Double_t alpha = 0.992105;//100 MHz cut off frequency
  for(int i = 0; i < int(signal.size()); i++){
    if(i == 0) FilteredSignal.push_back(alpha*signal[i]);
    else{
      FilteredSignal.push_back(alpha*signal[i] + (1.-alpha)*FilteredSignal[i-1] );
    }
  }

  //cout<<"Filter Size "<<FilteredSignal.size()<<", Dsigi size "<<signal.size()<<endl;
  return FilteredSignal;
}

std::vector<Double_t> anaCompass::SimpleHighPassFilter(std::vector<Double_t> signal, Double_t alpha){
  std::vector<Double_t> FilteredSignal;
  //Double_t alpha = 0.992105;//100 MHz cut off frequency
  for(int i = 0; i < int(signal.size()); i++){
    if(i == 0) FilteredSignal.push_back(signal[i]);
    else{
      FilteredSignal.push_back( alpha*( signal[i]-signal[i-1] + FilteredSignal[i-1] )   );
    }
  }

  //cout<<"Filter Size "<<FilteredSignal.size()<<", Dsigi size "<<signal.size()<<endl;
  return FilteredSignal;

}

std::vector<Double_t> anaCompass::MovingAverageFilter(std::vector<Double_t> signal,Int_t N)
{

  std::vector<Double_t> filter;
  Int_t N2 = std::floor(N/2);
  for(int i = N2; i < int(signal.size())-N2; i++){
    Double_t sum = 0;
    for(int j = i-N2; j <= i+N2; j++){
      sum += signal[j];
    }
    sum /= N;
    filter.push_back(sum);
  }

  for(int i = 0; i < N2 ; i++){
    std::vector<Double_t>::iterator it;
    it = filter.begin();
    filter.insert(it,0.);
    filter.push_back(0);
  }
  return filter;
}

void anaCompass::plotWave(Int_t ientry, Int_t pmtNum, hitMap pmtHits ) {
  TString histName;
  TString peaksName; 
  TString sumName; 
  hitMapIter hitIter=pmtHits.begin();
  TPmtHit phit0 = hitIter->second;
  unsigned nhits = pmtHits.size();

  // baseline
  histName.Form("BaseEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hbase = (TH1F*) hBaseline[pmtNum]->Clone(histName);

   // baseline fit
  histName.Form("BaseFitEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hbaseFit = (TH1F*) hBaselineWMA[pmtNum]->Clone(histName);

  // title with first peak
  histName.Form("RawEv%i_PMT_%i",ientry,pmtNum);
  TH1F* rhist = (TH1F*) hPMTRaw[pmtNum]->Clone(histName);

  // title with first peak
  histName.Form("DerEv%i_PMT_%i",ientry,pmtNum);
  TH1F* dhist = (TH1F*) hPMTDerivative[pmtNum]->Clone(histName);
  

  // title with first peak
  histName.Form("WaveEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hist = (TH1F*) hPMTSignal[pmtNum]->Clone(histName);

  // fill peaks
  peaksName.Form("PeaksEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hpeaks = (TH1F*) hPMTSignal[pmtNum]->Clone(peaksName);
  hpeaks->Reset();
  for (hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
    TPmtHit phiti = hitIter->second;
    //printf(" adding to plot (%u,%u) \n", phiti.firstBin,phiti.lastBin);
    for(Int_t ibin=phiti.firstBin; ibin<=phiti.lastBin; ++ibin) hpeaks->SetBinContent(ibin, hist->GetBinContent(ibin));
  }
  
  if(isSimulation) {
    histName.Form("SimHitsEv%i_PMT_%i",ientry,pmtNum);
    TH1F* hist = (TH1F*) hPMTSim[pmtNum]->Clone(histName);
  }

  printf("plotWave: %s %.0f \n",hbase->GetName(),hbase->GetEntries());
}

// summed waves
void anaCompass::sumWave(Int_t ipmt) {
  for (Int_t ibin=0; ibin<hSum[ipmt]->GetNbinsX(); ++ibin) {
    hSum[ipmt]->SetBinContent(ibin,  hSum[ipmt]->GetBinContent(ibin) + hPMTSignal[ipmt]->GetBinContent(ibin)  );
  }
}
// background estimator (class TSpectrum).
std::vector<Double_t> anaCompass::getTBaseline(TH1D* hPMTRaw, TH1D* hBaseline, Double_t& ave, Double_t& aveSigma) 
{
  hBaseline->Reset();
  int nbins = hPMTRaw->GetNbinsX();
  // normalize to zero at start of waveform
  Int_t max = hPMTRaw->FindBin(0.4E-6);
  // fill source array
  for (int i = 0; i < nbins; i++) source[i]=hPMTRaw->GetBinContent(i + 1);
  spec->Background(source,nbins,200,TSpectrum::kBackIncreasingWindow,TSpectrum::kBackOrder2,kTRUE,TSpectrum::kBackSmoothing15,kFALSE);
  std::vector<Double_t> vbase;
  std::vector<Double_t> vsort;
  Double_t norm=0;
  for (int i = 0; i < nbins; i++) {
    vbase.push_back(source[i]);
    vsort.push_back(source[i]);
    //if(i<max) norm += source[i];
  }
  norm /= Double_t(max);
  for (int i = 0; i < nbins; i++) {
    vbase[i] -= norm;
    vsort[i] -= norm;
    hBaseline->SetBinContent(i + 1,vbase[i]*microSec);
  }
  std::sort(vsort.begin(),vsort.end());
  ave = vbase[0.5*vsort.size()];
  aveSigma  = vsort[0.16*vsort.size()];
  aveSigma = std::abs(aveSigma-ave);
  //printf(" bin %i norm is %f ave %f sigma %f  \n ",max,norm,ave,aveSigma);
  return vbase;
}
//
// background estimator removing pulses 
std::vector<Double_t> anaCompass::getBaseline(std::vector<Double_t> digi, hitMap pmtHits, TH1D* hBaselineFit , TH1D* hBaseline, 
    Double_t& ave, Double_t& aveSigma) 
{
  std::vector<Double_t> vbase=digi;
  for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
        TPmtHit phiti = hitIter->second;
        for(int iz = phiti.firstBin; iz <= phiti.lastBin; ++iz)  vbase[iz]=0;
  }
  
  hBaseline->Reset();
  hBaselineFit->Reset();
  // get estimate of noise, baseline.
  std::vector<Double_t> vwork;
  for(unsigned is=0; is< vbase.size() ; ++is) if(vbase[is]!=0) vwork.push_back(vbase[is]);
  getAverage(vwork,ave,aveSigma);

  // fill baseline
  for (int i = 0; i <  hBaseline->GetNbinsX(); i++) {
    hBaselineFit->SetBinContent(i + 1,vbase[i]*microSec);
  }

  //Access the fit resuts
  hBaselineFit->Fit("pol1","OQ");
  TF1 *fpoly=NULL;
  fpoly= hBaselineFit->GetFunction("pol1");
  if(!fpoly->IsValid()) {
    hBaselineFit->Fit("pol0","OQ");
    fpoly= hBaselineFit->GetFunction("pol0");
  }
  fpoly->SetLineColor(kGreen);
  ave += fpoly->GetParameter(0); //value of 1st parameter
  // fill baseline function from fit.
  for (int i = 0; i <  hBaselineFit->GetNbinsX(); i++) {
    Double_t xbin = hBaselineFit->GetXaxis()->GetBinCenter(i+1);
    Double_t vbin = fpoly->Eval(xbin);  
    vbase[i]=vbin;
  }

  // fill histogram
  for (int i = 0; i < hBaseline->GetNbinsX(); i++) hBaseline->SetBinContent(i + 1,vbase[i]);
 
  // return vector
  return vbase;
}


void anaCompass::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma) 
{
  // get first estimate of noise, baseline.
  std::vector<Double_t> vsort;
  for(unsigned is=0; is< digi.size() ; ++is) vsort.push_back(digi[is]);
  std::sort(vsort.begin(),vsort.end());
  ave  = vsort[0.5*vsort.size()];
  sigma  = vsort[0.16*vsort.size()];
  sigma = std::abs(sigma-ave);
}
//implimented as in zugec et al, arXiv:1601.04512v1
std::vector<Double_t> anaCompass::differentiate(std::vector<Double_t> v, unsigned nstep)
{
  std::vector<Double_t> d;
  unsigned nsamples = v.size();
  Double_t sump=0;
  Double_t summ=0;
  d.push_back(0); // first entry is zero
  for(unsigned i=1; i<nsamples; ++i) {
    unsigned i2 = 2*i;
    unsigned max = TMath::Min( nstep, i);
    max = TMath::Min(max, nsamples - 1 -i);
    // beginning, middle, end cases
    if(i<=nstep && i2 <= nsamples -1 ) {
      sump = sump - v[i] + v[i2-1] + v[i2];
      summ = summ + v[i-1];
    }
    else if(i>nstep && i+nstep <= nsamples -1){
      sump = sump - v[i] + v[i+nstep];
      summ = summ + v[i-1]-v[i-1-nstep];
    }
    else if(i+nstep> nsamples-1 && i2 > nsamples-1) {
      sump = sump - v[i];
      summ = summ + v[i-1] - v[i2-nsamples-1] - v[i2-nsamples];
    }
    d.push_back(sump-summ);
  }
  return d;
}
peakType anaCompass::derivativePeaks(std::vector<Double_t> v,  Int_t nsum, Double_t rms,std::vector<Int_t>& peakKind ) 
{
  peakType peakList;
  peakKind.clear();
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  std::vector<double> crossingTime;

  unsigned vsize = v.size();
  Double_t cut = derivativeSigma*rms;
  Double_t ncut = -derivativeSigma*rms;
  // find all crossings
  for( unsigned ibin=1; ibin< vsize; ++ibin ) {
    Double_t u = double(ibin)*timeUnit*microSec;
    if( v[ibin]>cut && v[ibin-1] <cut ) {
      crossings.push_back(UPCROSS);
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( v[ibin]<cut && v[ibin-1] > cut ) {
      crossings.push_back(UPCROSS);
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( v[ibin]<ncut && v[ibin-1] > ncut ) {
      crossings.push_back(DOWNCROSS);
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    } else if( v[ibin]>ncut && v[ibin-1] < ncut ) {
      crossings.push_back(DOWNCROSS);
      crossingBin.push_back(ibin);
      crossingTime.push_back(u);
    }
  }

  if(crossings.size()<4) return peakList;

  // label found crossings, intially all false
  std::vector<bool> crossingFound;
  crossingFound.resize(crossings.size());
  for(unsigned jc=0; jc<crossings.size(); ++jc)  {
    //printf(" crossing %i bin %i time %f type %i \n",jc,crossingBin[jc],crossingTime[jc],crossings[jc]);
    crossingFound[jc]=false;
  }


  // parse crossings to make pairs 
  unsigned ip =0; 
  //printf("crossings %zu %E %E \n",crossings.size(),cut,ncut);
  while ( ip<= crossings.size() -4 ) {
    // UP UP DOWN DOWN 
    if(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS&&crossings[ip+2]==DOWNCROSS&&crossings[ip+3]==DOWNCROSS) {
      // check zero crossings
      unsigned nzero=0;
      for(unsigned ibin = crossingBin[ip+1]; ibin<crossingBin[ip+2]; ++ibin) if(v[ibin]>0&&v[ibin+1]<0) ++nzero;
      if(nzero<=3) {
        //printf(" peak %i time %f (%i %i %i %i ) nzero %u \n",ip,crossingTime[ip],crossings[ip],crossings[ip+1],crossings[ip+2],crossings[ip+3],nzero);
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+3]) );
        peakKind.push_back(0);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+3]-crossingBin[ip]),double(0));//sigma:d0:step:dstep
        crossingFound[ip]=true;
        crossingFound[ip+1]=true;
        crossingFound[ip+2]=true;
        crossingFound[ip+3]=true;
        ip=ip+4;
      } else 
        ++ip;
    } else {
      ++ip;
    }
  }
  //printf("peaks with fours %zu  \n",peakList.size());

  // pick up UP UP, DOWN DOWN cases
  ip =0; 
  while ( ip<= crossings.size() -2 ) {  
    if(!crossingFound[ip]&&!crossingFound[ip+1]) {
      // UP UP
      if(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS) {
        //printf(" up up peak %i time %f (%i %i )\n",ip,crossingTime[ip],crossings[ip],crossings[ip+1]); 
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
        peakKind.push_back(1);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(1));//sigma:d0:step:dstep
        // DOWN DOWN 
      } else if(crossings[ip]==DOWNCROSS&&crossings[ip+1]==DOWNCROSS) {
        //printf(" down down peak %i time %f (%i %i )\n",ip,crossingTime[ip],crossings[ip],crossings[ip+1]); 
        peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
        peakKind.push_back(2);
        ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(2));//sigma:d0:step:dstep
      }
    }
    ip=ip+2;
  }
  //printf("peaks with twos %zu  \n",peakList.size());
 
  // extend pulses to zero derivative
  /*
  for(unsigned ip=0; ip<peakList.size(); ++ip)  {
    // high direction
    unsigned high = std::get<1>(peakList[ip]);
    unsigned next = vsize;
    if(ip<peakList.size()-1) next = std::get<1>(peakList[ip+1]);
    for(unsigned kp= high; kp < next ; ++kp ) {
      std::get<1>(peakList[ip])=kp;
      if( v[kp] <0 && v[kp+1]>0 ) break;
    }

    // low direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned prev = 1;
    if(ip>0) prev = TMath::Max(prev,std::get<0>(peakList[ip-1]));
    for(unsigned kp= low; kp > prev ; --kp ) {
      std::get<0>(peakList[ip])=kp;
      if( v[kp] >0 && v[kp-1]<0 ) break;
    }
  }
  */

  // return list
  return peakList;
}
void anaCompass::trimPeaks(peakType& peakList, std::vector<Double_t> sv)
{

  unsigned svsize=sv.size();
  for(unsigned ip=0; ip<peakList.size(); ++ip)  {
    // high direction
    unsigned low = std::get<0>(peakList[ip]);
    unsigned high = std::get<1>(peakList[ip]);
    Double_t vmax=-9999;
    unsigned kmax=0;
    for(unsigned kp=low; kp<=high; ++kp) {
      if(sv[kp]>vmax) {
        vmax=sv[kp];
        kmax=kp;
      }
    }
    for(unsigned kp= kmax; kp <= high ; ++kp ) {
      std::get<1>(peakList[ip])=kp;
      if(kp==svsize) break;
      if(sv[kp+1]<0) break;
    }
    // low direction
    for(unsigned kp= kmax; kp >= low ; --kp ) {
      std::get<0>(peakList[ip])=kp;
      if(kp==0) break;
      if( sv[kp-1]<0) break;
    }
  }
}




//Weighted Moving Average baseline from Zugec et al 
std::vector<Double_t> anaCompass::getBaselineWMARecursive(Double_t ave, std::vector<Double_t> sig, std::vector<Double_t> weight,Int_t NWindow)
{
  // default baseline is zero
  std::vector<Double_t> baseline(sig.size(),ave);
  
   // starting values 
  Double_t KW=0;
  Double_t CW=0;
  Double_t SW=0;
  Double_t KS=0;
  Double_t CS=0;
  Double_t SS=0;

  // recursive sums
  Double_t kcos = TMath::Cos( TMath::Pi()/Double_t(NWindow));
  Double_t ksin = TMath::Sin( TMath::Pi()/Double_t(NWindow));
  //printf(" kcos %f ksin %f \n",kcos,ksin);
  Int_t resetBin = 100;

  for(int iw=0; iw<int(sig.size()); ++iw) {
    // reset
    if(iw%resetBin ==  0) {
      int low =  TMath::Max(0,iw-NWindow);
      int high = TMath::Min(iw+NWindow,int(sig.size())-1);
      KW=0;CW=0;SW=0;KS=0; CS=0;SS=0;
      for(int jp =low; jp < high; ++jp) {
        Double_t cosj = TMath::Cos(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        Double_t sinj = TMath::Sin(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
        KW+= weight[jp];
        KS+= sig[jp]*weight[jp];
        CW+= weight[jp]*cosj;
        CS+= sig[jp]*weight[jp]*cosj;
        SW+= weight[jp]*sinj;
        SS+= sig[jp]*weight[jp]*sinj;
      }
    } else {
      //save previous iteraton
      int i1 = iw+NWindow;
      int i2 = iw-1-NWindow;
      Double_t CWlast=CW;
      Double_t CSlast=CS;

      // case sums 
      CW = kcos*CWlast+ksin*SW;
      CS = kcos*CSlast+ksin*SS;
      SW = kcos*SW-ksin*CWlast;
      SS = kcos*SS-ksin*CSlast;

      if(iw  <= NWindow && iw + NWindow <= sig.size() - 1){
        KW += weight[i1];
        KS += sig[i1]*weight[i1];
        CW += -1.*weight[i1];
        CS += -1.*sig[i1]*weight[i1];
      } else if(iw > NWindow && iw + NWindow <= sig.size() - 1){
        KW += -1.*weight[i2] + weight[i1];
        KS += -1.* sig[i2]*weight[i2] +  sig[i1]*weight[i1];
        CW +=  kcos*weight[i2] - weight[i1];
        CS +=  kcos*sig[i2]*weight[i2] - sig[i1]*weight[i1];
        SW += -1.*ksin*weight[i2] ;
        SS += -1.*ksin*sig[i2]*weight[i2] ;
      } else if(iw  > NWindow && iw + NWindow > sig.size() - 1){
        KW += -1.*weight[i2];
        KS +=  -1.*sig[i2]*weight[i2];
        CW +=  kcos*weight[i2];
        CS +=  kcos*sig[i2]*weight[i2];
        SW += -1.*ksin*weight[i2];
        SS += -1.*ksin*sig[i2]*weight[i2];
      }  // sum baseline 
    }
    if( TMath::Abs(KW+CW)>1E-12) baseline[iw]=(KS+CS)/(KW+CW);
  }
  return baseline;
}

//Weights Weighted Moving Average baseline from Zugec et al 
std::vector<Double_t> anaCompass::getBaselineWeights(unsigned arraySize, peakType peakList,Int_t& maxwidth)
{
  std::vector<Double_t> weight(arraySize,0);
  if(peakList.size() == 0) return weight;

  // print peaks
  // for(unsigned ip=0; ip<peakList.size(); ++ip) printf(" ip %u alpha %u  beta  %u \n",ip,std::get<0>(peakList[ip]) , std::get<0>(peakList[ip]) );
  
  // construct the weights
  int alpha;
  int beta;
  maxwidth=0;
  for(int ip=0; ip<=int(peakList.size()); ++ip) {
    if(ip==0) beta = -1;
    else beta = std::get<1>(peakList[ip-1]);
    if(ip==peakList.size()) alpha = weight.size();
    else alpha = std::get<0>(peakList[ip]);
    int width = std::get<1>(peakList[ip])-std::get<0>(peakList[ip]);
    if(width>maxwidth) maxwidth=width;
    for(int jp=beta+1; jp<alpha; ++jp)  weight[jp]=alpha-beta-1;
    //if(jp==alpha-1||jp==beta+1) printf(" ip %u jp %u beta %i alpha %i  w %.1f P %lu \n",ip,jp,beta,alpha, weight[jp], sig.size());
  }
  return weight;
}

//Weighted Moving Average baseline from Zugec et al 
std::vector<Double_t> anaCompass::getBaselineWMA(Double_t ave, std::vector<Double_t> sig, std::vector<Double_t> weight,Int_t NWindow)
{
  // default baseline is zero
  std::vector<Double_t> basenon(sig.size(),ave);

  // non recursive sums 
  for(int iw=0; iw<int(sig.size()); ++iw) {
    // non-recursive 
    Double_t KWnon=0;
    Double_t CWnon=0;
    Double_t KSnon=0;
    Double_t CSnon=0;
    int jlow =  TMath::Max(0,iw-NWindow);
    int jhigh = TMath::Min(iw+NWindow,int(sig.size())-1);
    for(int jp=jlow; jp<=jhigh; ++jp) {
      Double_t cosj = TMath::Cos(  Double_t(jp-iw)*TMath::Pi()/Double_t(NWindow) );
      KWnon+= weight[jp];
      KSnon+= sig[jp]*weight[jp];
      CWnon+= weight[jp]*cosj;
      CSnon+= sig[jp]*weight[jp]*cosj;
    }
    basenon[iw]= (KSnon+CSnon)/(KWnon+CWnon);
  }
  // end of non recursive baseline
  return basenon;
}

//Weighted Moving Average
std::vector<Double_t> anaCompass::getBaselineWMANeil(Double_t ave,std::vector<Double_t> sig, std::vector<Double_t> weight,Int_t N)
{
  // default baseline is ave
  std::vector<Double_t> Baseline(sig.size(),ave);
  Double_t pi = TMath::Pi();
  /*
  //obnoxious way of making things faster
  //i.e. run times does not depend on how big the window, N, is. 
  //It only depends on total number of points
  */
  Double_t K_sw = 0;
  Double_t C_sw = 0;
  Double_t S_sw = 0;
  Double_t K_w = 0;
  Double_t C_w = 0;
  Double_t S_w = 0;
  Double_t kc = cos(pi/N);
  Double_t ks = sin(pi/N);
  Int_t resetBin = 100;
  if(sig.size() == 500) resetBin = 20;
  for(int i = 0;i< sig.size();i++){
    Int_t hiWindow = std::min(i+N,(Int_t)sig.size() -1);
    Int_t loWindow = std::max(0,i-N);
    if(i%resetBin ==  0){
      K_sw = 0;C_sw = 0;S_sw = 0;K_w = 0;C_w = 0;S_w = 0;
      for(int j = loWindow;j<hiWindow;j++){
        K_sw += sig[j]*weight[j];
        C_sw += sig[j]*weight[j]*cos((j-i)*pi/N);
        S_sw += sig[j]*weight[j]*sin((j-i)*pi/N);
        K_w += weight[j];
        C_w += weight[j]*cos((j-i)*pi/N);
        S_w += weight[j]*sin((j-i)*pi/N);
      }
    }
    else{
      if(i - N <= 0 && i + N <= sig.size() - 1){
        K_sw = K_sw + sig[i+N]*weight[i+N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw - sig[i+N]*weight[i+N];
        S_sw = kc*S_sw - ks*tempC_sw ;

        K_w = K_w + weight[i+N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w - weight[i+N];
        S_w = kc*S_w - ks*tempC_w ;
      }
      else if(i - N > 0 && i + N <= sig.size() - 1){
        K_sw = K_sw - sig[i-1-N]*weight[i-1-N] + sig[i+N]*weight[i+N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw + kc*sig[i-1-N]*weight[i-1-N] - sig[i+N]*weight[i+N];
        S_sw = kc*S_sw - ks*tempC_sw - ks*sig[i-1-N]*weight[i-1-N];

        K_w = K_w - weight[i-1-N] + weight[i+N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w + kc*weight[i-1-N] - weight[i+N];
        S_w = kc*S_w - ks*tempC_w - ks*weight[i-1-N];
      }
      else if(i - N > 0 && i + N > sig.size() - 1){
        K_sw = K_sw - sig[i-1-N]*weight[i-1-N];
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw + kc*sig[i-1-N]*weight[i-1-N];
        S_sw = kc*S_sw - ks*tempC_sw - ks*sig[i-1-N]*weight[i-1-N];

        K_w = K_w - weight[i-1-N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w + kc*weight[i-1-N];
        S_w = kc*S_w - ks*tempC_w-ks*weight[i-1-N];
      }
      else if(i - N <=0 && i + N > sig.size() - 1){
        K_sw = K_sw;
        Double_t tempC_sw = C_sw;
        C_sw = kc*C_sw + ks*S_sw;
        S_sw = kc*S_sw - ks*tempC_sw;

        K_w = K_w - weight[i-1-N];
        Double_t tempC_w = C_w;
        C_w = kc*C_w + ks*S_w;
        S_w = kc*S_w - ks*tempC_w;
      }
      else{ cout<<"you messed up the logic on the moving baseline"<<endl;}
    }
    //cout<<"ientry"<<ientry<<" "<<i<<"...K_sw "<<K_sw<<", C_sw "<<C_sw<<", S_sw "<<S_sw<<", K_w "<<K_w<<", C_w "<<C_w<<"..."<< (K_sw+C_sw)/(K_w+C_w)<<endl;
    Baseline[i]=(K_sw+C_sw)/(K_w+C_w);
  }
  return Baseline;
}


void anaCompass::clear()
{
  Channel=0;
  Timestamp=0;
  Board=0;
  Energy=0;
  EnergyShort=0;
  Flags=0;
  Samples->Reset();
}

void anaCompass::Init()
{
  Samples =0;
  fCurrent = -1;
  fChain->SetBranchAddress("Channel", &Channel, &b_Channel);
  fChain->SetBranchAddress("Timestamp", &Timestamp, &b_Timestamp);
  fChain->SetBranchAddress("Board", &Board, &b_Board);
  fChain->SetBranchAddress("Energy", &Energy, &b_Energy);
  fChain->SetBranchAddress("EnergyShort", &EnergyShort, &b_EnergyShort);
  fChain->SetBranchAddress("Flags", &Flags, &b_Flags);
  fChain->SetBranchAddress("Samples", &Samples, &b_Samples);
}

Long64_t anaCompass::getEvent(Long64_t event) 
{
  if(fChain->LoadTree(event)<0) return -1;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  Long64_t nb = fChain->GetEvent(event);
  return nb;
}
