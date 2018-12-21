/* MG revised */
#include "anaRun.hh"

anaRun::anaRun(TString tag, Int_t maxEvents)
{
  printf(" starting anaRun tag %s \n",tag.Data());
  int printInterval=10;
  firstChargeCut=0.1;
  ran = new TRandom3();
  Double_t simHitMatchTime=0;
  Double_t sigma=0;
  fsigma=sigma;
  if(sigma==0) fsigma=5;
  windowSize=15;
  nSigma=5;
  aveWidth=20;
  spec = new TSpectrum();
  //Int_t irunStop = irunStart;
  TString outFileName ; outFileName.Form("%s_Ev_%i_derivative.root",tag.Data(),maxEvents);
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  ntCal =  new TNtuple("ntCal","ntuple Cal","iev:ipmt:base:sigma:dbase:dsigma");
  ntHit =  new TNtuple("ntHit","ntuple Hit","npmt:nhits:order:istart:time:tstart:q:nwidth:qmax");
  ntNHit = new TNtuple("ntNHit","negative ntuple Hit","npmt:nhits:order:istart:time:tstart:q:nwidth:qmax");
  ntDer =  new TNtuple("ntDer"," deriviative ","sigma:d0:kover:type");
  ntEvent= new TNtuple("ntEvent","ntuple Event","entry:n0:n1:t00:t01:t10:t11:qp0:qp1:q00:q01:q10:q11:qsum0:qsum1");
  ntPulse= new TNtuple("ntPulse"," pulse ","sum:shigh:slow:nsamp:kover:qlow:qhigh:klow:khigh");
  ntWave = new TNtuple("ntWave"," wave ","event:v:d");

  hPeakNWidth = new TH1I("PeakNWidth","PeakNWidth",100,0,100);
  hAllNWidth = new TH1I("AllNWidth","AllNWidth",100,0,100);
  hNWidthCut = new TH1I("NWidthCut","NWidth cut ",100,0,100);
  hDerWidth = new TH1I("DerWidth","sliding window pulse number Width",2*maxHalfLength,0,2*maxHalfLength);
  hDerAfter = new TH1D("DerAfter"," after cut derivative",1000,-.004,.004);
  hSlideHigh = new TH1D("SlideHigh"," sliding window high threshold ",1000,0,1);
  hSlideLow  = new TH1D("SlideLow"," sliding window low threshold",1000,0,.1);
  hSlideQSum = new TH1D("SlideQSum"," sliding window Q sum",100,0,100);
  hWeight = new TH1D("Weight","baseline weight",2*baseLineHalfWindow,0,2*baseLineHalfWindow);
  hWeightOne = new TH1D("WeightOne","single baseline weight",2*baseLineHalfWindow,0,2*baseLineHalfWindow);

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


  //for(Int_t irun = irunStart; irun<irunStop +1; irun++){
  // open ouput file and make some histograms
  TString fileName; fileName.Form("rootData/%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  TFile *fin = new TFile(fileName,"readonly");
  if(fin->IsZombie()) {
    printf(" couldnt open file %s\n",fileName.Data());
    return;
  }
  else 
    printf("  found file %s \n",fileName.Data() );

  // get pmtTree from file 
  fin->GetObject("pmtTree",pmtTree);
  Long64_t nentries = pmtTree->GetEntries();
  printf(" number of entries is %lld derivative smoothing = %i \n",nentries,windowSize);

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
  pmtSimulation=NULL;
  isSimulation=false;
  if(pmtTree->FindBranch("pmtSimulation")) {
    pmtSimulation = new TPmtSimulation();
    if(pmtTree->SetBranchAddress("pmtSimulation", &pmtSimulation)==TTree::kMatch) isSimulation=true;
  }
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
  int nHists=0;
  int nHistsMax=300;
  Float_t qpmt[NPMT][MAXHIT];
  Float_t tpmt[NPMT][MAXHIT];
  Int_t   npmtHit[NPMT];
  Double_t qsum[NPMT];
  Double_t qped[NPMT];
  Int_t nped[NPMT];
  Int_t totalHits[NPMT]={0,0};

  // loop over entries
  if(maxEvents>0) nentries=maxEvents;
  for (UInt_t ientry=0; ientry<nentries; ientry++) {
    pmtTree->GetEntry(ientry);
    if(pmtEvent->time.size() == 0) continue;
    nSamples = pmtEvent->time.size();
    for(int ipmt=0; ipmt<NPMT ; ++ ipmt) {
      npmtHit[ipmt]=0;
      qped[ipmt]=0;
      qsum[ipmt]=0;
      nped[ipmt]=0;
      for(int ihit=0; ihit<MAXHIT ; ++ihit) {
        qpmt[ipmt][ihit]=0;
        tpmt[ipmt][ihit]=0;
      }
    }

    // check how many pmts we have
    int gotPMT = 0;
    if(pmtEvent->volt1.size()>0)  ++gotPMT;
    if(pmtEvent->volt2.size()>0)  ++gotPMT;
    if(ientry==0) printf(" .... events %lld samples %i PMT0 %zu PMT1 %zu \n",pmtTree->GetEntries(),nSamples,pmtEvent->volt1.size(),pmtEvent->volt2.size());
    
    if(gotPMT<1) return;

    TString name; name.Form("%s_Ev%i",tag.Data(),ientry);
    // define pmt signal histograms
    if(ientry==0) {
      source = new Double_t[nSamples];
      Double_t maxLife = pmtEvent->time[nSamples-1]*1E6;
      timeUnit=pmtEvent->time[1]-pmtEvent->time[0];
      simHitMatchTime=100.0*timeUnit;
      printf(" \n\n ***** setting time unit %E maxLife %f \n",timeUnit,maxLife);
      if(isSimulation) printf(" \n\n ***** setting matching time  %E \n",simHitMatchTime);
      if(isSimulation) ntSimMatch = new TNtuple("ntSimMatch"," hit matches ","pmt:nsim:nhit:match:nnot:nmiss");
      for(int ipmt=0; ipmt<gotPMT; ++ipmt) {
        hPMTRaw[ipmt] = new TH1D(Form("PMTRaw%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hPMTSignal[ipmt] = new TH1D(Form("PMTSignal%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hPMTDerivative[ipmt] = new TH1D(Form("PMTDeriv%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hSum[ipmt] = new TH1D(Form("SumPmt%i",ipmt),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hSum[ipmt]->GetXaxis()->SetTitle(" seconds ");
        hBaseline[ipmt]  = new TH1D(Form("Baseline%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hBaselineFit[ipmt]  = new TH1D(Form("BaselineFit%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hNHits[ipmt] = new TH1D(Form("NHits%i",ipmt),Form(" number of hits PMT %i ",ipmt),10,0,50);
        hNegNHits[ipmt] = new TH1D(Form("NegNHits%i",ipmt),Form(" number of neg hits PMT %i ",ipmt),10,0,50);
        hLife[ipmt] = new TH1D(Form("Life%i",ipmt),Form(" lifetime PMT %i ",ipmt),2000,-maxLife,maxLife);
        hLife[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
        hNLife[ipmt] = new TH1D(Form("NLife%i",ipmt),Form(" negative pulse lifetime PMT %i ",ipmt),2000,-maxLife,maxLife);
        hNLife[ipmt]->GetXaxis()->SetTitle(" micro-seconds ");
         //hLife[ipmt]->Sumw2();
        if(isSimulation) {
          hPMTSim[ipmt] = new TH1D(Form("PMTSim%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
          hPMTSimHitMatch[ipmt] = new TH1D(Form("PMTSimHitMatch%i_%s",ipmt,tag.Data()),"",1000,0,100*simHitMatchTime);
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
    
    std::vector<std::vector<Double_t> > baselineDigi; baselineDigi.resize(NPMT);

    // loop over PMT 1
    for(unsigned isample = 0; isample < pmtEvent->volt1.size(); isample++){
      Double_t volt0 = pmtEvent->volt1[isample];
      Double_t digi0 = -1.0*(double(volt0));
      ddigi[0].push_back(digi0);
      ndigi[0].push_back(-digi0);
      hPMTRaw[0]->SetBinContent(isample+1,digi0);
    }

    // loop over PMT 2
    for(unsigned isample = 0; isample < pmtEvent->volt2.size(); isample++){
      Double_t volt1 = pmtEvent->volt2[isample];
      Double_t digi1 = -1.0*(double(volt1));
      ddigi[1].push_back(digi1);
      ndigi[1].push_back(-digi1);
      hPMTRaw[1]->SetBinContent(isample+1,digi1);
    }

     
    
    //FFT
    if(ientry<2) {
       for(int ipmt=0; ipmt<2; ++ipmt) if(ddigi[ipmt].size()>0) FFT(ipmt,ientry,ddigi[ipmt]);
    }
  
    // derivative PMT
    Double_t derAve[NPMT], derSigma[NPMT];
    for(int j = 0; j < gotPMT; j++){
      deriv[j] = differentiate( ddigi[j],windowSize);
      nderiv[j] = differentiate( ndigi[j],windowSize);
      for(unsigned isample = 0; isample < deriv[j].size(); isample++){
        hPMTDerivative[j]->SetBinContent(isample+1,deriv[j][isample]);
        ntWave->Fill(ientry,ddigi[j][isample],deriv[j][isample]);
        getAverage(deriv[j],derAve[j],derSigma[j]);
      }
    }

      
    
    for(int pmtNum = 0 ; pmtNum < gotPMT; pmtNum++){
  
      /* 
      ** puluse finding and hit making
      */
      Double_t minDev = 0*sDev[pmtNum];
      Double_t maxDev = fsigma*sDev[pmtNum];
      Double_t firstTime, firstCharge;
      peakType peakList = derivativePeaks(deriv[pmtNum],windowSize,derSigma[pmtNum]);
      hitMap  pmtHits = makeHits(peakList, ddigi[pmtNum],maxDev,firstTime,firstCharge);

      // negative pulses
      Double_t nfirstTime, nfirstCharge;
      peakType npeakList = derivativePeaks(nderiv[pmtNum],windowSize,derSigma[pmtNum]);
      hitMap  npmtHits = makeHits(npeakList, ndigi[pmtNum],maxDev,nfirstTime,nfirstCharge);

     // get baseline
      baselineDigi[pmtNum] = getBaseline(ddigi[pmtNum],pmtHits,hBaselineFit[pmtNum],hBaseline[pmtNum],baseline[pmtNum],sDev[pmtNum]); 
      ntCal->Fill(ientry,pmtNum,baseline[pmtNum],sDev[pmtNum],derAve[pmtNum],derSigma[pmtNum]);
      if(ientry%printInterval==0) printf(" ...... %i ave %.4E sdev0 %.4E dave %.4E dsigma %.4E \n",
          ientry,baseline[pmtNum],sDev[pmtNum],derAve[pmtNum],derSigma[pmtNum]);

      // subtract baseline
      maxSample[pmtNum]=0;
      hPMTSignal[pmtNum]->Reset();
      if(isSimulation) hPMTSim[pmtNum]->Reset();
      for(int i = 0; i < int(ddigi[pmtNum].size()); i++){
        if(ddigi[pmtNum][i]>maxSample[pmtNum])  maxSample[pmtNum] = ddigi[pmtNum][i];
        ///ddigi[pmtNum][i] = ddigi[pmtNum][i]-baselineDigi[pmtNum][i];// + sDev[pmtNum];
        //ndigi[pmtNum][i] = ndigi[pmtNum][i]+baselineDigi[pmtNum][i];// + sDev[pmtNum];
        hPMTSignal[pmtNum]->SetBinContent(i,ddigi[pmtNum][i]-baselineDigi[pmtNum][i]);
      }


      unsigned nhits = pmtHits.size();
      hNHits[pmtNum]->Fill(nhits);
      unsigned negnhits = npmtHits.size();
      hNegNHits[pmtNum]->Fill(negnhits);

 
      if(ientry%printInterval==0) printf(" pmt  %i peaktime %lu nhits %u \n ",pmtNum,peakList.size(),nhits);

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


      hQFirst->Fill(firstTime,firstCharge);
      //if(firstCharge<firstChargeCut&&phit0.qsum>firstChargeCut) 
      if(firstTime==1E9&&phit0.qsum>firstChargeCut) 
        printf(" WARNING NO FIRST PULSE event %i  pmt %i pulses %i qhit %f time %f first %f charge %f \n",
            ientry,pmtNum,nhits,phit0.qsum,phit0.startTime*1E6,firstTime,firstCharge);

      int hitCount=0;
      for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
        TPmtHit phiti = hitIter->second;
        Double_t phitTime =  phiti.startTime*1E6-firstTime; 
        hQStart->Fill(phitTime,phiti.qsum);
        Int_t nwidth = phiti.lastBin - phiti.firstBin +1;
        Int_t istartBin =  hLife[pmtNum]->FindBin(phitTime); 
        qsum[pmtNum] += phiti.qsum;
        ntHit->Fill(pmtNum,nhits,hitCount++,istartBin, phiti.startTime*1E6, phitTime,phiti.qsum,nwidth,phiti.qpeak);
        if(phitTime!=0&&firstCharge>firstChargeCut) hLife[pmtNum]->SetBinContent( istartBin, hLife[pmtNum]->GetBinContent(istartBin)+phiti.qsum);
      }

      // compare hits and simulation 
      if(isSimulation) {
        std::vector<Double_t> startTime=pmtSimulation->startTime;
        // vector to hold best matches
        std::vector<Double_t> hitMatch(pmtHits.size(),100.);
        std::vector<Int_t> hitMatchNumber(pmtHits.size(),-1);
        for(unsigned isim = 0 ; isim < startTime.size(); ++ isim ) {
          hPMTSim[pmtNum]->Fill(startTime[isim]);
          unsigned ihit=0;
          for (hitMapIter hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
            TPmtHit phiti = hitIter->second;
            if(abs(phiti.peakt*timeUnit-startTime[isim])<hitMatch[ihit]) {
              hitMatch[ihit] = abs(phiti.peakt*timeUnit-startTime[isim]);
              hitMatchNumber[ihit]=isim;
            }
            ++ihit;
          }
        }
        unsigned nmatch=0;
        unsigned nnot=0;
        for(unsigned ihit=0; ihit<hitMatch.size(); ++ihit) {
          if(hitMatch[ihit]<simHitMatchTime) ++nmatch;
          else ++nnot;
          hPMTSimHitMatch[pmtNum]->Fill(hitMatch[ihit]);
          //printf(" %u (%E) %E number %i nmatch %i nnot %i \n",ihit, simHitMatchTime,hitMatch[ihit],hitMatchNumber[ihit], nmatch,nnot );
        }

        // count missed
        unsigned nmiss=0;
        for(unsigned isim = 0 ; isim < startTime.size(); ++ isim ) {
          bool matched=false;
          for(unsigned ihit=0; ihit< hitMatchNumber.size(); ++ihit) if(int(isim)==hitMatchNumber[ihit]) matched=true;
          if(!matched) ++nmiss;
        }

        ntSimMatch->Fill(float(pmtNum),float(startTime.size()),float(pmtHits.size()),float(nmatch),float(nnot),float(nmiss));
        if(pmtNum==0) simMatchStats->fill(startTime.size(),pmtHits.size(),nmatch,nnot,nmiss);
        //printf(" %i PMT%i ngen %zu  nhits %lu nmatches %u  not %u \n",ientry,pmtNum,startTime.size(),pmtHits.size(),nmatch,nnot);
        //if(ientry%printInterval==0) simMatchStats->print();
        simMatchStats->print();
      }
      if(nHists<nHistsMax) {
        plotWave(ientry,pmtNum,pmtHits );
        ++nHists;
      }

      // summed wave forms only if passes first charge cut
      if(firstCharge>firstChargeCut) sumWave(pmtNum);

      // negative pulse ntuple
      // hitMapIter hitIter;
      hitMapIter nhitIter;
      nhitIter=npmtHits.begin();
      TPmtHit nphit0 = nhitIter->second;

      hitCount=0;
      for (hitMapIter nhitIter=npmtHits.begin(); nhitIter!=npmtHits.end(); ++nhitIter) {
        TPmtHit phiti = nhitIter->second;
        Int_t nwidth = phiti.lastBin - phiti.firstBin +1;
        Double_t phitTime =  phiti.startTime*1E6-nfirstTime;
        hNegQStart->Fill(phitTime,phiti.qsum);
        Int_t istartBin =  hNLife[pmtNum]->FindBin(phitTime);
        ntNHit->Fill(pmtNum,nhits,hitCount++,istartBin, phiti.startTime*1E6,phitTime,phiti.qsum,nwidth,phiti.qpeak);
        if(phitTime!=0&&nfirstCharge>firstChargeCut) hNLife[pmtNum]->SetBinContent( istartBin, hNLife[pmtNum]->GetBinContent(istartBin)+phiti.qsum);
      }
      if(nfirstTime==1E9&&nphit0.qsum>firstChargeCut) 
        printf(" WARNING NO NEGATIVE FIRST PULSE event %i  pmt %i pulses %i qhit %f time %f first %f charge %f \n",
            ientry,pmtNum,negnhits,nphit0.qsum,nphit0.startTime*1E6,nfirstTime,nfirstCharge);



      // save for npmtHit ntuple
      qpmt[pmtNum][0]=phit0.qsum;
      tpmt[pmtNum][0]=phit0.startTime*1E6;    // milli-sceconds 
      qpmt[pmtNum][1]=phit1.qsum;
      tpmt[pmtNum][1]=phit1.startTime*1E6;

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

      // add noise hit
      // Int_t nwidth0 = phit0.lastBin - phit0.firstBin +1;
      // ntHit->Fill(pmtNum,nhits,hitCount++,0,qped[pmtNum],0,nwidth0,qped[pmtNum]/Double_t(nwidth0));
      //ntHit->Fill(pmtNum,nhits,hitCount++,0,qped[pmtNum],0,nwidth0,sDev[pmtNum]*Double_t(aveWidth));

    }

    //cout << ddigi[0].size() << "   " <<  nped[0]  << "  " << nped[1] << " " << qped[0] << "  " << qped[1] <<  endl;

    if(npmtHit[0]>0||npmtHit[1]>0) 
      ntEvent->Fill(ientry,npmtHit[0],npmtHit[1],tpmt[0][0],tpmt[0][1],tpmt[1][0],
          tpmt[1][1],qped[0],qped[1],qpmt[0][0],qpmt[0][1],qpmt[1][0],qpmt[1][1],qsum[0],qsum[1]);
  }
  printf(" total hits %i, %i out of %lld events\n",totalHits[0],totalHits[1],nentries);

  outfile->Purge();
  outfile->Write();

  return;
}


hitMap anaRun::makeHits(peakType peakList, std::vector<Double_t> ddigi,Double_t sigma, Double_t& firstTime, Double_t& firstCharge) 
{

  firstTime=1E9;
  firstCharge=0;
  hitMap pmtHits;
  if(peakList.size()<1) return pmtHits;
  Double_t qmax=0;
  
  for(unsigned ip=0; ip<peakList.size(); ++ip) {
    unsigned klow  = std::get<0>(peakList[ip]);
    unsigned khigh = std::get<1>(peakList[ip]);
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
    phit.peakMaxTime=pmtEvent->time[peakt];
    phit.peakt=peakt;
    phit.startTime=pmtEvent->time[klow];
    phit.peakWidth=pmtEvent->time[khigh] - pmtEvent->time[klow];

    /* just use the biggest pulse */
    if(qsum>qmax) {
      qmax=qsum;
      firstTime=phit.startTime*1E6;
      firstCharge = qsum;
    }

    pmtHits.insert ( std::pair<Double_t,TPmtHit>(qsum,phit) );
    hPeakNWidth->Fill(phit.lastBin-phit.firstBin+1);
  }
  if(firstCharge<firstChargeCut&&pmtHits.size()>0&&qmax>firstChargeCut) printf("\t WARNING XXXXX NO FIRST PULSE pulses %i max %f \n",int(pmtHits.size()),qmax);
  //printf("\t XXXXX nhits %i first %f charge %f \n",int(pmtHits.size()),firstTime,firstCharge);

  return  pmtHits;
}

std::vector<std::complex<double> > anaRun::FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal)
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

std::vector< Double_t > anaRun::inverseFFT(Int_t ipmt,Int_t ievent,std::vector<std::complex<double> > VectorComplex,std::vector<Double_t> sum)
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

std::vector<Double_t> anaRun::BubbleSort(std::vector<Double_t> A){
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

std::vector<Double_t> anaRun::SimpleLowPassFilter(std::vector<Double_t> signal, Double_t alpha){
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

std::vector<Double_t> anaRun::SimpleHighPassFilter(std::vector<Double_t> signal, Double_t alpha){
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

std::vector<Double_t> anaRun::MovingAverageFilter(std::vector<Double_t> signal,Int_t N)
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

void anaRun::plotWave(Int_t ientry, Int_t pmtNum, hitMap pmtHits ) {
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
  TH1F* hbaseFit = (TH1F*) hBaselineFit[pmtNum]->Clone(histName);

  // title with first peak
  histName.Form("RawEv%i_PMT_%i",ientry,pmtNum);
  TH1F* rhist = (TH1F*) hPMTRaw[pmtNum]->Clone(histName);

  // title with first peak
  histName.Form("DerEv%i_PMT_%i",ientry,pmtNum);
  TH1F* dhist = (TH1F*) hPMTDerivative[pmtNum]->Clone(histName);
  

  // title with first peak
  histName.Form("WaveEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hist = (TH1F*) hPMTSignal[pmtNum]->Clone(histName);
  printf("plotWave: %s %.0f \n",hbase->GetName(),hbase->GetEntries());

  // fill peaks
  peaksName.Form("PeaksEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hpeaks = (TH1F*) hPMTSignal[pmtNum]->Clone(peaksName);
  hpeaks->Reset();
  for (hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
    TPmtHit phiti = hitIter->second;
    for(Int_t ibin=phiti.firstBin; ibin<=phiti.lastBin; ++ibin) hpeaks->SetBinContent(ibin, hist->GetBinContent(ibin));
  }
  
  if(isSimulation) {
    histName.Form("SimHitsEv%i_PMT_%i",ientry,pmtNum);
    TH1F* hist = (TH1F*) hPMTSim[pmtNum]->Clone(histName);
  }

}

// summed waves
void anaRun::sumWave(Int_t ipmt) {
  for (Int_t ibin=0; ibin<hSum[ipmt]->GetNbinsX(); ++ibin) {
    hSum[ipmt]->SetBinContent(ibin,  hSum[ipmt]->GetBinContent(ibin) + hPMTSignal[ipmt]->GetBinContent(ibin)  );
  }
}
// background estimator (class TSpectrum).
std::vector<Double_t> anaRun::getTBaseline(TH1D* hPMTRaw, TH1D* hBaseline, Double_t& ave, Double_t& aveSigma) 
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
    hBaseline->SetBinContent(i + 1,vbase[i]);
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
std::vector<Double_t> anaRun::getBaseline(std::vector<Double_t> digi, hitMap pmtHits, TH1D* hBaselineFit , TH1D* hBaseline, 
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
    hBaselineFit->SetBinContent(i + 1,vbase[i]);
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


void anaRun::getAverage(std::vector<Double_t> digi, Double_t& ave, Double_t& sigma) 
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
std::vector<Double_t> anaRun::differentiate(std::vector<Double_t> v, unsigned nstep)
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
peakType anaRun::derivativePeaks(std::vector<Double_t> v, Int_t nsum, Double_t rms) 
{
  peakType peakList;
  std::vector<unsigned> crossings;
  std::vector<unsigned> crossingBin;
  unsigned vsize = v.size();
  Double_t cut = 3.5*rms;
  Double_t ncut = -3.5*rms;
  // find all crossings
  for( unsigned ibin=1; ibin< vsize; ++ibin ) {
    if( v[ibin]>cut && v[ibin-1] <cut ) {
      crossings.push_back(UPCROSS);
      crossingBin.push_back(ibin);
    } else if( v[ibin]<cut && v[ibin-1] > cut ) {
      crossings.push_back(UPCROSS);
      crossingBin.push_back(ibin);
    } else if( v[ibin]<ncut && v[ibin-1] > ncut ) {
      crossings.push_back(DOWNCROSS);
      crossingBin.push_back(ibin);
    } else if( v[ibin]>ncut && v[ibin-1] < ncut ) {
      crossings.push_back(DOWNCROSS);
      crossingBin.push_back(ibin);
    }
  }

  // parse crossings to make pairs 
  unsigned ip =0; 
  //printf("crossings %zu \n",crossings.size());
  while( ip< crossings.size() -2 ) {
    if( ip<crossings.size()-3&&crossings[ip]==UPCROSS && crossings[ip+1]==UPCROSS && crossings[ip+2]==UPCROSS) {
      //printf(" peak %i (%i %i )\n",ip,crossings[ip],crossings[ip+1]); 
      peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
      ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(0));//sigma:d0:step:dstep
      ip = ip+2;
    } else if(ip<crossings.size()-4&&(crossings[ip]==UPCROSS&&crossings[ip+1]==UPCROSS&&crossings[ip+2]==DOWNCROSS&&crossings[ip+3]==DOWNCROSS)) {
     // printf(" peak %i (%i %i %i %i )\n",ip,crossings[ip],crossings[ip+1],crossings[ip+2],crossings[ip+3]); 
      peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+3]) );
      ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+3]-crossingBin[ip]),double(1));//sigma:d0:step:dstep
      ip=ip+4;
    } else if(ip==crossings.size()-2&&crossings[ip]==UPCROSS && crossings[ip+1]==UPCROSS) {
      //printf(" peak %i (%i %i )\n",ip,crossings[ip],crossings[ip+1]); 
      peakList.push_back( std::make_pair(crossingBin[ip],crossingBin[ip+1]) );
      ntDer->Fill(rms,v[crossingBin[ip]],double(crossingBin[ip+1]-crossingBin[ip]),double(0));//sigma:d0:step:dstep
      ++ip;
    } else {
      ++ip;
    }
  }

  // extend pulses to zero derivative
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
  // return list
  return peakList;
}

