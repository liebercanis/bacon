/* MG revised */
#include "anaRun.hh"

anaRun::anaRun(TString tag, Int_t maxEvents)
{
  printf(" starting anaRun tag %s \n",tag.Data());
  firstChargeCut=0.1;
  ran = new TRandom3();
  Double_t sigma=0;
  fsigma=sigma;
  if(sigma==0) fsigma=5;
  windowSize=3;
  nSigma=5;
  aveWidth=20;
  spec = new TSpectrum();
  //Int_t irunStop = irunStart;
  TString outFileName ; outFileName.Form("baconRunAna_%s_Ev_%i_window_baseline%i.root",tag.Data(),maxEvents,1);
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  ntCal = new TNtuple("ntCal","ntuple Cal","iev:ipmt:base:sigma");
  ntHit =   new TNtuple("ntHit","ntuple Hit","npmt:nhits:order:istart:time:tstart:q:nwidth:qmax");
  ntNHit =  new TNtuple("ntNHit","negative ntuple Hit","npmt:nhits:order:istart:time:tstart:q:nwidth:qmax");
  ntEvent = new TNtuple("ntEvent","ntuple Event","entry:n0:n1:t00:t01:t10:t11:qp0:qp1:q00:q01:q10:q11:qsum0:qsum1");
  ntPulse = new TNtuple("ntPulse"," pulse ","sum:shigh:slow:nsamp:kover:qlow:qhigh:klow:khigh");

  hPeakNWidth = new TH1I("PeakNWidth","PeakNWidth",100,0,100);
  hAllNWidth = new TH1I("AllNWidth","AllNWidth",100,0,100);
  hNWidthCut = new TH1I("NWidthCut","NWidth cut ",100,0,100);
  hSlideWidth = new TH1I("SlideWidth","sliding window pulse number Width",100,0,100);
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
  TString fileName; fileName.Form("rootData/baconRun_%s.root",tag.Data());
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
  printf(" number of entries is %lld sliding window peakfinding nSigma = %i \n",nentries,nSigma);

  // set up memory for reading
  pmtEvent = new TPmtEvent();
  pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);

  //switch to output file
  outfile->cd();

  //Event Vectors
  Int_t movingN1 = 0;
  Int_t movingN2 = 0;
  Double_t runningNoise1 = 0;
  Double_t runningNoise2 = 0;
  std::vector<Double_t> runningNoise(2);
  std::vector<Int_t> movingN(2);
  cout << " XXXXXXXX movingN size " << movingN.size() << endl;
  printf(" peak finding parameters: \n");
  printf(" \t sigma %.2f minLength %i  \n",fsigma, minLength); 

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
  int printInterval=10;
  if(maxEvents>0) nentries=maxEvents;
  for (UInt_t ientry=0; ientry<nentries; ientry++) {
    pmtTree->GetEntry(ientry);
    if(pmtEvent->time.size() == 0) continue;
    nSamples = pmtEvent->time.size();
    if(ientry==0) printf(" .... nSamples %i ....  ",nSamples);
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

    TString name; name.Form("%s_Ev%i",tag.Data(),ientry);

    int gotPMT = 1;

    // define pmt signal histograms
    if(ientry==0) {
      source = new Double_t[nSamples];
      Double_t maxLife = pmtEvent->time[nSamples-1]*1E6;
      printf(" \n\n ***** setting maxLife %f \n",maxLife);
      for(int ipmt=0; ipmt<gotPMT; ++ipmt) {
        hPMTRaw[ipmt] = new TH1D(Form("PMTRaw%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
        hPMTSignal[ipmt] = new TH1D(Form("PMTSignal%i_%s",ipmt,tag.Data()),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
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
      }
      // initialize fft 
      fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
      fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
    }


    std::vector<std::vector<Double_t> > ddigi; ddigi.resize(NPMT);
    std::vector<std::vector<Double_t> > ndigi; ndigi.resize(NPMT);
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

    // get baseline
    for(int ipmt=0; ipmt<gotPMT; ++ipmt) {
      baselineDigi[ipmt] = getBaseline(ddigi[ipmt],hBaselineFit[ipmt],hBaseline[ipmt],baseline[ipmt],sDev[ipmt]); 
      ntCal->Fill(ientry,ipmt,baseline[ipmt],sDev[ipmt]);
    }

    // subtract baseline
    for(int j = 0; j < gotPMT; j++){
      maxSample[j]=0;
      hPMTSignal[j]->Reset();
      for(int i = 0; i < int(ddigi[j].size()); i++){
        if(ddigi[j][i]>maxSample[j])  maxSample[j] = ddigi[j][i];
        ddigi[j][i] = ddigi[j][i]-baselineDigi[j][i];// + sDev[j];
        ndigi[j][i] = ndigi[j][i]+baselineDigi[j][i];// + sDev[j];
        hPMTSignal[j]->SetBinContent(i,ddigi[j][i]);
      }
    }

    if(ientry%printInterval==0) printf(" ...... %i fsigma %f fsigma*sdev0 %f \n",ientry,fsigma,sDev[0]);
    for(int pmtNum = 0 ; pmtNum < gotPMT; pmtNum++){
  
      /* 
      ** puluse finding and hit making
      */
      Double_t minDev = 0*sDev[pmtNum];
      Double_t maxDev = fsigma*sDev[pmtNum];
      Double_t firstTime, firstCharge;
      std::vector<Int_t> peakTime = windowPeaks(ddigi[pmtNum],windowSize,sDev[pmtNum]);
      hitMap  pmtHits = findHits(peakTime, ddigi[pmtNum],maxDev,firstTime,firstCharge);

      // negative pulses
      Double_t nfirstTime, nfirstCharge;
      std::vector<Int_t> npeakTime = windowPeaks(ndigi[pmtNum],windowSize,sDev[pmtNum]);
      hitMap  npmtHits = findHits(npeakTime, ndigi[pmtNum],maxDev,nfirstTime,nfirstCharge);


      unsigned nhits = pmtHits.size();
      hNHits[pmtNum]->Fill(nhits);
      unsigned negnhits = npmtHits.size();
      hNegNHits[pmtNum]->Fill(negnhits);

      if(ientry%printInterval==0) printf(" pmt  %i peaktime %lu nhits %u  ",pmtNum,peakTime.size(),nhits);

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


      if(nHists<nHistsMax) {
        plotWave(ientry,pmtNum,pmtHits );
        ++nHists;
      }

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
      ntEvent->Fill(ientry,npmtHit[0],npmtHit[1],tpmt[0][0],tpmt[0][1],tpmt[1][0],tpmt[1][1],qped[0],qped[1],qpmt[0][0],qpmt[0][1],qpmt[1][0],qpmt[1][1],qsum[0],qsum[1]);
  }
  printf(" total hits %i, %i out of %lld events\n",totalHits[0],totalHits[1],nentries);

  outfile->Purge();
  outfile->Write();

  return;
}

std::vector<Int_t> anaRun::findPeaks(std::vector<Double_t> v, Double_t threshold, Double_t sthreshold) 
{
  // Produces a list of peaks above the threshold
  std::vector<Int_t> peakTime;
  Int_t klow=0;
  Int_t khigh=0;
  Int_t kover=0;
  Int_t vsize = Int_t(v.size());

  //printf(" findPeaks \n");
  for( Int_t ibin=0; ibin< vsize; ++ibin ) {
    //if(ibin%100==0) printf("iiiiiii  ibin %i v %f threshold %f \n",ibin,v[ibin],threshold);
    if( v[ibin]>threshold && v[ibin+1] && v[ibin+2]>threshold) {// starting possible new hit
      // consider this a "seed" and find full hit
      klow=ibin;
      for(Int_t k=ibin-1; k>=max(0,ibin-maxHalfLength); --k) {
        if(v[k]<=sthreshold) break;
        klow=k;
      }
      khigh=ibin;
      for(Int_t k=ibin+1; k<min(ibin+maxHalfLength,vsize); ++k) {
        if(v[k]<=sthreshold) break;
        khigh=k;
      }
      kover = khigh-klow+1;
      // found good pulse
      hAllNWidth->Fill(kover);
      if(kover>minLength) {
        double qsum=0;
        for(Int_t k=klow ; k<= khigh; ++k) {
          peakTime.push_back(k);
          qsum +=v[k];
        }
        //printf(" ***** kover = %i peakTime  qsum %f sthreshod %f ibin %i klow %i khigh %i ", kover, qsum,sthreshold,ibin,klow,khigh);
        //for(Int_t k=klow ; k<= khigh; ++k) printf(" v[%i]=%f ",k,v[k]);
        //printf("\n");
      } 
      // skip to end of sthreshold search 
      ibin=khigh;
    }
  }
  

  //for(UInt_t ih=0; ih<peakTime.size(); ++ih) printf("  %i t= %i ADC= %f\n",ih,peakTime[ih],v[peakTime[ih]]);
  return peakTime;
}


/* 
  A sliding integration window of 12 ns (3 samples) is used to scan the calibrated waveforms to identify the pulse. 
  The pulse is extracted whenever the integral exceeds 5 times the RMS of noise samples times the square root of the 
  number of samples in the window. The boundaries of the pulse region is defined when the sliding window integral drops 
  below the RMS times divided by the square root of the number of samples. 
*/

std::vector<Int_t> anaRun::windowPeaks(std::vector<Double_t> v, Int_t nsum, Double_t rms) 
{
  std::vector<Int_t> peakTime;
  Int_t vsize = Int_t(v.size());

  //printf(" windowPeaks \n");
  for( Int_t ibin=0; ibin< vsize; ++ibin ) {
    Int_t klow=0;
    Int_t khigh=0;
    Int_t kover=0;
    Int_t nsamples = nsum;
    Double_t slow  = rms/sqrt(Double_t(nsum));
    // starting bin has to be greater than baseline
    if( v[ibin]<slow ) continue;

    Double_t sum = 0;
    Int_t imax=0;
    Double_t vmax=0;
    for(Int_t isum = ibin; isum < min(ibin+nsum,vsize) ; ++isum) {
      sum+= v[isum];
      if(v[isum]>vmax) {
        vmax = v[isum];
        imax = isum;
      }
    }

    Double_t shigh = Double_t(nSigma)*rms*sqrt(Double_t(nsum));
    hSlideHigh->Fill(shigh);

    if( sum>shigh)  {// starting hit
      klow=imax;
      for(Int_t k=imax-1; k>=max(0,imax-maxHalfLength); --k) { // sum to low edge
        Double_t tempSum = sum+v[k];
        slow  = rms/sqrt(Double_t(nsamples+1));
        if(sum<slow||v[k]<0) break;
        ++nsamples;
        sum=tempSum;
        hSlideLow->Fill(slow);
        klow=k;
      }
      khigh=imax;
      for(Int_t k=imax+1; k<min(imax+maxHalfLength,vsize); ++k) {  // sum to high edge
        Double_t tempSum = sum+v[k];
        slow  = rms/sqrt(Double_t(nsamples+1));
        if(sum<slow||v[k]<0) break;
        sum=tempSum;
        ++nsamples;
        hSlideLow->Fill(slow);
        khigh=k;
      }
      kover = khigh-klow+1;
      hSlideWidth->Fill(kover);
      if(kover<nsum) continue;
      // found pulse
      for(Int_t k=klow ; k<= khigh; ++k) peakTime.push_back(k);
      if(v[klow]<0||v[khigh]<0) printf(" ***** kover = %i peakTime  qsum %f ibin %i k[%i]= %E k[%i]= %E \n)", kover, sum,ibin,klow,v[klow],khigh,v[khigh]);
      //for(Int_t k=klow ; k<= khigh; ++k) printf(" v[%i]=%f ",k,v[k]);
      //printf("\n");
      hSlideQSum->Fill(sum);
      // skip to end of sthreshold search 
      ibin=khigh;
      ntPulse->Fill(sum,shigh,slow,float(nsamples),float(kover),v[klow],v[khigh],float(klow),float(khigh));
    }
  }
  //
  return peakTime;
}


hitMap anaRun::findHits(std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma, Double_t& firstTime, Double_t& firstCharge) 
{

  firstTime=1E9;
  firstCharge=0;
  hitMap pmtHits;
  if(peakTime.size()<1) return pmtHits;

  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(Int_t it=nlast; it>0; --it) {
    bool makeHit=false;
    if(peakTime[it]-peakTime[it-1]!=1||(it==1&&hitTime.size()>=minLength)) makeHit=true;
    //if(peakTime[it]-peakTime[it-1]!=1) makeHit=true;

    if(makeHit) {
      hitTime.push_back(peakTime[it]);
      hitList.push_back(hitTime);
      hitTime.clear();
      continue;
    }
    hitTime.push_back(peakTime[it]);
  }

  Int_t nhits=0;
  Double_t qmax=0;

  for(UInt_t il=0; il<hitList.size(); ++il) {
    hitTime=hitList[il];
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0;
    Double_t qsum = 0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      //printf(" \t ih = %i time(us)  %f sample %f  \n ",hitTime[ih],1e6*pmtEvent->time[hitTime[ih] ],ddigi[hitTime[ih]]);
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
      }
      qsum += ddigi[hitTime[ih]];
    }

    /* remove this cut 
    if (qpeak < sigma){
      continue;
    }
    */

    TPmtHit phit;
    phit.peakBin=Int_t(hitTime.size());
    phit.qsum=qsum;
    phit.qpeak=qpeak;
    phit.firstBin = hitTime[hitTime.size() - 1];
    phit.lastBin = hitTime[0];
    phit.peakMaxTime=pmtEvent->time[peakt];
    phit.peakt=peakt;
    phit.startTime=pmtEvent->time[hitTime[hitTime.size() - 1] ];
    phit.peakWidth=pmtEvent->time[hitTime[0] ] - pmtEvent->time[hitTime[hitTime.size() - 1] ];

    /* just use the biggest pulse */
    if(qsum>qmax) {
      qmax=qsum;
      firstTime=phit.startTime*1E6;
      firstCharge = qsum;
    }

    /** find first light pulse > cut  
    if(qsum>firstChargeCut&&phit.startTime*1E6<firstTime) {
      firstTime=phit.startTime*1E6;
      firstCharge = qsum;
      //printf(" \t setting firstTime %f firstCharge %f \n ",firstTime,firstCharge);
    }
    */

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
    hfft->SetBinContent(i,hfft->GetBinContent(i)+std::abs(c));
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

  /* 
     for(int j = 0 ; j <dsigi.size(); j++){
     for(int i = 0; i < 4; i++){ 
     dsigi[j] = MovingAverageFilter(dsigi[j],5);
     }
     }
     for(int j = 0; j < dsigi.size(); j++){
     for(int i = 0; i < dsigi[j].size() ; i++){
     if(j == 0)
     hPMTSignalFiltered1->SetBinContent(i,dsigi[j][i]);
     else
     hPMTSignalFiltered2->SetBinContent(i,dsigi[j][i]);
     }
     }
     */
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
  histName.Form("WaveEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hist = (TH1F*) hPMTSignal[pmtNum]->Clone(histName);
 
  printf("plotWave: %s %f \n",hbase->GetName(),hbase->GetEntries());

  // fill peaks
  peaksName.Form("PeaksEv%i_PMT_%i",ientry,pmtNum);
  TH1F* hpeaks = (TH1F*) hPMTSignal[pmtNum]->Clone(peaksName);
  hpeaks->Reset();
  for (hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
    TPmtHit phiti = hitIter->second;
    for(Int_t ibin=phiti.firstBin; ibin<=phiti.lastBin; ++ibin) hpeaks->SetBinContent(ibin, hist->GetBinContent(ibin));
  }

  /*/ sum hist title 
  sumName.Form("SumEv_%i_PMT_%i_Q_%.2E",int(ientry),pmtNum, phit0.qsum);
  TH1F* hsum = (TH1F*) hPMTSignal[pmtNum]->Clone(sumName);
  hsum->Reset();
  Double_t qsum=0;
  for (Int_t ibin=0; ibin<hsum->GetNbinsX(); ++ibin) {
    qsum +=  hpeaks->GetBinContent(ibin);
    hsum->SetBinContent(ibin, qsum);
  }
  */

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

// background estimator removing pulses 
std::vector<Double_t> anaRun::getBaseline(std::vector<Double_t> digi, TH1D* hBaselineFit , TH1D* hBaseline, Double_t& ave, Double_t& aveSigma) 
{
  std::vector<Double_t> vbase;
  vbase.resize(digi.size());
  hBaseline->Reset();
  hBaselineFit->Reset();
  // get first estimate of noise, baseline.
  Double_t vave, vsigma;
  getAverage(digi,vave,vsigma);

  // subtract baseline estimate
  std::vector<Double_t> vsignal;
  for(unsigned is=0; is< digi.size() ; ++is) vsignal.push_back(digi[is]-vave);

  std::vector<Int_t> peakTime = windowPeaks(vsignal,windowSize,vsigma);

  for(unsigned is=0; is< vsignal.size(); ++is) {
    vbase[is]=vsignal[is]+vave;  // add back the average to the baseline
    for(unsigned ip=0; ip< peakTime.size(); ++ip) if(peakTime[ip]==int(is)) vbase[is]=ran->Gaus(vave,vsigma);
  }

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
  ave = vave+fpoly->GetParameter(0); //value of 1st parameter
  // fill baseline function from fit.
  for (int i = 0; i <  hBaselineFit->GetNbinsX(); i++) {
    Double_t xbin = hBaselineFit->GetXaxis()->GetBinCenter(i+1);
    Double_t vbin = fpoly->Eval(xbin);  
    vbase[i]=vbin;
  }

  // recalculate ave,sigma
  getAverage(vbase,ave,aveSigma);

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

  /*
  // summed 
  Double_t sum=0,sum2=0;
  for(unsigned is=0; is< digi.size() ; ++is) {
    sum  += digi[is];
  }

  sum  /= sum/Double_t(digi.size());
  for(unsigned is=0; is< digi.size() ; ++is) {
    sum2 += pow(digi[is]-sum,2.0);
  }

  sum2 /= sum2/Double_t(digi.size());

  Double_t rms = sum2;
  if(rms>0) rms = sqrt(rms);

  printf(" ave %E (%E) rms %E (%E) \n",ave,sum,sigma,rms);
  */

}

