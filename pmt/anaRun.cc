/* MG revised */
#include "anaRun.hh"

anaRun::anaRun(TString tag)
{
  Double_t maxBaseline = 0.01;
  //Int_t irunStop = irunStart;
  TString outFileName ; outFileName.Form("baconRunAna_%s.root",tag.Data());
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());

  ntCal = new TNtuple("ntCal","ntuple Cal","irun:ientry:baseline0:baseline1:sdev0:sdev1");
  ntHit = new TNtuple("ntHit","ntuple Hit","npmt:nhits:order:nbins:q:start:width:qmax");
  ntEvent = new TNtuple("ntEvent","ntuple Event","entry:n0:n1:t00:t01:t10:t11:q00:q01:q10:q11");


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
  printf(" number of entries is %lld peakfinding sigma = %.2f minlength = %i \n",nentries,sigma,minLength);

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
  TH1D *hBaseline1  = new TH1D(Form("Baseline1_%s",tag.Data()),"",100,-maxBaseline,maxBaseline);
  TH1D *hBaseline2  = new TH1D(Form("Baseline2_%s",tag.Data()),"",100,-maxBaseline,maxBaseline);
  TH1D *hNHits = new TH1D("NHits"," number of hits ",10,0,10);

  printf(" peak finding parameters: \n");
  printf(" \t sigma %.2f minLength %i  \n",sigma, minLength); 

  double maxSample[2];
  int printInterval=1000;
  int nHists=0;
  int nHistsMax=100;
  TH1D* hPMTSignal[2];
  Float_t qpmt[2][2];
  Float_t tpmt[2][2];
  Int_t   npmtHit[2];
  // loop over entries
  for (UInt_t ientry=0; ientry<nentries; ientry++) {
    if(ientry%printInterval==0) printf(" ... %i ",ientry);
    pmtTree->GetEntry(ientry);
    if(pmtEvent->time.size() == 0) continue;
    nSamples = pmtEvent->time.size();
    hBaseline1->Reset();
    hBaseline2->Reset();
    for(int ipmt=0; ipmt<2 ; ++ ipmt) {
      npmtHit[ipmt]=0;
      for(int ihit=0; ihit<2 ; ++ihit) {
        qpmt[ipmt][ihit]=0;
        tpmt[ipmt][ihit]=0;
      }
    }

    TString name; name.Form("run%s_ev%i",tag.Data(),ientry);

    /*
       if( ientry == 0){ 
       hPMTSum1 = new TH1D("PMT_Sum_0" + name,"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
       hPMTSum2 = new TH1D("PMT_Sum_1" + name,"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
       }
       */

    if(ientry==0) {
      hPMTSignal[0] = new TH1D("PMTSignal1_" + name + TString("_baseline_subtracted"),"",
          nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
      hPMTSignal[1] = new TH1D("PMTSignal2_" + name + TString("_baseline_subtracted"),"",
          nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
    }

    // initialize fft 
    fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
    fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
    std::vector<Double_t>  sdigi1,ddigi1,fftdigi1;
    std::vector<Double_t>  sdigi2,ddigi2,fftdigi2;
    Double_t baseline1 = 0,sum1 = 0,T1 = 0,T2 = 0;
    Double_t baseline2 = 0,sum2 = 0;


    //For cal ntuple
    //Int_t tStart = hPMTSignal1->FindBin(-50e-9);
    //Int_t tStop = hPMTSignal1->FindBin(50e-9);

    Int_t tStart = 0;//hPMTSignal1->FindBin(0e-9);
    Int_t tStop = nSamples;//hPMTSignal1->FindBin(500e-9);
    Int_t tStartZeroWidth = hPMTSignal[0]->FindBin(0e-9);
    Int_t tStopZeroWidth = hPMTSignal[0]->FindBin(50e-9);

    if(ientry == 0)  cout<<"StartTime = "<<hPMTSignal[0]->GetBinCenter(tStart)<<", StopTime = "<< hPMTSignal[0]->GetBinCenter(tStop)<<endl;


    std::vector<std::vector<Double_t> > ddigi;ddigi.resize(NPMT);

    deltaT = 0;
    // loop over samples 
    for(Int_t isam = 0; isam < nSamples; isam++){
      Double_t volt0 = pmtEvent->volt1[isam];
      Double_t volt1 = pmtEvent->volt2[isam];
      Double_t digi0 = -1.0*(double(volt0));
      Double_t digi1 = -1.0*(double(volt1));
      ddigi[0].push_back(digi0);
      ddigi[1].push_back(digi1);
      hBaseline1->Fill(digi0);
      hBaseline2->Fill(digi1);
      if(isam > 0) deltaT += std::fabs(pmtEvent->time[isam] - pmtEvent->time[isam - 1]);
    }

    deltaT /= (nSamples - 1);

    //Q is quite mode! very important
    TFitResultPtr ptr1 = hBaseline1->Fit("gaus","Q");
    baseline[0]=0;
    sDev[0]=0;
    if(!hBaseline1->GetFunction("gaus")) {
      printf(" \t \t cannot fit 1 event %u  \n",ientry);
      (TH1F*) hBaseline1->Clone(Form("%s-event%u-nofit",hBaseline1->GetName(),ientry));
    } else {
      baseline[0] =  hBaseline1->GetFunction("gaus")->GetParameter(1);
      sDev[0] = hBaseline1->GetFunction("gaus")->GetParameter(2);
    }
    TFitResultPtr ptr2 = hBaseline2->Fit("gaus","Q");
    baseline[1]=0;
    sDev[1]=0;
    if(!hBaseline2->GetFunction("gaus")) {
      printf(" \t \t cannot fit 2 event %u \n",ientry);
      (TH1F*) hBaseline2->Clone(Form("%s-event%u-nofit",hBaseline2->GetName(),ientry));
    } else {
      baseline[1] = hBaseline2->GetFunction("gaus")->GetParameter(1);
      sDev[1] = hBaseline2->GetFunction("gaus")->GetParameter(2);
    }
    ntCal->Fill(0,ientry,baseline[0],baseline[1],sDev[0],sDev[1]);
    if(std::fabs(baseline[0]) > maxBaseline ||std::fabs(baseline[1]) > maxBaseline ){
      printf(" WARNING LARGE BASELINE  event %u base 1 %f base 2 %f \n",ientry,baseline[0],baseline[1]);
    }

    for(int j = 0; j < ddigi.size();j++){
      maxSample[j]=0;
      hPMTSignal[j]->Reset();
      for(int i = 0; i < ddigi[j].size(); i++){
        if(ddigi[j][i]>maxSample[j])  maxSample[j] = ddigi[j][i];
        ddigi[j][i] = ddigi[j][i]-baseline[j];// + sDev[j];
        hPMTSignal[j]->SetBinContent(i,ddigi[j][i]);
      }
    }

    if(ientry%printInterval==0) printf(" deltaT %.2E baseline %f %f sDev %f %f  ",deltaT,baseline[0],baseline[1],sDev[0],sDev[1]);


    Double_t T0_temp = 1;
    std::vector<bool> nohitFlag;nohitFlag.resize(2);
    

    for(int pmtNum = 0 ; pmtNum < NPMT; pmtNum++){
      //for(double sigma = 3 ; sigma <= 4; sigma += 0.2){

      Double_t minDev = 0*sDev[pmtNum], maxDev = sigma*sDev[pmtNum];
      std::vector<Int_t> peakTime = findPeaks(ddigi[pmtNum],maxDev,minDev); //,tStart,tStop);
      hitMap  pmtHits = findHits( peakTime, ddigi[pmtNum],maxDev);
      unsigned nhits = pmtHits.size();
      if(ientry%printInterval==0) printf(" pmt  %i peaktime %lu nhits %u  ",pmtNum,peakTime.size(),nhits);

      hNHits->Fill(nhits);
      if(nhits<1) continue;

      npmtHit[pmtNum]=nhits;
      singlePhoton++;
      TString pmtSignalName;
      hitMapIter hitIter;
      hitIter=pmtHits.begin();
      TPmtHit phit0 = hitIter->second;
      TPmtHit phit1;
      // next hit if it is there
      if(pmtHits.size()>1) {
        ++hitIter;
        phit1 = hitIter->second;
      }

      if(nHists<nHistsMax) {
        // title with first peak
        pmtSignalName.Form("Ev_%i_PMT_%i_charge_%.2f_peak_%.2f_peakt_%.2f_start_%.2f_width_%.2f_nhit_%i",
            int(ientry),pmtNum, 
            deltaT*phit0.qsum*1E10, 
            phit0.qpeak*1E10,
            phit0.peakBin, 
            phit0.startTime,
            phit0.peakWidth,
            nhits);
        (TH1F*) hPMTSignal[pmtNum]->Clone(pmtSignalName);
        cout << pmtSignalName << endl;
        ++nHists;
      } 
      int hitCount=0;
      for (hitIter=pmtHits.begin(); hitIter!=pmtHits.end(); ++hitIter) {
        TPmtHit phiti = hitIter->second;
        ntHit->Fill(pmtNum,nhits,hitCount++,phiti.peakBin,deltaT*phiti.qsum*1E10,phiti.startTime*1E6,phiti.peakWidth*1E6,phiti.qpeak*1E10);
      }

      // save for npmtHit ntuple
      qpmt[pmtNum][0]=deltaT*phit0.qsum*1E10;
      tpmt[pmtNum][0]=phit0.startTime*1E6;
      qpmt[pmtNum][1]=deltaT*phit1.qsum*1E10;
      tpmt[pmtNum][1]=phit1.startTime*1E6;

    }
    if(ientry%printInterval==0) printf(" single Photon %.0f \n",singlePhoton);
    
    if(npmtHit[0]>0||npmtHit[1]>0) 
      ntEvent->Fill(ientry,npmtHit[0],npmtHit[1],tpmt[0][0],tpmt[0][1],tpmt[1][0],tpmt[1][1],qpmt[0][0],qpmt[0][1],qpmt[1][0],qpmt[1][1]);

    }
    cout<<singlePhoton<<" pmts with hits out of "<<nentries<<endl;
    //}//irun loop

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
    if( v[ibin]>threshold) {// starting possible new hit
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


hitMap anaRun::findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma) 
{

  hitMap pmtHits;
  if(peakTime.size()<1) return pmtHits;

  std::vector<Int_t> hitTime;
  std::vector<std::vector<Int_t> > hitList;
  UInt_t nlast = peakTime.size()-1;
  for(Int_t it=nlast; it>0; --it) {
    bool makeHit=false;
    if(peakTime[it]-peakTime[it-1]!=1||(it==1&&hitTime.size()>=minLength)) makeHit=true;

    if(makeHit) {
      hitTime.push_back(peakTime[it]);
      hitList.push_back(hitTime);
      hitTime.clear();
      continue;
    }
    hitTime.push_back(peakTime[it]);
  }

  Int_t nhits=0;

  for(UInt_t il=0; il<hitList.size(); ++il) {
    hitTime=hitList[il];
    Double_t qhit=0;
    UInt_t peakt=0;
    Double_t qpeak=0,qsum = 0;
    Double_t qUnpeak=0;
    for(UInt_t ih=0; ih<hitTime.size(); ++ih) {
      //printf(" \t ih = %i time(us)  %f sample %f  \n ",hitTime[ih],1e6*pmtEvent->time[hitTime[ih] ],ddigi[hitTime[ih]]);
      if(ddigi[hitTime[ih]]>qpeak) {
        peakt=hitTime[ih];
        qpeak = ddigi[hitTime[ih]];
      }
      qsum += ddigi[hitTime[ih]];
    }
    if (qpeak < sigma){
      continue;
    }
    TPmtHit phit;
    phit.peakBin=Int_t(hitTime.size());
    phit.qsum=qsum;
    phit.qpeak=qpeak;
    phit.peakMaxTime=pmtEvent->time[peakt];
    phit.peakt=peakt;
    phit.startTime=pmtEvent->time[hitTime[hitTime.size() - 1] ];
    phit.peakWidth=pmtEvent->time[hitTime[0] ] - pmtEvent->time[hitTime[hitTime.size() - 1] ];
    pmtHits.insert ( std::pair<Double_t,TPmtHit>(qsum,phit) );
  }

  return  pmtHits;
}

std::vector<std::complex<double> > anaRun::FFT(Int_t ipmt,Int_t ievent,std::vector<Double_t> signal)
{
  /*
     for(int i = 0; i < ddigi.size(); i++){
     std::vector<std::complex<double> > fftPair = FFT(i,ientry,ddigi[i]);
     std::vector<Double_t> vecTemp = inverseFFT(i,ientry,fftPair,sum);
     }
     */  
  std::vector<std::complex<double> > VectorComplex;
  for(int is =0; is<nSamples; ++is) {
    if(ipmt==0) fFFT->SetPoint(is, signal[is]);
    else fFFT->SetPoint(is, signal[is]);
  }
  fFFT->Transform();
  if(ipmt == 0){
    hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),nSamples/2,0,nSamples/2);
    //hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),nSamples,0,nSamples);
  }
  else if(ipmt == 1){
    hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),nSamples/2,0,nSamples/2);
    //hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),nSamples,0,nSamples);
  }

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
  if(ipmt == 0){
    hIfft = new TH1D(Form("inverseFFTChannelOne_%i",ievent),Form("InverseFFT_One_%i",ievent),nSamples,pmtEvent->time[0],pmtEvent->time[nSamples - 1]);
    //hfft = new TH1D(Form("FFTChannelOne_%i",ievent),Form("FFT Channel One %i",ievent),nSamples,0,nSamples);
  }
  else if(ipmt == 1){
    hIfft = new TH1D(Form("inverseFFTChannelTwo_%i",ievent),Form("InverseFFT_Two_%i",ievent),nSamples,pmtEvent->time[0],pmtEvent->time[nSamples - 1]);
    //hfft = new TH1D(Form("FFTChannelTwo_%i",ievent),Form("FFT Channel Two %i",ievent),nSamples,0,nSamples);
  }
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
  for(int i = 0; i < signal.size(); i++){
    if(i == 0) FilteredSignal.push_back(alpha*signal[i]);
    else{
      FilteredSignal.push_back(alpha*signal[i] + (1.-alpha)*FilteredSignal[i-1] );
    }
  }

  //cout<<"Filter Size "<<FilteredSignal.size()<<", Ddigi size "<<signal.size()<<endl;
  return FilteredSignal;
}

std::vector<Double_t> anaRun::SimpleHighPassFilter(std::vector<Double_t> signal, Double_t alpha){
  std::vector<Double_t> FilteredSignal;
  //Double_t alpha = 0.992105;//100 MHz cut off frequency
  for(int i = 0; i < signal.size(); i++){
    if(i == 0) FilteredSignal.push_back(signal[i]);
    else{
      FilteredSignal.push_back( alpha*( signal[i]-signal[i-1] + FilteredSignal[i-1] )   );
    }
  }

  //cout<<"Filter Size "<<FilteredSignal.size()<<", Ddigi size "<<signal.size()<<endl;
  return FilteredSignal;

}

std::vector<Double_t> anaRun::MovingAverageFilter(std::vector<Double_t> signal,Int_t N)
{

  /* 
     for(int j = 0 ; j <ddigi.size(); j++){
     for(int i = 0; i < 4; i++){ 
     ddigi[j] = MovingAverageFilter(ddigi[j],5);
     }
     }
     for(int j = 0; j < ddigi.size(); j++){
     for(int i = 0; i < ddigi[j].size() ; i++){
     if(j == 0)
     hPMTSignalFiltered1->SetBinContent(i,ddigi[j][i]);
     else
     hPMTSignalFiltered2->SetBinContent(i,ddigi[j][i]);
     }
     }
     */
  std::vector<Double_t> filter;
  Int_t N2 = std::floor(N/2);
  for(int i = N2; i < unsigned(signal.size())-N2; i++){
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



