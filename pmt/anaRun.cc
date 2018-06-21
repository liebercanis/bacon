/* MG revised */
#include "anaRun.hh"

anaRun::anaRun(TString tag)
{
  //Int_t irunStop = irunStart;
  TString outFileName ; outFileName.Form("baconRunAna_%s.root",tag.Data());
  TFile *outfile = new TFile(outFileName,"recreate");
  printf(" opening output file %s \n",outFileName.Data());
  
  ntupleCal = new TNtuple("ntcal","ntuple Cal","irun:ientry:baseline0:baseline1:sdev0:sdev1");
  ntupleRun = new TNtuple("ntrun","ntuple Run","run:npmt:entry:nhits:nbins:q:start:width:qmax");


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
    cout << " number of entries is " << nentries << endl;

    // set up memory for reading
    pmtEvent = new TPmtEvent();
    pmtTree->SetBranchAddress("pmtEvent", &pmtEvent);
    
    //switch to output file
    outfile->cd();

    Double_t maxBaseline = 0.05;
    TH1D *hBaseline1  = new TH1D(Form("Baseline1_%s",tag.Data()),"",100,-maxBaseline,maxBaseline);
    TH1D *hBaseline2  = new TH1D(Form("Baseline2_%s",tag.Data()),"",100,-maxBaseline,maxBaseline);
    TH1D *hNHits = new TH1D("NHits"," number of hits ",10,0,10);

    printf(" peak finding parameters: maxBaseline %.4f sigma %.4f \n",maxBaseline,sigma);

    double maxSample[2];
    int printInterval=1000;
    int nHists=0;
    int nHistsMax=100;
    int singlePhoton =0;

    // loop over entries
    for (UInt_t ientry=0; ientry<nentries; ientry++) {
      if(ientry%printInterval==0) printf(" ... %i ",ientry);
      pmtTree->GetEntry(ientry);
      if(pmtEvent->time.size()<1) continue;
      nSamples = pmtEvent->time.size();
      hBaseline1->Reset();
      hBaseline2->Reset();

      TString name; name.Form("run%s_ev%i",tag.Data(),ientry);
      TH1D* hPMTSignal[2];

      hPMTSignal[0] = new TH1D("PMTSignal1_" + name + TString("_baseline_subtracted"),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);
      hPMTSignal[1] = new TH1D("PMTSignal2_" + name + TString("_baseline_subtracted"),"",nSamples,pmtEvent->time[0],pmtEvent->time[nSamples-1]);

      // initialize fft 
      fFFT = TVirtualFFT::FFT(1, &nSamples, "R2C M K");
      fInverseFFT = TVirtualFFT::FFT(1, &nSamples, "C2R M K");
      std::vector<Double_t>  sdigi1,ddigi1,fftdigi1;
      std::vector<Double_t>  sdigi2,ddigi2,fftdigi2;

      std::vector<std::vector<Double_t> > ddigi;  
      ddigi.resize(NPMT);

      // loop over samples 
      for(Int_t isam = 0; isam < nSamples; isam++){
        Double_t volt1 = pmtEvent->volt1[isam];
        Double_t volt2 = pmtEvent->volt2[isam];
        Double_t digi0 = -1.0*(double(volt1));
        Double_t digi1 = -1.0*(double(volt2));
        ddigi[0].push_back(digi0);
        ddigi[1].push_back(digi1);
        hBaseline1->Fill(digi0);
        hBaseline2->Fill(digi1);
      }

      //Q is quite mode! very important
      hBaseline1->Fit("gaus","Q");
      if(!hBaseline1->GetFunction("gaus") ) continue;
      baseline[0] = hBaseline1->GetFunction("gaus")->GetParameter(1);
      sDev[0] = hBaseline1->GetFunction("gaus")->GetParameter(2);
      if(std::fabs(baseline[0]) > maxBaseline){
        cout << " WARNING LARGE BASELINE 1 " << ientry << "  " << baseline[0] << endl;
        delete hPMTSignal[0];
        delete hPMTSignal[1];
        continue; 
      }

      hBaseline2->Fit("gaus","Q");
      if(!hBaseline2->GetFunction("gaus") ) continue;
      baseline[1] = hBaseline2->GetFunction("gaus")->GetParameter(1);
      sDev[1] = hBaseline2->GetFunction("gaus")->GetParameter(2);
      if(std::fabs(baseline[1]) > maxBaseline){
        cout << " WARNING LARGE BASELINE 2 " << ientry << "  " << baseline[1] << endl;
        delete hPMTSignal[0];
        delete hPMTSignal[1];
        continue; 
      }

      for(int j = 0; j < ddigi.size();j++){
        maxSample[j]=0;
        for(int i = 0; i < ddigi[j].size(); i++){
          if(ddigi[j][i]>maxSample[j])  maxSample[j] = ddigi[j][i];
          ddigi[j][i] = ddigi[j][i]-baseline[j];// + sDev[j];
          hPMTSignal[j]->SetBinContent(i,ddigi[j][i]);
        }
      }

      if(ientry%printInterval==0) printf("  baseline %f %f sDev %f %f  ",baseline[0],baseline[1],sDev[0],sDev[1]);
      ntupleCal->Fill(0,ientry,baseline[0],baseline[1],sDev[0],sDev[1]);

      double deltaT = 1.0;

      for(int pmtNum = 0 ; pmtNum < NPMT; pmtNum++){
        qhitMax.clear();
        peakMaxTime.clear();
        peakBin.clear();
        qSum.clear();
        startTime.clear();
        peakWidth.clear();
        peakNbins.clear();

        Double_t minDev = 0*sDev[pmtNum], maxDev = sigma*sDev[pmtNum];
        std::vector<Int_t> peakTime = findPeaks(ddigi[pmtNum],maxDev,minDev); //,tStart,tStop);
        Int_t nhits = findHits( peakTime, ddigi[pmtNum],maxDev);
        if(ientry%printInterval==0) printf(" pmt  %i peaktime %lu nhits %i  ",pmtNum,peakTime.size(),nhits);
        hNHits->Fill(nhits);

        if(nhits<1) { 
          delete hPMTSignal[pmtNum];
          continue;
        }
        ++singlePhoton;

        if(nHists<nHistsMax){ 
          TString pmtSignalName;
          // title with first peak
          pmtSignalName.Form("PMT_%i_nhits_%i_charge_%.2f_start_%.2f_nbins_%i_nhit_%i",
              pmtNum, nhits, deltaT*qSum[0]*1E10, startTime[0], peakNbins[0], nhits);
          hPMTSignal[pmtNum]->SetTitle(pmtSignalName);
          cout << " hist count  " << ++nHists << "  " << pmtSignalName << endl;
        } else 
          delete hPMTSignal[pmtNum];

        for(unsigned ip=0; ip<qSum.size(); ++ip) ntupleRun->Fill(0,pmtNum,ientry,nhits,peakNbins[ip],deltaT*qSum[ip]*1E10,startTime[ip],peakWidth[ip]*1e6,qhitMax[ip]*1E10);

      }
      if(ientry%printInterval==0) printf("\n");
    }
    cout<<singlePhoton<<" events with hits out of "<<nentries<< " ntupleRun " << ntupleRun->GetEntries() << endl;
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


Int_t anaRun::findHits( std::vector<Int_t> peakTime, std::vector<Double_t> ddigi,Double_t sigma) 
{

  if(peakTime.size()<1) return 0;
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
    peakNbins.push_back(Int_t(hitTime.size()));
    qSum.push_back(qsum);
    qhitMax.push_back(qpeak);
    peakMaxTime.push_back(pmtEvent->time[peakt]);
    peakBin.push_back(peakt);
    startTime.push_back(pmtEvent->time[hitTime[hitTime.size() - 1] ] );
    peakWidth.push_back( pmtEvent->time[hitTime[0] ] - pmtEvent->time[hitTime[hitTime.size() - 1] ] );
    ++nhits;
  }

  return  nhits;
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



