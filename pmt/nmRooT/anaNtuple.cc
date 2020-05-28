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
#include <TRandom3.h>

void anaNtuple(){
  
  Double_t irunStart = 0;
  TRandom3 *ran = new TRandom3();
  //TString filename = TString("cosmicMuonsThinPanel1PMTSmallRadius.root");
  //TString filename = TString("BACoNAnaOutputAllMaps.root");
  TString filename = TString("BACoNAnaOutput.root");
  TFile *f = new TFile(filename);
  cout<<"opening file "<<filename<<endl;
  if(f->IsZombie()){
    cout<<"cannot open file"<<endl;
  }
  TTree *t1 = (TTree*)f->Get("ntupleLAr");

  Float_t scint,wls,ceren,edep,nDetected,nDetectedOpticalMap,XenonDoped,uPanel1,uPanel2,muonTrigger,
          XenonDoped1000m,XenonDoped50m,XenonDoped20m,XenonDoped10m,XenonDoped5m,XenonDoped3m,XenonDoped1m;
  t1->SetBranchAddress("scint",&scint);
  t1->SetBranchAddress("wls",&wls);
  t1->SetBranchAddress("ceren",&ceren);
  t1->SetBranchAddress("edep",&edep);
  t1->SetBranchAddress("nDetected",&nDetected);
  //t1->SetBranchAddress("nDetectedOpticalMap",&nDetectedOpticalMap);
  //t1->SetBranchAddress("XenonDoped",&XenonDoped);
  t1->SetBranchAddress("LAr0.60m",&nDetectedOpticalMap);
  t1->SetBranchAddress("XenonDoped1000m",&XenonDoped1000m);
  t1->SetBranchAddress("XenonDoped50m",&XenonDoped50m);
  t1->SetBranchAddress("XenonDoped20m",&XenonDoped20m);
  t1->SetBranchAddress("XenonDoped10m",&XenonDoped10m);
  t1->SetBranchAddress("XenonDoped5m",&XenonDoped5m);
  t1->SetBranchAddress("XenonDoped3m",&XenonDoped3m);
  t1->SetBranchAddress("XenonDoped1m",&XenonDoped1m);
  t1->SetBranchAddress("uPanel1",&uPanel1);
  t1->SetBranchAddress("uPanel2",&uPanel2);
  t1->SetBranchAddress("muonTrigger",&muonTrigger);
  //scaling data
  cout<<"looking a background data"<<endl;
  TFile *f2 = new TFile("runAna9999_9999_DS_2.root");
  TH1D *hBackground;
  f2->GetObject("Charge_Run9999",hBackground);
  hBackground->SetNameTitle("XeD10ppmBackground","XeD10ppmBackground");
  Double_t maxBackground = hBackground->GetMaximum();
  //f2->Close();
  cout<<"Looking a muon data"<<endl;
  TFile *f3 = new TFile("runAna20220_20222_DS_2.root");
  TH1D* hXeDMuons;
  f3->GetObject("TotalChargeAverage",hXeDMuons);
  Double_t maxMuon = hXeDMuons->GetMaximum();
  //f3->Close();
  hBackground->Scale(maxMuon/maxBackground);

  
  cout<<"Doing a background subtraction, maxBackground = "<< maxBackground<<", max Muon "<< maxMuon<<endl;
  TH1D* hSubtraction = new TH1D("backgroundSubtractedXe10ppm","backgroundSubtractedXe10ppm",1000,10,10000);
  hSubtraction->Add(hXeDMuons,1);
  hSubtraction->Add(hBackground,-1);
  for(int i = 0; i < hSubtraction->GetNbinsX();i++){
    if(hSubtraction->GetBinContent(i+1) < 0)
      hSubtraction->SetBinContent(i+1,0);
  }
  hSubtraction->GetXaxis()->SetRangeUser(10,3000);
  //cout<<"Integral "<<hSubtraction->Integral()<<endl;
  hSubtraction->Scale(1./hSubtraction->Integral());

  //XeD 10 ppm 1 ppm N2
  //TFile* f88 = new TFile("runAna30101_30101_DS_2.root");
  TFile* f88 = new TFile("runAna30203_30208_DS_2.root");
  TH1D * hXeDN2Muon;
  f88->GetObject("TotalChargeAverage",hXeDN2Muon);
  hXeDN2Muon->Add(hXeDMuons,1);
  hBackground->Scale(hXeDN2Muon->GetMaximum()/hBackground->GetMaximum());
  hXeDN2Muon->Add(hBackground,-1);

  hXeDN2Muon->SetTitle("1ppmN2in10ppmXed");
  hXeDN2Muon->Scale(1./hXeDN2Muon->Integral());
  hXeDN2Muon->GetXaxis()->SetRangeUser(10,3000);


  //Liquid Argon Pure
  TF1 * fFit = new TF1("pureFit","[0]*exp(x*[1])+[2]*TMath::Landau(x,[3],[4])+[5]*exp(x*[6])",10,3000);
  fFit->SetParameter(0,3.64420e+04);
  fFit->SetParLimits(0,3.24420e+04,3.84420e+04);  

  fFit->SetParameter(1,-4.27476e-02);
  fFit->SetParLimits(1,-4.47476e-02,-4.07476e-02);
  
  fFit->SetParameter(2,7.53751e+03);
  fFit->SetParLimits(2,7.33751e+03,7.73751e+03);
  
  fFit->SetParameter(3,3.26046e+02);
  fFit->SetParLimits(3,300,350);
  
  fFit->SetParameter(4,7.14339e+01);
  fFit->SetParLimits(4,6.94339e+01,7.34339e+01);
  
  fFit->SetParameter(5,1.16342e+03);
  fFit->SetParLimits(5,0.916342e+03,1.36342e+03);
  
  fFit->SetParameter(6,-1.08294e-03);
  fFit->SetParLimits(6,-1.38294e-03,-0.88294e-03);

  TF1 * fExponetial = new TF1("expoPureLAr","[0]*exp(x*[1])",10,3000);

  TFile * f4 = new TFile("Pure_LAr.root");
  TH1D* hPureLAr = new TH1D("PureLArData","PureLArData",300,10,3000);
  f4->GetObject("TotalChargeAverage",hPureLAr);
  //hPureLAr->Fit("pureFit","0","",75,10000);
  hPureLAr->Fit("expoPureLAr","","",800,3000);
  hPureLAr->SetNameTitle("PureLAr","PureLAr");

  //TF1 * fBackground = new TF1("pureBackground","[0]*exp(x*[1])+[2]*exp(x*[3])");
  TF1 * fBackground = new TF1("pureBackground","[0]*exp(x*[1])");
  /*
  fBackground->SetParameter(0,fFit->GetParameter(0));
  cout<<"para 0 "<<fBackground->GetParameter(0)<<endl;;
  fBackground->SetParameter(1,fFit->GetParameter(1));
  cout<<"para 1 "<<fBackground->GetParameter(1)<<endl;
  fBackground->SetParameter(0,fFit->GetParameter(5));
  cout<<"para 2 "<<fBackground->GetParameter(2)<<endl;
  fBackground->SetParameter(1,fFit->GetParameter(6));
  cout<<"para 3 "<<fBackground->GetParameter(3)<<endl;

  */
  fBackground->SetParameter(0,fExponetial->GetParameter(0));
  fBackground->SetParameter(1,fExponetial->GetParameter(1));

  Double_t maxPureLAr = hPureLAr->GetMaximum();
  cout<<"maximum for Pure LAr is "<<maxPureLAr<<endl;

  for(int i = 0; i < hPureLAr->GetNbinsX() && 0; i ++){
    double binCenter = hPureLAr->GetBinCenter(i+1);
    double val = fBackground->Eval(binCenter);
    double binVal = hPureLAr->GetBinContent(i+1);
    //cout<<"val "<<val<<", binval "<<binVal<<", binCenter "<<binCenter<<", i "<<i<<endl;
    if(binVal-val > 0)
      hPureLAr->SetBinContent(i+1,binVal-val);
    else
      hPureLAr->SetBinContent(i+1,0);
  }

  //1 ppm N2 in Pure LAr
  TFile * f420 = new TFile("N2_1ppm_Doping.root");
  TH1D* hN2LAr = new TH1D("N2LArData","N2LArData",300,10,3000);
  f420->GetObject("TotalChargeAverage",hN2LAr);
  hN2LAr->SetNameTitle("N2LArData","N2LArData");
  hN2LAr->Scale(maxPureLAr/hN2LAr->GetMaximum());

  fFit->SetParameter(0,2.11145e+04);
  fFit->SetParLimits(0,1.911145e+04,2.311145e+04);  

  fFit->SetParameter(1,-4.31004e-02);
  fFit->SetParLimits(1,-2.79754e-01,-2.79754e-03);
  
  fFit->SetParameter(2,6.86702e+03);
  fFit->SetParLimits(2,1.03379e+03,1.03379e+04);
  
  fFit->SetParameter(3,2.77491e+02);
  fFit->SetParLimits(3,200,350);
  
  fFit->SetParameter(4,7.14339e+01);
  fFit->SetParLimits(4,20,80);
  
  fFit->SetParameter(5,8.61522e+02);
  fFit->SetParLimits(5,6.61522e+02,10.61522e+02);
  
  fFit->SetParameter(6,-1.17005e-03);
  fFit->SetParLimits(6,-1.58282e-03,-.88282e-03);

  hN2LAr->Fit("pureFit","","",30,4000);


  fBackground->SetParameter(0,fFit->GetParameter(0));
  cout<<"para 0 "<<fBackground->GetParameter(0)<<endl;;
  fBackground->SetParameter(1,fFit->GetParameter(1));
  cout<<"para 1 "<<fBackground->GetParameter(1)<<endl;
  fBackground->SetParameter(2,fFit->GetParameter(5));
  cout<<"para 2 "<<fBackground->GetParameter(2)<<endl;
  fBackground->SetParameter(3,fFit->GetParameter(6));
  cout<<"para 3 "<<fBackground->GetParameter(3)<<endl;

  for(int i = 0; i < hN2LAr->GetNbinsX(); i ++){
    double binCenter = hN2LAr->GetBinCenter(i+1);
    double val = fBackground->Eval(binCenter);
    double binVal = hN2LAr->GetBinContent(i+1);
    //cout<<"val "<<val<<", binval "<<binVal<<", binCenter "<<binCenter<<", i "<<i<<endl;
    if(binVal-val > 0)
      hN2LAr->SetBinContent(i+1,binVal-val);
    else
      hN2LAr->SetBinContent(i+1,0);
  }


  TFile * outFile = new TFile("ntupleOut.root","recreate");

  TH1D * hXeD[8];
  hXeD[0] = new TH1D("XeD_Before_Smearing","XeD_Before_Smearing",1000,10,10000);
  hXeD[0]->SetLineColor(2);
  hXeD[1] = new TH1D("XeD_After_Smearing1000m_","XeD_After_Smearing1000m",100,10,15000);
  hXeD[2] = new TH1D("XeD_After_Smearing50m_","XeD_After_Smearing50m",100,10,15000);
  hXeD[3] = new TH1D("XeD_After_Smearing20m_","XeD_After_Smearing20m",100,10,15000);
  hXeD[4] = new TH1D("XeD_After_Smearing10m_","XeD_After_Smearing10m",100,10,15000);
  hXeD[5] = new TH1D("XeD_After_Smearing5m_","XeD_After_Smearing5m",300,10,3000);
  hXeD[6] = new TH1D("XeD_After_Smearing3m_","XeD_After_Smearing3m",100,10,15000);
  hXeD[7] = new TH1D("XeD_After_Smearing1m_","XeD_After_Smearing1m",100,10,15000);
  TH1D * hLAr[2];
  hLAr[0] = new TH1D("LAr_Before_Smearing","LAr_Before_Smearing",100,10,10000);
  hLAr[0]->SetLineColor(2);
  hLAr[1] = new TH1D("LAr_After_Smearing_","LAr_After_Smearing",100,10,10000);

  Int_t NEntries = (Int_t)t1->GetEntries();
  for(int i = 0; i < NEntries; i++){
      t1->GetEntry(i);

      if(muonTrigger != 1) continue;

      Double_t nPhotonsSigma = 0.3125;
      Double_t scale = 604./3911;//4235;// 600./4400.;
      Double_t nPhotonsXeD1000m = 0,nPhotonsXeD50m = 0,nPhotonsXeD20m = 0,nPhotonsXeD10m = 0,nPhotonsXeD5m = 0,nPhotonsXeD3m = 0,nPhotonsXeD1m = 0,nPhotons = 0;

      nDetectedOpticalMap = scale*nDetectedOpticalMap;
      nPhotons = ran->Gaus(nDetectedOpticalMap,sqrt(nDetectedOpticalMap) + nDetectedOpticalMap*nPhotonsSigma);
      hLAr[0]->Fill(nDetectedOpticalMap);
      hLAr[1]->Fill(nPhotons);

      XenonDoped1000m = scale*XenonDoped1000m;
      //nPhotonsXeD1000m = ran->Gaus(XenonDoped1000m,sqrt(XenonDoped1000m) + XenonDoped1000m*nPhotonsSigma);
      nPhotonsXeD1000m = ran->Gaus(XenonDoped1000m,XenonDoped1000m*nPhotonsSigma);

      XenonDoped50m = scale*XenonDoped50m;
      //nPhotonsXeD50m = ran->Gaus(XenonDoped10m,sqrt(XenonDoped50m) + XenonDoped50m*nPhotonsSigma);
      nPhotonsXeD50m = ran->Gaus(XenonDoped10m,XenonDoped50m*nPhotonsSigma);

      XenonDoped20m = scale*XenonDoped20m;
      //nPhotonsXeD20m = ran->Gaus(XenonDoped20m,sqrt(XenonDoped20m) + XenonDoped20m*nPhotonsSigma);
      nPhotonsXeD20m = ran->Gaus(XenonDoped20m,XenonDoped20m*nPhotonsSigma);

      XenonDoped10m = scale*XenonDoped10m;
      //nPhotonsXeD10m = ran->Gaus(XenonDoped10m,sqrt(XenonDoped10m) + XenonDoped10m*nPhotonsSigma);
      nPhotonsXeD10m = ran->Gaus(XenonDoped10m,XenonDoped10m*nPhotonsSigma);

      XenonDoped5m = scale*XenonDoped5m;
      //nPhotonsXeD5m = ran->Gaus(XenonDoped5m,sqrt(XenonDoped5m) + XenonDoped5m*nPhotonsSigma);
      nPhotonsXeD5m = ran->Gaus(XenonDoped5m,XenonDoped5m*nPhotonsSigma);

      XenonDoped3m = scale*XenonDoped3m;
      //nPhotonsXeD3m = ran->Gaus(XenonDoped3m,sqrt(XenonDoped3m) + XenonDoped3m*nPhotonsSigma);
      nPhotonsXeD3m = ran->Gaus(XenonDoped3m, XenonDoped3m*nPhotonsSigma);

      XenonDoped1m = scale*XenonDoped1m;
      //nPhotonsXeD1m = ran->Gaus(XenonDoped1m,sqrt(XenonDoped1m) + XenonDoped1m*nPhotonsSigma);
      nPhotonsXeD1m = ran->Gaus(XenonDoped1m, XenonDoped1m*nPhotonsSigma);
      
      hXeD[0]->Fill(XenonDoped);
      hXeD[1]->Fill(nPhotonsXeD1000m);
      hXeD[2]->Fill(nPhotonsXeD50m);
      hXeD[3]->Fill(nPhotonsXeD20m);
      hXeD[4]->Fill(nPhotonsXeD10m);
      hXeD[5]->Fill(nPhotonsXeD5m);
      hXeD[6]->Fill(nPhotonsXeD3m);
      hXeD[7]->Fill(nPhotonsXeD1m);
  }
  TF1 *fLandau = new TF1("LandauFit","[0]*TMath::Landau(x,[1],[2])",10,12000);
  Double_t mean = 0, meanErr = 0,sigma = 0, sigmaErr = 0, maxBin = 0, maxVal = 0;
  Double_t  atten[7] = {1000,50,20,10,5,3,1};
  Double_t ratio[7],ratioErr[7];
  hXeD[0]->Scale(1./hXeD[0]->Integral());
  
  //hXeD[0]->GetXaxis()->SetRangeUser(10,3000);
  hXeD[1]->Scale(1./hXeD[1]->Integral());

  maxVal = hXeD[1]->GetMaximum();
  maxBin = hXeD[1]->GetBinCenter(hXeD[1]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[1]->Fit("LandauFit","","0Q",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"1000m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;
  ratio[0] = sigma/mean;
  ratioErr[0] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));

  /////////////////
  hXeD[2]->Scale(1./hXeD[2]->Integral());

  maxVal = hXeD[2]->GetMaximum();
  maxBin = hXeD[2]->GetBinCenter(hXeD[2]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[2]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"50m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[1] = sigma/mean;
  ratioErr[1] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));


  //////////////////////
  hXeD[3]->Scale(1./hXeD[3]->Integral());

  maxVal = hXeD[3]->GetMaximum();
  maxBin = hXeD[3]->GetBinCenter(hXeD[3]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[3]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"20m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[2] = sigma/mean;
  ratioErr[2] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));

  ////////////////
  hXeD[4]->Scale(1./hXeD[4]->Integral());
  maxVal = hXeD[4]->GetMaximum();
  maxBin = hXeD[4]->GetBinCenter(hXeD[4]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[4]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"10m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[3] = sigma/mean;
  ratioErr[3] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));

  ////////////////


  hXeD[5]->Scale(1./hXeD[5]->Integral());
  maxVal = hXeD[5]->GetMaximum();
  maxBin = hXeD[5]->GetBinCenter(hXeD[5]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[5]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"5m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[4] = sigma/mean;
  ratioErr[4] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));

  ////////////////


  hXeD[6]->Scale(1./hXeD[6]->Integral());
  maxVal = hXeD[6]->GetMaximum();
  maxBin = hXeD[6]->GetBinCenter(hXeD[6]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[6]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"3m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[5] = sigma/mean;
  ratioErr[5] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));


  ////////////////

  hXeD[7]->Scale(1./hXeD[7]->Integral());
  maxVal = hXeD[6]->GetMaximum();
  maxBin = hXeD[6]->GetBinCenter(hXeD[7]->GetMaximumBin());
  fLandau->SetParLimits(0,0,100*maxVal);
  fLandau->SetParameter(0,10*maxVal);
  fLandau->SetParLimits(1,0.5*maxBin,5*maxBin);
  fLandau->SetParameter(1,maxBin);
  fLandau->SetParLimits(2,0.1*maxBin,5*maxBin);
  fLandau->SetParameter(2,0.30*maxBin);
  hXeD[7]->Fit("LandauFit","Q0","",10,12000);

  mean = fLandau->GetParameter(1);
  meanErr = fLandau->GetParError(1);
  sigma = fLandau->GetParameter(2);
  sigmaErr = fLandau->GetParError(2);
  cout<<"1m "<<mean<<"+/-"<<meanErr<<", sigma "<<sigma<<"+/-"<<sigmaErr<<
  ", ratio "<<sigma/mean<<"+/-"<<sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2))<<endl;

  ratio[6] = sigma/mean;
  ratioErr[6] = sqrt(pow(sigmaErr/mean,2)+pow(sigma*meanErr/(mean*mean),2));

  ////////////////
  Double_t  attenErr[7] = {0,0,0,0,0,0,0};
  TGraphErrors * grAtt = new TGraphErrors(7,atten,ratio,attenErr,ratioErr);
  grAtt->Write();

  hLAr[0]->Scale(1./hLAr[0]->Integral());
  //hLAr[0]->GetXaxis()->SetRangeUser(10,3000);
  //hLAr[1]->Scale(1./hLAr[1]->Integral());
  hPureLAr->GetXaxis()->SetRangeUser(50,3000);
  hLAr[1]->Scale(hPureLAr->GetMaximum()/hLAr[1]->GetMaximum());
  //hLAr[1]->GetXaxis()->SetRangeUser(10,3000);
  //f->Close();
  
   TCanvas * c = new TCanvas("c1","c1");
  c->cd();
  //c->SetLogy();
  //hPureLAr->Scale(1./hPureLAr->Integral());
  double max0 = hLAr[1]->GetMaximum();
  double max1 = hPureLAr->GetMaximum();
  hN2LAr->Scale(1./hN2LAr->Integral());
  //hPureLAr->Scale(max0/max1);
  hPureLAr->GetXaxis()->SetRangeUser(10,3000);
  //hPureLAr->GetYaxis()->SetRangeUser(0,max0*1.1);
  hPureLAr->SetLineColor(2);
  hPureLAr->Draw();
  hLAr[1]->Draw("same");
  cout<<"all files done"<<endl;
  hXeDMuons->Write();
  hBackground->Write();
  hSubtraction->Write();
  hPureLAr->Write();
  hN2LAr->Write();
  hXeDN2Muon->Write();

  outFile->Write();
  //outFile->Close();

}
