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
#include <algorithm>
//
#include "TPmtEvent.hxx"
#include "TPmtSimulation.hxx"
#include "TPmtRun.hxx"

void AfterPulse(){

 Float_t bins[] = {1300 ,1400, 1500, 1600,1700 };
 Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1; // or just = 9
 TH1F* hCRatio1PMT1 = new TH1F("CRatio1_1","CRatio1_1", binnum, bins);
 TH1F* hCRatio1PMT2 = new TH1F("CRatio1_2","CRatio1_2", binnum, bins);
 TH1F* hCRatio2PMT1 = new TH1F("CRatio2_1","CRatio2_1", binnum, bins);
 TH1F* hCRatio2PMT2 = new TH1F("CRatio2_2","CRatio2_2", binnum, bins);
 
 TH1F* hTRatio1 = new TH1F("TRatio1","TRatio1", binnum, bins);
 TH1F* hTRatio2 = new TH1F("TRatio2","TRatio2", binnum, bins);
 TH1F* hTRatio3 = new TH1F("TRatio3","TRatio3", binnum, bins);
 TH1F* hTRatio4 = new TH1F("TRatio4","TRatio4", binnum, bins);


  int startNum = 1018;
  int stopNum  = 1021;
  for(int z = startNum; z <= stopNum ; z++){
    TString fileDir = "/home/nmcfadde/RooT/PMT/processedData/DS2/";
    //TString fileName = "SimAnaResults.FillRun_1019.root";
    TString fileName = Form("SimAnaResults.FillRun_%i.root",z);
    cout<<"Opening new data set "<<fileName<<endl;

    TFile * infile = new TFile(fileDir+fileName);

    TTree * processedTree = new TTree();
    infile->GetObject("processedTree",processedTree);

    TPmtEvent* pmtEvent = new TPmtEvent();
    processedTree->SetBranchAddress("processedEvent",&pmtEvent);

    TPmtRun *pmtRun = new TPmtRun();
    processedTree->SetBranchAddress("pmtRun",&pmtRun);

    Int_t NEvents = (Int_t)processedTree->GetEntries();

    TH1D *hStartTimes[2];
    TH1D *hChargeRatio1[2];
    TH1D *hChargeRatio2[2];
    for(int ientry = 0; ientry< NEvents; ientry++){
      processedTree->GetEntry(ientry);
      Int_t nPMTs = pmtRun->nPMTs;

      for(int j = 0; j < nPMTs; j++){ 
        if(ientry == 0){
          hStartTimes[j] = new TH1D(Form("PulseStartTimes%i_%i",j,z),Form("PulseStartTimes%i_%i",j,z),10000,0,10e-6);
          hChargeRatio1[j] = new TH1D(Form("ChargeRatio1_%i_%i",j,z),Form("ChargeRatio1_%i_%i",j,z),1000,0,1);
          hChargeRatio2[j] = new TH1D(Form("ChargeRatio2_%i_%i",j,z),Form("ChargeRatio2_%i_%i",j,z),1000,0,1);
        }
        Int_t nPulses = pmtRun->charge[j].size();
        Double_t c0 = 0,c1 = 0,c2 = 0;
        for(int i = 0; i < nPulses;i++){
          Double_t peakTime   = pmtRun->peakTimes[j][i];
          Double_t charge     = pmtRun->charge[j][i];
          hStartTimes[j]->Fill(peakTime,charge);
          //hStartTimes[j]->Fill(peakTime);
          if(peakTime > 0.6e-6 && peakTime < 1.25e-6) c0 += charge;
          if(peakTime > 1.5e-6 && peakTime < 2e-6) c1 += charge;
          if(peakTime > 2e-6 && peakTime < 2.5e-6) c2 += charge;
        }
        hChargeRatio1[j]->Fill(c1/c0);
        hChargeRatio2[j]->Fill(c2/c0);
      }
    }

    TCanvas *c0 = new TCanvas(Form("Time Distribution%i",z),Form("Time Distribution%i",z));
    c0->cd();
    c0->SetLogy();
    for(int j = 0 ; j < 2; j++){
      Double_t max0 = 0,max1 = 0,max2 = 0,max3 = 0;

      for(int i = 0; i < hStartTimes[j]->GetNbinsX();i++){
        Double_t val = hStartTimes[j]->GetBinContent(i);
        if(i > 900 && i < 1450){
          if(val > max0) max0 = val;
        }
        else if(i > 1450 && i < 2000){
          if(val > max1) max1 = val;
        }
        else if(i > 2000 && i < 2500){
          if(val > max2) max2 = val;
        }
      }
      hTRatio1->SetBinContent(stopNum-z+1,max1/max0);
      hTRatio2->SetBinContent(stopNum-z+1,max2/max0);
      if(j == 0){
      cout<<"pmtNum "<<j<<"\techo 1 R = "<<max1/max0<<", echo 2 R = "<<max2/max0<<endl;
        hStartTimes[j]->Draw();
        hTRatio1->SetBinContent(stopNum-z+1,max1/max0);
        hTRatio2->SetBinContent(stopNum-z+1,max2/max0);
      }
      else{
        cout<<"pmtNum "<<j<<"\techo 1 R = "<<max1/max0<<", echo 2 R = "<<max2/max0<<endl;
        hStartTimes[j]->SetLineColor(2);
        hStartTimes[j]->Draw("same");
        hTRatio3->SetBinContent(stopNum-z+1,max1/max0);
        hTRatio4->SetBinContent(stopNum-z+1,max2/max0);
      }
    }
    TCanvas *c1 = new TCanvas(Form("Time Distribution1_%i",z),Form("Time Distribution1_%i",z));
    c1->cd();
    c1->SetLogy();
    hChargeRatio1[0]->Draw();
    hCRatio1PMT1->SetBinContent(stopNum-z+1,hChargeRatio1[0]->GetMean());
    hChargeRatio1[1]->SetLineColor(2);
    hChargeRatio1[1]->Draw("same");
    hCRatio1PMT2->SetBinContent(stopNum-z+1,hChargeRatio1[1]->GetMean());

    TCanvas *c2 = new TCanvas(Form("Time Distribution2_%i",z),Form("Time Distribution2_%i",z));
    c2->cd();
    c2->SetLogy();
    hChargeRatio2[0]->Draw();
    hCRatio2PMT1->SetBinContent(stopNum-z+1,hChargeRatio2[0]->GetMean());
    hChargeRatio2[1]->SetLineColor(2);
    hCRatio2PMT2->SetBinContent(stopNum-z+1,hChargeRatio2[1]->GetMean());
    hChargeRatio2[1]->Draw("same");



  }

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  hCRatio1PMT1->SetLineColor(3);
  hCRatio1PMT1->Draw();
  hCRatio1PMT2->Draw("same");
  hCRatio1PMT2->SetLineColor(2);
  hCRatio2PMT1->Draw("same");
  hCRatio2PMT1->SetLineColor(3);
  hCRatio2PMT2->Draw("same");
  hCRatio2PMT2->SetLineColor(2);

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->cd();

  hTRatio1->SetLineColor(3);
  hTRatio1->Draw();
  hTRatio2->SetLineColor(2);
  hTRatio2->Draw("same");
  hTRatio3->SetLineColor(4);
  hTRatio3->Draw("same");
  hTRatio4->SetLineColor(6);
  hTRatio4->Draw("same");



}
