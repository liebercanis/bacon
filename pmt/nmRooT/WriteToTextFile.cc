#include <iostream>
#include <fstream>
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


using namespace std;

void WriteToTextFile () {

  //TString fName = "XeD_10ppm_N2_10ppm_Doping";
  TString fName = "Pure_LAr";
  TFile *fIn = new TFile(fName+TString(".root"));
  TH1D *h;
  fIn->GetObject("TripletTotalAverage",h);
  
  Int_t nEntries = h->GetEntries();
  ofstream myfile;
  myfile.open (fName+TString(".txt") );
  myfile << fName + TString("\n");
  
  for(int i = 0; i < h->GetNbinsX();i++){
    Double_t binContent = h->GetBinContent(i+1);
    Double_t binCenter = h->GetBinCenter(i+1);
    myfile<<binContent<<" "<<binCenter<<" \n";
  }

  myfile.close();
}
