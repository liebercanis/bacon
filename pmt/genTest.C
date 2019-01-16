 void genTest(Int_t max=10000) {

   TRandom *rand = new TRandom2();

  Double_t startTime = 0;
  Double_t stopTime = 4e-6;
  TH1D* hTest = new TH1D("test","test",1000,startTime,stopTime);

  Double_t tau = 1.0E-6;

  Int_t counter = 0,failCounter = 0;

  for(Int_t counter=0; counter< max; ++counter)  {
    Double_t x = -tau*TMath::Log(1-rand->Rndm());
    hTest->Fill(x);
  }

  TCanvas* can = new TCanvas("genTest","genTest");
  hTest->Draw();

}


