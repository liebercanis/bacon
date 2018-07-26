void plots(Int_t ifile = 51,Double_t sigma=3) 
{
  TString fileName ; fileName.Form("baconRunAna_%i_0_%.0fSigma.root",ifile,sigma);

  cout << " looking for file " << fileName << endl;
  TFile *_file0 = TFile::Open(fileName);

  if(_file0->IsZombie()) {
    cout << " not finding file " << fileName << endl;
    return;
  }

  TNtuple *ntEvent;
  _file0->GetObject("ntEvent",ntEvent);

   TNtuple *ntHit;
   _file0->GetObject("ntHit",ntHit);

   TNtuple *ntNHit;
   _file0->GetObject("ntNHit",ntNHit);



   printf(" total number of events is %lli hits %lli nhits %lli \n",ntEvent->GetEntries(),ntHit->GetEntries(),ntNHit->GetEntries());


  TH1F* hq0sum = new TH1F("QSum0"," sum charge pmt 0 ",260,-.01,.5);
  TH1F* hq1sum = new TH1F("QSum1"," sum charge pmt 1 ",260,-.01,.5);
  TH1F* hq0ped = new TH1F("QPed0"," ped pmt 0 ",260,-.01,.5);
  TH1F* hq1ped = new TH1F("QPed1"," ped pmt 1 ",260,-.01,.5);
  
  TH1F* hq0   = new TH1F("allQ0"," all pulses pmt 0 ",320,-.01,.5);
  TH1F* hq1   = new TH1F("allQ1"," all pulses pmt 1 ",320,-.01,.5);

  TH1F* hqc0   = new TH1F("allQC0"," all pulses pmt 0 ",320,-.01,.5);
  TH1F* hqc1   = new TH1F("allQC1"," all pulses pmt 1 ",320,-.01,.5);



  TH1F* hnq0   = new TH1F("allNQ0"," all neg pulses pmt 0 ",320,-.01,.5);
  TH1F* hnq1   = new TH1F("allNQ1"," all neg pulses pmt 1 ",320,-.01,.5);
  hnq0->SetLineColor(kRed); 
  hnq1->SetLineColor(kRed); 

  TH1F* hcq00 = new TH1F("Charge00"," pulse 0 charge pmt 0 ",260,-.01,.5);
  TH1F* hcq01 = new TH1F("Charge01"," pulse 1 charge pmt 0 ",260,-.01,.5);
  TH1F* hcq10 = new TH1F("Charge10"," pulse 0 charge pmt 1 ",260,-.01,.5);
  TH1F* hcq11 = new TH1F("Charge11"," pulse 1 charge pmt 1 ",260,-.01,.5);


  TString all0Name;
  all0Name.Form("qsumPmt0-Run%iSigma%.0f",ifile,sigma);
  TCanvas *call0= new TCanvas(all0Name,all0Name);
  call0->SetLogy();
  ntHit->Draw("q>>allQ0","npmt==0","");
  ntNHit->Draw("q>>allNQ0","npmt==0","sames");
  ntHit->Draw("q>>allQC0","npmt==0&&q>.001","sames");
  call0->Print(".pdf");

  TString all1Name;
  all1Name.Form("qsumPmt1-Run%iSigma%.0f",ifile,sigma);
  TCanvas *call1= new TCanvas(all1Name,all1Name);
  call1->SetLogy();
  ntHit->Draw("q>>allQ1","npmt==1","");
  ntNHit->Draw("q>>allNQ1","npmt==1","sames");
  ntHit->Draw("q>>allQC1","npmt==1&&q>.001","sames");
  call1->Print(".pdf");


  return;


  TString canqsum0Name;
  canqsum0Name.Form("qsumPmt0-Run%i",ifile);
  TCanvas *cqsum0= new TCanvas(canqsum0Name,canqsum0Name);
  cqsum0->SetLogy();
  //ntEvent->Draw("qp0>>QPed0","n0>0","");
  ntEvent->Draw("qsum0>>QSum0","n0>0","");
  cqsum0->Print(".pdf");



  TString canqsum1Name;
  canqsum1Name.Form("qsumPmt1-Run%i",ifile);
  TCanvas *cqsum1= new TCanvas(canqsum1Name,canqsum1Name);
  cqsum1->SetLogy();
  //ntEvent->Draw("qp1>>QPed1","n1>0","");
  ntEvent->Draw("qsum1>>QSum1","n1>0","");
  cqsum1->Print(".pdf");



  
  TString cant0Name;
  cant0Name.Form("startpmt0-Run%i",ifile);
  TCanvas *ct0 = new TCanvas(cant0Name,cant0Name);
  ntEvent->Draw("t00","n0>0","");
  ntEvent->Draw("t01","n0>1","sames");
  ct0->Print(".pdf");


  TString cant1Name;
  cant1Name.Form("startpmt1-Run%i",ifile);
  TCanvas *ct1 = new TCanvas(cant1Name,cant1Name);
  ntEvent->Draw("t10","n1>0","");
  ntEvent->Draw("t11","n1>1","sames");
  ct1->Print(".pdf");




  TString canq0Name;
  canq0Name.Form("qpmt0-Run%i",ifile);
  TCanvas *cq0 = new TCanvas(canq0Name,canq0Name);
  cq0->SetLogy();
  hq0sum->Draw();
  ntEvent->Draw("q00>>Charge00","n0>0","sames");
  ntEvent->Draw("q01>>Charge01","n0>1","sames");
  //TPaveStats *st = (TPaveStats*) hcq01->GetListOfFunctions()->FindObject("stats");
  //if(st) st->SetX1NDC(0); //new x start position
  //if(st) st->SetX2NDC(0); //new x end position
  cq0->Print(".pdf");


  TString canq1Name;
  canq1Name.Form("qpmt1-Run%i",ifile);
  TCanvas *cq1= new TCanvas(canq1Name,canq1Name);
  cq1->SetLogy();
  ntEvent->Draw("q10>>Charge10","n1>0","");
  ntEvent->Draw("q11>>Charge11","n1>1","sames");
  cq1->Print(".pdf");

  
  TString cq01OneName;
  cq01OneName.Form("qpm0to1-One-Run%i",ifile);
  TCanvas *cq01One = new TCanvas(cq01OneName,cq01OneName);
  cq01One->SetLogz();
  ntEvent->Draw("q00:q10","n0>0&&n1>0","contz");
  cq01One->Print(".png");

  TString cq01TwoName;
  cq01TwoName.Form("qpm0to1-Two-Run%i",ifile);
  TCanvas *cq01Two = new TCanvas(cq01TwoName,cq01TwoName);
  cq01Two->SetLogz();
  ntEvent->Draw("q01:q11","n0>1&&n1>1","contz");
  cq01Two->Print(".png");

}
