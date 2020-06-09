typedef struct{
  // Declaration of leaf types
  Float_t         run;
  Float_t         nev;
  Float_t         mean;
  Float_t         meanErr;
  Float_t         sig;
  Float_t         sigErr;
} ntMembers;

enum {MAXFILE=5};


void ntupRun()
{
  printf(" duh \n");

  TString fname[MAXFILE];
  TFile *fin[MAXFILE];
  TFile *fout = new TFile("runs.root","RECREATE");

  TString binTag("-bins-1250");
  fname[0]=TString(Form("life-10000_10100_DS_2%s.root",binTag.Data()));
  fname[1]=TString(Form("life-20000_20021_DS_2%s.root",binTag.Data()));
  fname[2]=TString(Form("life-20022_20041_DS_2%s.root",binTag.Data()));
  fname[3]=TString(Form("life-20042_20062_DS_2%s.root",binTag.Data()));
  fname[4]=TString(Form("life-20063_29999_DS_2%s.root",binTag.Data()));

  ntMembers val;
  TTree *tr = new TTree("trRun","");
  tr->Branch("trRun",&val,"run:nev:mean:meanErr:sig:sigErr");

  TNtuple *nt;
  Float_t         run;
  Float_t         nev;
  Float_t         mean;
  Float_t         meanErr;
  Float_t         sig;
  Float_t         sigErr;

  vector<float> vrun;
  vector<float> vrunErr;
  vector<float> lan;
  vector<float> lanErr;

  for(int i =0; i<MAXFILE; ++i ) {
    fin[i] = new TFile(fname[i],"READONLY");
    fin[i]->GetObject("ntRun",nt);
    printf(" file %i)  %s  #runs=  %lli  \n",i,fname[i].Data(),nt->GetEntries());
    nt->SetBranchAddress("run",&run);
    nt->SetBranchAddress("nev",&nev);
    nt->SetBranchAddress("mean",&mean);
    nt->SetBranchAddress("meanErr",&meanErr);
    nt->SetBranchAddress("sig",&sig);
    nt->SetBranchAddress("sigErr",&sigErr);
    for(Long64_t j=0; j<nt->GetEntries(); ++j) {
      nt->GetEntry(j);
      float crun=run;
      printf(" \t run %lld number %.0f events %.0f \n",j,run,nev);
      if(crun<20000) crun = run-10100+20000;
      vrun.push_back(crun);
      vrunErr.push_back(0);
      lan.push_back(mean);
      lanErr.push_back(meanErr);
      //
      val.run=run;
      val.nev=nev;
      val.mean=mean;
      val.meanErr=meanErr;
      val.sig=sig;
      val.sigErr=sigErr;
      tr->Fill();
    }
   //printf(" %i %lld \n",i,tr->GetEntries());
  }
  fout->Write();
  
  TGraphErrors *glandau = new TGraphErrors(vrun.size(), &vrun[0],&lan[0],&vrunErr[0],&lanErr[0]);
  TCanvas *clandau = new TCanvas("landau","landau");
  clandau->SetGridx(); clandau->SetGridy();
  glandau->SetTitle("landau");
  glandau->SetMarkerColor(kRed);
  glandau->SetMarkerStyle(22);
  glandau->SetMarkerSize(1);
  glandau->GetXaxis()->SetTitle(" run");
  glandau->GetYaxis()->SetTitle(" landau mean");
  glandau->Draw("ap");
  clandau->Print(".png");

}


