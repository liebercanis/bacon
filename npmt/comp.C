TH1D* histNorm(TH1D* hist, double nev, double meanSpe=1, bool setErr=true)
{

  TH1D* hnew = (TH1D*) hist->Clone(Form("%s%s",hist->GetName(),"Norm"));
  for(Int_t ibin=0; ibin < hist->GetNbinsX() ; ++ibin) {
    hnew->SetBinContent(ibin, hist->GetBinContent(ibin)/nev*meanSpe);
    if(setErr) hnew->SetBinError(ibin,hist->GetBinError(ibin)/nev*meanSpe);
    else hnew->SetBinError(ibin,sqrt(hist->GetBinContent(ibin)));
  }
  return hnew;
}



void comp() {
  enum {MAXFILE=5};
  enum {PFITS=5};
  TF1 *fp[MAXFILE];
  TH1D* hLife[MAXFILE];
  TH1D* hLifeNorm[MAXFILE];

  TFile *cfile = new TFile("comps.root","RECREATE");

  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    hLife[ifile]=NULL;
    fp[ifile]=NULL;
  }
  TString tag("DS_2");

  TString runTag[MAXFILE];
  runTag[0]=TString("00_PPM");
  runTag[1]=TString("01_PPM");
  runTag[2]=TString("02_PPM");
  runTag[3]=TString("05_PPM");
  runTag[4]=TString("10_PPM");

  TString runRange[MAXFILE];
  runRange[0]=TString("10000_10100");
  runRange[1]=TString("20000_20021");
  runRange[2]=TString("20022_20041");
  runRange[3]=TString("20042_20062");
  runRange[4]=TString("20063_29999");

  
  double runPPM[MAXFILE]={0,1,2,5,10};
  
  TFile *fin[MAXFILE];

  TString fname[MAXFILE];
  TString trigTag("All");
  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    fname[ifile]=TString(Form("life-%s_DS_2-%s-bins-400.root",runRange[ifile].Data(),trigTag.Data()));
    cout << ifile << ")  " << "  " <<fname[ifile] << endl;
  }

  TH1D* hist;
  TH1D* hist2;
  TH1D* hist3;
  TH1D* hCharge;
  vector<double> vEvents;
  vector<double> vSpe;
  TNtuple *ntRun;
  for(int ifile =0; ifile<MAXFILE; ++ifile ) {
    fin[ifile] = new TFile(fname[ifile],"READONLY");
    fin[ifile]->GetObject("LifeAll",hist);
    fin[ifile]->GetObject("EventSumCut",hist2);
    fin[ifile]->GetObject("EventSumNeil",hist3);
    fin[ifile]->GetObject("ntRun",ntRun);

    // runs for this file
    Float_t frun,fnev,fmean,fmeanErr,fsig,fsigErr,fmpv,fmpvErr,fwid,fwidErr;
    ntRun->SetBranchAddress("run",&frun);
    ntRun->SetBranchAddress("nev",&fnev);
    ntRun->SetBranchAddress("mean",&fmean);
    ntRun->SetBranchAddress("meanErr",&fmeanErr);
    ntRun->SetBranchAddress("sig",&fsig);
    ntRun->SetBranchAddress("sigErr",&fsigErr);
    ntRun->SetBranchAddress("mpv",&fmpv);
    ntRun->SetBranchAddress("mpvErr",&fmpvErr);
    ntRun->SetBranchAddress("wid",&fwid);
    ntRun->SetBranchAddress("widErr",&fwidErr);
    ntRun->GetEntry(0);
    int runFirst, runLast;
    runFirst = frun;
    ntRun->GetEntry( ntRun->GetEntries()-1);
    runLast = frun;

    if(!hist) continue;
    //cout << "file " << fin[ifile]->GetName() << "  " << hist->GetName() << " " << hist2->GetName() << " " << hist3->GetName() <<endl;
    hLife[ifile]= (TH1D*) hist->Clone(Form("LifeAll%i",ifile));
    hCharge= (TH1D*) hist3->Clone(Form("ChargeAll%i",ifile));
    TObject * fobj= hLife[ifile]->GetListOfFunctions()->FindObject("modelFit0");
    //hLife[ifile]->GetListOfFunctions()->ls();
    if(fobj) cout << " remove  fobj " << fobj->GetName() << endl;
    hLife[ifile]->GetListOfFunctions()->RecursiveRemove(fobj);
    hLife[ifile]->SetTitle(Form("life-%s-pass-%s",runTag[ifile].Data(),trigTag.Data()));
    double numEvents =  hCharge->GetEntries();
    vEvents.push_back(numEvents);
    double numPulses = hLife[ifile]->GetEntries();
    double sumSPE = hLife[ifile]->Integral();
    hLifeNorm[ifile] = histNorm( hLife[ifile],numEvents);
    //hLife[ifile]->SetNormFactor(numEvents);
    printf(" file %3i has %3lld runs from %i to %i events %0.f binwidth %f pulses %f SPE/event %f  \n",
        ifile,ntRun->GetEntries(), runFirst, runLast, numEvents, hLife[ifile]->GetBinWidth(1), numPulses, sumSPE/numEvents );
    vSpe.push_back(sumSPE);
    cfile->Append(hLife[ifile]);
    cfile->Append(hLifeNorm[ifile]);
    cfile->Append(hCharge);
  }
  for(int ifile=0; ifile<MAXFILE; ++ifile) if(!hLife[ifile]) return;

  TFile *file1 = new TFile("histData.root","readonly");

  //     Life-XeD_10ppm_N2_10ppm_Doping;
  TH1D* hNLife[MAXFILE];
  TString Nnames[MAXFILE];
  TString Mnames[MAXFILE];
  Nnames[0]=TString("Life-Pure_LAr");
  Nnames[1]=TString("Life-XeD_1ppm_Doping");
  Nnames[2]=TString("Life-XeD_2ppm_Doping");
  Nnames[3]=TString("Life-XeD_5ppm_Doping");
  Nnames[4]=TString("Life-XeD_10ppm_Doping");
  Mnames[0]=TString("NLifePure");
  Mnames[1]=TString("NLife1ppm");
  Mnames[2]=TString("NLife2ppm");
  Mnames[3]=TString("NLife5ppm");
  Mnames[4]=TString("NLife10ppm");

  TH1D *hista;

  for(int ifile =0; ifile<MAXFILE; ++ifile )  {
    file1->GetObject(Nnames[ifile].Data(),hista);
    double total =  hista->Integral();
    printf(" file %i integral %f bw %f \n",ifile,total,hista->GetBinWidth(1));
    hista->SetName(Mnames[ifile]);
    hNLife[ifile] = histNorm(hista,1.0,vSpe[ifile],false);
    cfile->Append(hNLife[ifile]);
  }

  //Life-N2_1ppm_Doping
  //Life-XeD_10ppm_N2_1ppm_Doping

  cfile->Write();
  cfile->ls();

  TCanvas *can[MAXFILE];
  gStyle->SetOptStat(11111111);

  for(int ifile=0; ifile<MAXFILE; ++ifile) {
    can[ifile] = new TCanvas(Form("comp-%s",runTag[ifile].Data()),Form("comp-%s",runTag[ifile].Data()));
    hLife[ifile]->SetTitle( Form("compLife-%s",runTag[ifile].Data()) );
    hNLife[ifile]->SetMarkerColor(kRed);
    hNLife[ifile]->SetLineColor(kRed);
    hNLife[ifile]->SetMarkerStyle(22);
    hNLife[ifile]->SetMarkerSize(.4);
    hLife[ifile]->SetMarkerColor(kBlack);
    hLife[ifile]->SetLineColor(kBlack);
    hLife[ifile]->SetMarkerStyle(21);
    hLife[ifile]->SetMarkerSize(.4);
    can[ifile]->SetLogy();
    hLife[ifile]->SetStats(11111);
    hLife[ifile]->Draw("E");
    hNLife[ifile]->Draw("Psame");
    can[ifile]->Print(".png");
  }


}
