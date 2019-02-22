TFile *fin;
TString tag;
void next(Int_t iev=0);

void hlist() {
 
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if(hname.Contains("RawEv")) cout << hname << endl;
  }
}

//
void plot1(Int_t ievent) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("WaveEv%i_PMT_0",ievent);
  TString psearch;
  psearch.Form("PeaksEv%i_PMT_0",ievent);
  printf(" looking for %s %s \n",search.Data(),psearch.Data());

  TH1D* shist=NULL;
  TString ssearch;
  ssearch.Form("SimHitsEv%i_PMT_0",ievent);

  
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search) && !hname.Contains(psearch) ) hist = (TH1D*) obj;;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;;
    if( hname.Contains(ssearch)) shist = (TH1D*) obj;;
  }

  if(!(hist&&phist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(), phist->GetName() );

  phist->SetLineColor(kRed);
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3004);

  hist->SetLineColor(kBlue);
  hist->SetFillColor(kBlue);
  hist->SetFillStyle(0);

  if(shist) {
    shist->SetLineColor(kGreen);
    shist->SetFillColor(kGreen);
    shist->SetFillStyle(0);
  }

  TString canName;
  canName.Form("%s-Event%i",tag.Data(),ievent);

  TCanvas * can = new TCanvas(canName,canName);
  hist->Draw();
  phist->Draw("same");
  if(shist) shist->Draw("same");
  //if(shist) shist->Print("all");
  can->Print(".png");
}

//
void plotb(Int_t ievent) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("RawEv%i_PMT_0",ievent);
  TString psearch;
  psearch.Form("BaseFitEv%i_PMT_0",ievent);
  printf(" looking for %s %s \n",search.Data(),psearch.Data());

  
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search) && !hname.Contains(psearch) ) hist = (TH1D*) obj;;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;;
  }

  if(!(hist&&phist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(), phist->GetName() );

  hist->SetLineColor(kBlue);
  hist->SetFillColor(kBlue);
  hist->SetFillStyle(0);
  
  phist->SetLineColor(kGreen);
  phist->SetFillColor(kGreen);
  phist->SetFillStyle(0);


  TString canName;
  canName.Form("%s-Raw-Event%i",tag.Data(),ievent);

  TCanvas * can = new TCanvas(canName,canName);
  hist->Draw();
  phist->Draw("same");
  //hist->Draw("same");
  can->Print(".png");
}

//
void plotd(Int_t ievent) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;
  TH1D* dhist=NULL;
  TH1D* nhist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("RawEv%i_PMT_0",ievent);
  TString dsearch;
  dsearch.Form("DerEv%i_PMT_0",ievent);
  TString psearch;
  psearch.Form("PeaksEv%i_PMT_0",ievent);

  TString nsearch;
  nsearch.Form("NoiseEv%i_PMT_0",ievent);

  
  TH1D* shist=NULL;
  TString ssearch;
  ssearch.Form("SimHitsEv%i_PMT_0",ievent);

  
  printf(" looking for %s %s %s %s %s \n",search.Data(),dsearch.Data(),psearch.Data(),ssearch.Data(),nsearch.Data());

  
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString hname(obj->GetName());
    if( hname.Contains(search) && !hname.Contains(psearch) ) hist = (TH1D*) obj;;
    if( hname.Contains(psearch)) phist = (TH1D*) obj;
    if( hname.Contains(dsearch)) dhist = (TH1D*) obj;
    if( hname.Contains(ssearch)) shist = (TH1D*) obj;
    if( hname.Contains(nsearch)) nhist = (TH1D*) obj;
  }

  if(!(hist&&phist)) { 
    printf(" cannot find hist for event %i \n",ievent);
    return;
  }

  printf("found %s  %s \n",hist->GetName(), phist->GetName() );

  hist->SetLineColor(kBlue);
  hist->SetFillColor(kBlue);
  hist->SetFillStyle(0);

  dhist->SetLineColor(kBlack);
  dhist->SetFillColor(kBlack);
  dhist->SetFillStyle(0);
  
  phist->SetLineColor(kRed);
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3002);

  if(shist) {
    printf("found %s   \n",shist->GetName() );
    shist->SetLineColor(kGreen);
    shist->SetFillColor(kGreen);
    shist->SetFillStyle(0);
  }


  TString canName;
  canName.Form("%s-Der-Event%i",tag.Data(),ievent);

  TCanvas * can = new TCanvas(canName,canName);
  can->Divide(1,3);
  can->cd(1); hist->Draw();
  can->cd(2); dhist->Draw();
  can->cd(3); {
    phist->Draw();
    //nhist->Draw("same");
    if(shist) shist->Draw("same");
  }
  can->Print(".png");
}


void plotAll(int Max=0) {
 
  TList* list = fin->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list);
  TKey* key;
  TObject* obj;

  int plotted =0;
  while ( (key = (TKey*)next()) ) {
    obj = key->ReadObj() ;
    if (!obj->InheritsFrom("TH1D")) continue;
    TString name(obj->GetName());
    if(name.Contains("Wave")) {
      TString number = name(name.First('E')+2,name.First('P') - name.First('E') -3);
      Int_t iev = number.Atoi();
      cout << name << " number= " << number << endl;
      plot1(iev);
      ++plotted;
      if(Max>0 && plotted == Max) break;
    }
  }

}


void plotEv(TString theTag = "baconRun_run_1015_0_Ev_0_derivative", Int_t ievent=0) 
{
  tag = theTag;
  TString fileName;
  fileName.Form("%s.root",tag.Data());
  printf("looking for file %s\n",fileName.Data()) ;
  fin = TFile::Open(fileName);

  if(!fin) return;

  if (!fin->IsOpen()) {
    printf("<E> Cannot open input file %s\n",fileName.Data()) ;
    exit(1) ;
  }
  printf("opened input file %s\n",fileName.Data()) ;

  plot1(ievent);
}


