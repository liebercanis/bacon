TFile *fin;

void next(Int_t ifile=51,Int_t iev=0);

void plot1(Int_t ievent=0, Int_t ifile = 51,Int_t isigma=3) 
{
  TString fileName;
  fileName.Form("baconRunAna_%i_0_%iSigma.root",ifile,isigma);

  
  fin = TFile::Open(fileName);

  if(!fin) return;

  printf("looking for file %s\n",fileName.Data()) ;
  if (!fin->IsOpen()) {
    printf("<E> Cannot open input file %s\n",fileName.Data()) ;
    exit(1) ;
  }
  printf("opened input file %s\n",fileName.Data()) ;

  next(ifile,ievent);
}

void next(Int_t ifile,Int_t ievent) {

  TH1D* hist=NULL;
  TH1D* phist=NULL;

  if (!fin->IsOpen()) return;

  TString search;
  search.Form("Ev_%i",ievent);
  TString psearch;
  psearch.Form("PeaksEv_%i",ievent);
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

  phist->SetLineColor(kRed);
  phist->SetFillColor(kRed);
  phist->SetFillStyle(3002);

  TString canName;
  canName.Form("Run%i-Event%i",ifile,ievent);

  TCanvas * can = new TCanvas(canName,canName);
  hist->Draw();
  phist->Draw("same");
  can->Print(".pdf");
}
