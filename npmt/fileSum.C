
void fileSum()
{
  //double lifeBins = 1.0E4/8.0; // ns bins

  TString dirName("/data1/bacon/processedData/DS2");
  TSystemDirectory dir("dataFiles",dirName);
  TList *files = dir.GetListOfFiles();

  TIter next(files);
  TSystemFile *file;
  vector<string> fname;
  vector<string> ftag;
  while( (file = (TSystemFile*) next()) ) {
    string fullName = string(file->GetName());
    string exten  = fullName.substr( fullName.find_last_of(".")+1 );
    if(exten!=string("root")) continue;
    string name = fullName.substr( 0, fullName.find_last_of(".") );
    cout << fullName << endl;
    fname.push_back(fullName);
    ftag.push_back(name);
  }

  return;

  TFile* fout = new TFile("histData.root","recreate");

  cout << " found data " << fname.size() << " files " << endl;
  for(unsigned ifile=0 ; ifile < fname.size(); ++ifile) {
    cout << " full " << fname[ifile] << " tag " << ftag[ifile] << endl;
    //TH1D* hLife = (TH1D*) hist->Clone(Form("Life-%s",ftag[ifile].c_str()));
    //hLife->SetTitle( Form("Lifetime hist %s",ftag[ifile].c_str()) ) ;
  }
  fout->Write();
  fout->ls();
}
