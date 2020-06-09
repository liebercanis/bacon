
std::vector<double> vx;
std::vector<double> vy;


void read1(string fileName, TH1D* hist) 
{
  ifstream in;
  string full = string("dataFiles/")+fileName;
  in.open(full);
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.c_str());
    return;
  }

  cout << " opened " << fileName << endl;
  double t,v;
  char name[120];
  int nlines =0;
  while (in.good()) {
    if(nlines==0) in >> name;
    else in >> v >> t;
    //printf(" %E %E \n",t,v);
    t=t*1E6;
    int hitBin =  hist->FindBin(t); 
    hist->SetBinContent(hitBin, hist->GetBinContent(hitBin)+v);
    ++nlines;
  }
  in.close();
  cout << " read " << name << " lines  " << nlines << " entries  " << hist->GetEntries()  <<  endl;
}

TGraph *gread1(string fileName) 
{
  TGraph *gr = NULL;
  ifstream in;
  string full = string("dataFiles/")+fileName;
  in.open(full);
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.c_str());
    return gr;
  }

  vx.clear();
  vy.clear();
  cout << " opened " << fileName << endl;
  double t,v;
  char name[120];
  int nlines =0;
  while (in.good()) {
    if(nlines==0) in >> name;
    else in >> v >> t;
    vx.push_back(t);
    vy.push_back(v);
    ++nlines;
  }
  in.close();
  gr = new TGraph(vx.size(), &vx[0],&vy[0]);
  cout << " read " << name << " lines  " << nlines << " entries  " << gr->GetN()  <<  endl;
  //for(unsigned i=0; i<10; ++i) cout << vx[i+1] - vx[i] << endl;
  return gr;
}




void readData()
{
  //double lifeBins = 1.0E4/8.0; // ns bins
  //double maxLife=10.0;
  double lifeBins = 400;
  double maxLife=10.0;

  TH1D* hist = new TH1D("Life"," lifetime ",lifeBins,0,maxLife);
  hist->GetXaxis()->SetTitle(" micro-seconds ");
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(22);
  hist->SetMarkerSize(.2);

  TSystemDirectory dir("dataFiles","dataFiles");
  TList *files = dir.GetListOfFiles();

  TIter next(files);
  TSystemFile *file;
  vector<string> fname;
  vector<string> ftag;
  while( (file = (TSystemFile*) next()) ) {
    string fullName = string(file->GetName());
    string exten  = fullName.substr( fullName.find_last_of(".")+1 );
    if(exten!=string("txt")) continue;
    string name = fullName.substr( 0, fullName.find_last_of(".") );
    fname.push_back(fullName);
    ftag.push_back(name);
  }

  TFile* fout = new TFile("histData.root","recreate");

  cout << " found data " << fname.size() << " files " << endl;
  for(unsigned ifile=0 ; ifile < fname.size(); ++ifile) {
    cout << " full " << fname[ifile] << " tag " << ftag[ifile] << endl;
    TH1D* hLife = (TH1D*) hist->Clone(Form("Life-%s",ftag[ifile].c_str()));
    hLife->SetTitle( Form("Lifetime hist %s",ftag[ifile].c_str()) ) ;
    read1(fname[ifile],hLife);
    TGraph *gr = gread1(fname[ifile]);
    if(!gr) continue;
    gr->SetName( Form("gLife-%s",ftag[ifile].c_str()) ) ;
    gr->SetTitle( Form("Lifetime graph %s",ftag[ifile].c_str()) ) ;
    fout->Append(gr);
  }
  fout->Write();
  fout->ls();
    
}
