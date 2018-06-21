{
  printf("\n this is rootlogon for bacon \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  gSystem->AddIncludePath(" -I. -I./obj/");
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  gROOT->LoadMacro("util.C");
  //gROOT->LoadMacro("ntupleRun.C",&iload);
  int iload = gSystem->Load("$HOME/RnD/bacon/obj/libBaconRoot.so");
  printf(" loaded libPmtRoot = %i zero is success! \n",iload);
}

