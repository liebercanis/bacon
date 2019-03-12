{
  printf("\n this is rootlogon for bacon \n");
  TString arch=gSystem->GetBuildArch();
  cout << " arch is " << arch << endl; 
  gSystem->AddIncludePath(" -I. -I./obj/");
  printf(" include path %s \n\n",gSystem->GetIncludePath());
  gROOT->LoadMacro("util.C");
  //int iload = gSystem->Load("/home/admin/MGDO/tam/lib/libTAM.so");
  //printf(" loaded libTAM.so = %i zero is success! \n",iload);
}

