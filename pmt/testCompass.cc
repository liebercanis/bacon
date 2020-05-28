/* MG revised */
#include "testCompass.hh"
#include <TFitResult.h>

compassEvent* comp;

testCompass::testCompass(TString tag, Int_t maxEvents ) 
{
  printf(" starting testCompass tag %s \n",tag.Data());
  TString fileName; fileName.Form("rootData/DS4/compass_run_%s.root",tag.Data());
  printf(" looking for file %s\n",fileName.Data());
  fChain = new TChain("Data");
  fChain->Add(fileName.Data());
  fChain->GetListOfFiles()->ls();
  fChain->SetMakeClass(1);
  fChain->SetBranchStatus("Samples",0);
  comp= new compassEvent(fChain);

  cout << "  fChain with " << fChain->GetEntries() << endl;

  TObjArray  *branches = fChain->GetListOfBranches();


  for (Int_t i=0;i<fChain->GetNbranches();i++){
    TBranch *br = (TBranch*) branches->At(i);
    br->Print();
  }

  comp= new compassEvent(fChain);
  cout << "size " << comp->Samples->GetSize() << endl;

  fChain->SetBranchStatus("Samples",1);

  for(Long64_t entry=0; entry<10; ++entry) {
  comp->getEvent(entry);
  cout << entry << "ch " << comp->Channel << " " <<  comp->Samples->GetSize() << endl;
  }
}
