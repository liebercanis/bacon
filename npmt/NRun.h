//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 19 15:24:55 2020 by ROOT version 6.12/06
// from TTree ntRun/
// found on file: life-10000_10100_DS_2-All-bins-400.root
//////////////////////////////////////////////////////////

#ifndef NRun_h
#define NRun_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class NRun {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         run;
   Float_t         nev;
   Float_t         mean;
   Float_t         meanErr;
   Float_t         sig;
   Float_t         sigErr;
   Float_t         mpv;
   Float_t         mpvErr;
   Float_t         wid;
   Float_t         widErr;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_nev;   //!
   TBranch        *b_mean;   //!
   TBranch        *b_meanErr;   //!
   TBranch        *b_sig;   //!
   TBranch        *b_sigErr;   //!
   TBranch        *b_mpv;   //!
   TBranch        *b_mpvErr;   //!
   TBranch        *b_wid;   //!
   TBranch        *b_widErr;   //!

   NRun(TTree *tree=0);
   virtual ~NRun();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NRun_cxx
NRun::NRun(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("life-10000_10100_DS_2-All-bins-400.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("life-10000_10100_DS_2-All-bins-400.root");
      }
      f->GetObject("ntRun",tree);

   }
   Init(tree);
}

NRun::~NRun()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NRun::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NRun::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NRun::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("nev", &nev, &b_nev);
   fChain->SetBranchAddress("mean", &mean, &b_mean);
   fChain->SetBranchAddress("meanErr", &meanErr, &b_meanErr);
   fChain->SetBranchAddress("sig", &sig, &b_sig);
   fChain->SetBranchAddress("sigErr", &sigErr, &b_sigErr);
   fChain->SetBranchAddress("mpv", &mpv, &b_mpv);
   fChain->SetBranchAddress("mpvErr", &mpvErr, &b_mpvErr);
   fChain->SetBranchAddress("wid", &wid, &b_wid);
   fChain->SetBranchAddress("widErr", &widErr, &b_widErr);
   Notify();
}

Bool_t NRun::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NRun::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NRun::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NRun_cxx
