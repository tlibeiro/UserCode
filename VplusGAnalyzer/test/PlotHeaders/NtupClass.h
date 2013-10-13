//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 11 17:23:21 2013 by ROOT version 5.32/00
// from TTree ntup_VG/ntup_VG
// found on file: /uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput/VplusGTreeAnalyzer_Data.root
//////////////////////////////////////////////////////////

#ifndef NtupClass_h
#define NtupClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class NtupClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         run;
   Float_t         event;
   Float_t         lumi;
   Float_t         HT;
   Float_t         m_V_gr;
   Float_t         pt_V_gr;
   Float_t         pt_V_ungr;
   Float_t         eta_V_gr;
   Float_t         pt_G;
   Float_t         eta_G;
   Float_t         dr_VG;
   Float_t         dphi_VG;
   Float_t         m_VG;
   Float_t         m_V_fft;
   Float_t         pt_V_fft;
   Float_t         eta_V_fft;
   Float_t         maxdrjj_V_fft;
   Float_t         dr_VG_fft;
   Float_t         dphi_VG_fft;
   Float_t         m_VG_fft;
   Float_t         photonIso;
   Float_t         tau2tau1;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_m_V_gr;   //!
   TBranch        *b_pt_V_gr;   //!
   TBranch        *b_pt_V_ungr;   //!
   TBranch        *b_eta_V_gr;   //!
   TBranch        *b_pt_G;   //!
   TBranch        *b_eta_G;   //!
   TBranch        *b_dr_VG;   //!
   TBranch        *b_dphi_VG;   //!
   TBranch        *b_m_VG;   //!
   TBranch        *b_m_V_fft;   //!
   TBranch        *b_pt_V_fft;   //!
   TBranch        *b_eta_V_fft;   //!
   TBranch        *b_maxdrjj_V_fft;   //!
   TBranch        *b_dr_VG_fft;   //!
   TBranch        *b_dphi_VG_fft;   //!
   TBranch        *b_m_VG_fft;   //!
   TBranch        *b_photonIso;   //!
   TBranch        *b_tau2tau1;   //!

   NtupClass(TTree *tree=0);
   virtual ~NtupClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NtupClass_cxx
NtupClass::NtupClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput/VplusGTreeAnalyzer_Data.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput/VplusGTreeAnalyzer_Data.root");
      }
      f->GetObject("ntup_VG",tree);

   }
   Init(tree);
}

NtupClass::~NtupClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NtupClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NtupClass::LoadTree(Long64_t entry)
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

void NtupClass::Init(TTree *tree)
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
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("m_V_gr", &m_V_gr, &b_m_V_gr);
   fChain->SetBranchAddress("pt_V_gr", &pt_V_gr, &b_pt_V_gr);
   fChain->SetBranchAddress("pt_V_ungr", &pt_V_ungr, &b_pt_V_ungr);
   fChain->SetBranchAddress("eta_V_gr", &eta_V_gr, &b_eta_V_gr);
   fChain->SetBranchAddress("pt_G", &pt_G, &b_pt_G);
   fChain->SetBranchAddress("eta_G", &eta_G, &b_eta_G);
   fChain->SetBranchAddress("dr_VG", &dr_VG, &b_dr_VG);
   fChain->SetBranchAddress("dphi_VG", &dphi_VG, &b_dphi_VG);
   fChain->SetBranchAddress("m_VG", &m_VG, &b_m_VG);
   fChain->SetBranchAddress("m_V_fft", &m_V_fft, &b_m_V_fft);
   fChain->SetBranchAddress("pt_V_fft", &pt_V_fft, &b_pt_V_fft);
   fChain->SetBranchAddress("eta_V_fft", &eta_V_fft, &b_eta_V_fft);
   fChain->SetBranchAddress("maxdrjj_V_fft", &maxdrjj_V_fft, &b_maxdrjj_V_fft);
   fChain->SetBranchAddress("dr_VG_fft", &dr_VG_fft, &b_dr_VG_fft);
   fChain->SetBranchAddress("dphi_VG_fft", &dphi_VG_fft, &b_dphi_VG_fft);
   fChain->SetBranchAddress("m_VG_fft", &m_VG_fft, &b_m_VG_fft);
   fChain->SetBranchAddress("photonIso", &photonIso, &b_photonIso);
   fChain->SetBranchAddress("tau2tau1", &tau2tau1, &b_tau2tau1);
   Notify();
}

Bool_t NtupClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NtupClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NtupClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NtupClass_cxx
