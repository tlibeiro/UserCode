//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 24 14:33:54 2013 by ROOT version 5.32/00
// from TTree verifylhe/verifylhe
// found on file: /uscms_data/d3/tlibeiro/GenerateWGamma/CMSSW_5_3_8_patch1/src/GeneratorInterface/LHEInterface/test/verifyLHE.root
//////////////////////////////////////////////////////////

#ifndef VerifyLHEEventClass_h
#define VerifyLHEEventClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;
// Header file for the classes stored in the TTree if any.
// Fixed size dimensions of array or collections stored in the TTree if any.

class VerifyLHEEventClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain         *chain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         Wmass;
   Float_t         Wpt;
   Float_t         Weta;
   Float_t         Wphi;
   Float_t         photonpt;
   Float_t         photoneta;
   Float_t         photonphi;
   Float_t         photonWdr;
   Float_t         photonWdphi;
   Float_t         wdau1pt;
   Float_t         wdau1eta;
   Float_t         wdau1phi;
   Float_t         wdau2pt;
   Float_t         wdau2eta;
   Float_t         wdau2phi;
   Float_t         wdausdr;
   Float_t         num_W_daus;
   Float_t         num_radjets;
   Float_t         radjet1pt;
   Float_t         radjet2pt;
   Float_t         radjet3pt;
   Float_t         radjet4pt;

   // List of branches
   TBranch        *b_Wmass;   //!
   TBranch        *b_Wpt;   //!
   TBranch        *b_Weta;   //!
   TBranch        *b_Wphi;   //!
   TBranch        *b_photonpt;   //!
   TBranch        *b_photoneta;   //!
   TBranch        *b_photonphi;   //!
   TBranch        *b_photonWdr;   //!
   TBranch        *b_photonWdphi;   //!
   TBranch        *b_wdau1pt;   //!
   TBranch        *b_wdau1eta;   //!
   TBranch        *b_wdau1phi;   //!
   TBranch        *b_wdau2pt;   //!
   TBranch        *b_wdau2eta;   //!
   TBranch        *b_wdau2phi;   //!
   TBranch        *b_wdausdr;   //!
   TBranch        *b_num_W_daus;   //!
   TBranch        *b_num_radjets;   //!
   TBranch        *b_radjet1pt;   //!
   TBranch        *b_radjet2pt;   //!
   TBranch        *b_radjet3pt;   //!
   TBranch        *b_radjet4pt;   //!

   VerifyLHEEventClass(const vector<string>& infile, const string& outfile,TTree *tree=0);
   virtual ~VerifyLHEEventClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   TFile* f_out;
   TNtuple* ntup_VG;
   string outfile;
};

#endif

#ifdef VerifyLHEEventClass_cxx
VerifyLHEEventClass::VerifyLHEEventClass(const vector<string>& infiles, const string& ofile,TTree *tree) : fChain(0), outfile(ofile) 
{

	TChain * chain = new TChain("verifylheevents/verifylhe","");
	const unsigned numinput(infiles.size());
	for(unsigned i(0);i<numinput;++i)
		chain->Add(infiles[i].c_str());
	tree = chain;
	Init(tree);


 ostringstream ntup_vars; 
	vector<string> vars; 
  vars.push_back("m_V_gr");
  vars.push_back("pt_V_gr");
  vars.push_back("pt_V_ungr");
  vars.push_back("eta_V_gr");
  vars.push_back("pt_G");
	vars.push_back("eta_G");
	vars.push_back("dr_VG");
	vars.push_back("dphi_VG");
	vars.push_back("m_VG");
	const unsigned numvars(9);
	for(unsigned i(0);i<numvars;++i)
		if(!i) ntup_vars<<vars[i];
		else ntup_vars<<':'<<vars[i];

	ntup_VG = new TNtuple("ntup_VG","ntup_VG",ntup_vars.str().c_str());

}

VerifyLHEEventClass::~VerifyLHEEventClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t VerifyLHEEventClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t VerifyLHEEventClass::LoadTree(Long64_t entry)
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

void VerifyLHEEventClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("Wmass", &Wmass, &b_Wmass);
   fChain->SetBranchAddress("Wpt", &Wpt, &b_Wpt);
   fChain->SetBranchAddress("Weta", &Weta, &b_Weta);
   fChain->SetBranchAddress("Wphi", &Wphi, &b_Wphi);
   fChain->SetBranchAddress("photonpt", &photonpt, &b_photonpt);
   fChain->SetBranchAddress("photoneta", &photoneta, &b_photoneta);
   fChain->SetBranchAddress("photonphi", &photonphi, &b_photonphi);
   fChain->SetBranchAddress("photonWdr", &photonWdr, &b_photonWdr);
   fChain->SetBranchAddress("photonWdphi", &photonWdphi, &b_photonWdphi);
   fChain->SetBranchAddress("wdau1pt", &wdau1pt, &b_wdau1pt);
   fChain->SetBranchAddress("wdau1eta", &wdau1eta, &b_wdau1eta);
   fChain->SetBranchAddress("wdau1phi", &wdau1phi, &b_wdau1phi);
   fChain->SetBranchAddress("wdau2pt", &wdau2pt, &b_wdau2pt);
   fChain->SetBranchAddress("wdau2eta", &wdau2eta, &b_wdau2eta);
   fChain->SetBranchAddress("wdau2phi", &wdau2phi, &b_wdau2phi);
   fChain->SetBranchAddress("wdausdr", &wdausdr, &b_wdausdr);
   fChain->SetBranchAddress("num_W_daus", &num_W_daus, &b_num_W_daus);
   fChain->SetBranchAddress("num_radjets", &num_radjets, &b_num_radjets);
   fChain->SetBranchAddress("radjet1pt", &radjet1pt, &b_radjet1pt);
   fChain->SetBranchAddress("radjet2pt", &radjet2pt, &b_radjet2pt);
   fChain->SetBranchAddress("radjet3pt", &radjet3pt, &b_radjet3pt);
   fChain->SetBranchAddress("radjet4pt", &radjet4pt, &b_radjet4pt);
   Notify();
}

Bool_t VerifyLHEEventClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void VerifyLHEEventClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t VerifyLHEEventClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef VerifyLHEEventClass_cxx
