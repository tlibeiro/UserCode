//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  9 17:17:16 2013 by ROOT version 5.32/00
// from TTree VplusGTree/VplusGTree
// found on file: /uscms_data/d3/tlibeiro/VPlusGamma/condorfiles/fftjetTest/VplusGTree_1_fftjetTest.root
//////////////////////////////////////////////////////////

#ifndef tmpClass_h
#define tmpClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class tmpClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event_runNo;
   Int_t           event_evtNo;
   Int_t           event_lumi;
   Int_t           event_bunch;
   Float_t         event_HT;
   Float_t         GenPart_pt[10];
   Float_t         GenPart_eta[10];
   Float_t         GenPart_phi[10];
   Float_t         GenPart_m[10];
   Float_t         GenPart_id[10];
   Float_t         GenPart_id_mother[10];
   Float_t         FFTJetLA2R20_pt[6];
   Float_t         FFTJetLA2R20_eta[6];
   Float_t         FFTJetLA2R20_phi[6];
   Float_t         FFTJetLA2R20_m[6];
   Float_t         FFTJetLA2R20_splft[6];
   Float_t         FFTJetLA2R20_mglft[6];
   Float_t         FFTJetLA2R20_scale[6];
   Float_t         FFTJetLA3R20_pt[6];
   Float_t         FFTJetLA3R20_eta[6];
   Float_t         FFTJetLA3R20_phi[6];
   Float_t         FFTJetLA3R20_m[6];
   Float_t         FFTJetLA3R20_splft[6];
   Float_t         FFTJetLA3R20_mglft[6];
   Float_t         FFTJetLA3R20_scale[6];
   Float_t         FFTJetLA4R20_pt[6];
   Float_t         FFTJetLA4R20_eta[6];
   Float_t         FFTJetLA4R20_phi[6];
   Float_t         FFTJetLA4R20_m[6];
   Float_t         FFTJetLA4R20_splft[6];
   Float_t         FFTJetLA4R20_mglft[6];
   Float_t         FFTJetLA4R20_scale[6];
   Float_t         FFTJetLA2R30_pt[6];
   Float_t         FFTJetLA2R30_eta[6];
   Float_t         FFTJetLA2R30_phi[6];
   Float_t         FFTJetLA2R30_m[6];
   Float_t         FFTJetLA2R30_splft[6];
   Float_t         FFTJetLA2R30_mglft[6];
   Float_t         FFTJetLA2R30_scale[6];
   Float_t         FFTJetLA3R30_pt[6];
   Float_t         FFTJetLA3R30_eta[6];
   Float_t         FFTJetLA3R30_phi[6];
   Float_t         FFTJetLA3R30_m[6];
   Float_t         FFTJetLA3R30_splft[6];
   Float_t         FFTJetLA3R30_mglft[6];
   Float_t         FFTJetLA3R30_scale[6];
   Float_t         FFTJetLA4R30_pt[6];
   Float_t         FFTJetLA4R30_eta[6];
   Float_t         FFTJetLA4R30_phi[6];
   Float_t         FFTJetLA4R30_m[6];
   Float_t         FFTJetLA4R30_splft[6];
   Float_t         FFTJetLA4R30_mglft[6];
   Float_t         FFTJetLA4R30_scale[6];
   Float_t         FFTJetLA2R40_pt[6];
   Float_t         FFTJetLA2R40_eta[6];
   Float_t         FFTJetLA2R40_phi[6];
   Float_t         FFTJetLA2R40_m[6];
   Float_t         FFTJetLA2R40_splft[6];
   Float_t         FFTJetLA2R40_mglft[6];
   Float_t         FFTJetLA2R40_scale[6];
   Float_t         FFTJetLA3R40_pt[6];
   Float_t         FFTJetLA3R40_eta[6];
   Float_t         FFTJetLA3R40_phi[6];
   Float_t         FFTJetLA3R40_m[6];
   Float_t         FFTJetLA3R40_splft[6];
   Float_t         FFTJetLA3R40_mglft[6];
   Float_t         FFTJetLA3R40_scale[6];
   Float_t         FFTJetLA4R40_pt[6];
   Float_t         FFTJetLA4R40_eta[6];
   Float_t         FFTJetLA4R40_phi[6];
   Float_t         FFTJetLA4R40_m[6];
   Float_t         FFTJetLA4R40_splft[6];
   Float_t         FFTJetLA4R40_mglft[6];
   Float_t         FFTJetLA4R40_scale[6];
   Float_t         FFTJetLA2R50_pt[6];
   Float_t         FFTJetLA2R50_eta[6];
   Float_t         FFTJetLA2R50_phi[6];
   Float_t         FFTJetLA2R50_m[6];
   Float_t         FFTJetLA2R50_splft[6];
   Float_t         FFTJetLA2R50_mglft[6];
   Float_t         FFTJetLA2R50_scale[6];
   Float_t         FFTJetLA3R50_pt[6];
   Float_t         FFTJetLA3R50_eta[6];
   Float_t         FFTJetLA3R50_phi[6];
   Float_t         FFTJetLA3R50_m[6];
   Float_t         FFTJetLA3R50_splft[6];
   Float_t         FFTJetLA3R50_mglft[6];
   Float_t         FFTJetLA3R50_scale[6];
   Float_t         FFTJetLA4R50_pt[6];
   Float_t         FFTJetLA4R50_eta[6];
   Float_t         FFTJetLA4R50_phi[6];
   Float_t         FFTJetLA4R50_m[6];
   Float_t         FFTJetLA4R50_splft[6];
   Float_t         FFTJetLA4R50_mglft[6];
   Float_t         FFTJetLA4R50_scale[6];
   Int_t           NumPhotons;
   Float_t         Photon_pt[5];   //[NumPhotons]
   Float_t         Photon_e[5];   //[NumPhotons]
   Float_t         Photon_eta[5];   //[NumPhotons]
   Float_t         Photon_phi[5];   //[NumPhotons]
   Float_t         Photon_Theta[5];   //[NumPhotons]
   Float_t         Photon_Vx[5];   //[NumPhotons]
   Float_t         Photon_Vy[5];   //[NumPhotons]
   Float_t         Photon_Vz[5];   //[NumPhotons]
   Float_t         Photon_SC_Et[5];   //[NumPhotons]
   Float_t         Photon_SC_E[5];   //[NumPhotons]
   Float_t         Photon_SC_Eta[5];   //[NumPhotons]
   Float_t         Photon_SC_Phi[5];   //[NumPhotons]
   Float_t         Photon_SC_Theta[5];   //[NumPhotons]
   Float_t         Photon_SC_x[5];   //[NumPhotons]
   Float_t         Photon_SC_y[5];   //[NumPhotons]
   Float_t         Photon_SC_z[5];   //[NumPhotons]
   Float_t         PFisocharged03[5];   //[NumPhotons]
   Float_t         PFisophoton03[5];   //[NumPhotons]
   Float_t         PFisoneutral03[5];   //[NumPhotons]
   Float_t         Photon_HoverE[5];   //[NumPhotons]
   Float_t         Photon_HoverE2011[5];   //[NumPhotons]
   Float_t         Photon_SigmaIetaIeta[5];   //[NumPhotons]
   Int_t           Photon_hasPixelSeed[5];   //[NumPhotons]
   Int_t           Photon_passElecVeto[5];   //[NumPhotons]
   Int_t           Photon_Id2011[5];   //[NumPhotons]
   Int_t           Photon_Id2012[5];   //[NumPhotons]

   // List of branches
   TBranch        *b_event_runNo;   //!
   TBranch        *b_event_evtNo;   //!
   TBranch        *b_event_lumi;   //!
   TBranch        *b_event_bunch;   //!
   TBranch        *b_event_HT;   //!
   TBranch        *b_GenPart_pt;   //!
   TBranch        *b_GenPart_eta;   //!
   TBranch        *b_GenPart_phi;   //!
   TBranch        *b_GenPart_m;   //!
   TBranch        *b_GenPart_id;   //!
   TBranch        *b_GenPart_id_mother;   //!
   TBranch        *b_FFTJetLA2R20_pt;   //!
   TBranch        *b_FFTJetLA2R20_eta;   //!
   TBranch        *b_FFTJetLA2R20_phi;   //!
   TBranch        *b_FFTJetLA2R20_m;   //!
   TBranch        *b_FFTJetLA2R20_splft;   //!
   TBranch        *b_FFTJetLA2R20_mglft;   //!
   TBranch        *b_FFTJetLA2R20_scale;   //!
   TBranch        *b_FFTJetLA3R20_pt;   //!
   TBranch        *b_FFTJetLA3R20_eta;   //!
   TBranch        *b_FFTJetLA3R20_phi;   //!
   TBranch        *b_FFTJetLA3R20_m;   //!
   TBranch        *b_FFTJetLA3R20_splft;   //!
   TBranch        *b_FFTJetLA3R20_mglft;   //!
   TBranch        *b_FFTJetLA3R20_scale;   //!
   TBranch        *b_FFTJetLA4R20_pt;   //!
   TBranch        *b_FFTJetLA4R20_eta;   //!
   TBranch        *b_FFTJetLA4R20_phi;   //!
   TBranch        *b_FFTJetLA4R20_m;   //!
   TBranch        *b_FFTJetLA4R20_splft;   //!
   TBranch        *b_FFTJetLA4R20_mglft;   //!
   TBranch        *b_FFTJetLA4R20_scale;   //!
   TBranch        *b_FFTJetLA2R30_pt;   //!
   TBranch        *b_FFTJetLA2R30_eta;   //!
   TBranch        *b_FFTJetLA2R30_phi;   //!
   TBranch        *b_FFTJetLA2R30_m;   //!
   TBranch        *b_FFTJetLA2R30_splft;   //!
   TBranch        *b_FFTJetLA2R30_mglft;   //!
   TBranch        *b_FFTJetLA2R30_scale;   //!
   TBranch        *b_FFTJetLA3R30_pt;   //!
   TBranch        *b_FFTJetLA3R30_eta;   //!
   TBranch        *b_FFTJetLA3R30_phi;   //!
   TBranch        *b_FFTJetLA3R30_m;   //!
   TBranch        *b_FFTJetLA3R30_splft;   //!
   TBranch        *b_FFTJetLA3R30_mglft;   //!
   TBranch        *b_FFTJetLA3R30_scale;   //!
   TBranch        *b_FFTJetLA4R30_pt;   //!
   TBranch        *b_FFTJetLA4R30_eta;   //!
   TBranch        *b_FFTJetLA4R30_phi;   //!
   TBranch        *b_FFTJetLA4R30_m;   //!
   TBranch        *b_FFTJetLA4R30_splft;   //!
   TBranch        *b_FFTJetLA4R30_mglft;   //!
   TBranch        *b_FFTJetLA4R30_scale;   //!
   TBranch        *b_FFTJetLA2R40_pt;   //!
   TBranch        *b_FFTJetLA2R40_eta;   //!
   TBranch        *b_FFTJetLA2R40_phi;   //!
   TBranch        *b_FFTJetLA2R40_m;   //!
   TBranch        *b_FFTJetLA2R40_splft;   //!
   TBranch        *b_FFTJetLA2R40_mglft;   //!
   TBranch        *b_FFTJetLA2R40_scale;   //!
   TBranch        *b_FFTJetLA3R40_pt;   //!
   TBranch        *b_FFTJetLA3R40_eta;   //!
   TBranch        *b_FFTJetLA3R40_phi;   //!
   TBranch        *b_FFTJetLA3R40_m;   //!
   TBranch        *b_FFTJetLA3R40_splft;   //!
   TBranch        *b_FFTJetLA3R40_mglft;   //!
   TBranch        *b_FFTJetLA3R40_scale;   //!
   TBranch        *b_FFTJetLA4R40_pt;   //!
   TBranch        *b_FFTJetLA4R40_eta;   //!
   TBranch        *b_FFTJetLA4R40_phi;   //!
   TBranch        *b_FFTJetLA4R40_m;   //!
   TBranch        *b_FFTJetLA4R40_splft;   //!
   TBranch        *b_FFTJetLA4R40_mglft;   //!
   TBranch        *b_FFTJetLA4R40_scale;   //!
   TBranch        *b_FFTJetLA2R50_pt;   //!
   TBranch        *b_FFTJetLA2R50_eta;   //!
   TBranch        *b_FFTJetLA2R50_phi;   //!
   TBranch        *b_FFTJetLA2R50_m;   //!
   TBranch        *b_FFTJetLA2R50_splft;   //!
   TBranch        *b_FFTJetLA2R50_mglft;   //!
   TBranch        *b_FFTJetLA2R50_scale;   //!
   TBranch        *b_FFTJetLA3R50_pt;   //!
   TBranch        *b_FFTJetLA3R50_eta;   //!
   TBranch        *b_FFTJetLA3R50_phi;   //!
   TBranch        *b_FFTJetLA3R50_m;   //!
   TBranch        *b_FFTJetLA3R50_splft;   //!
   TBranch        *b_FFTJetLA3R50_mglft;   //!
   TBranch        *b_FFTJetLA3R50_scale;   //!
   TBranch        *b_FFTJetLA4R50_pt;   //!
   TBranch        *b_FFTJetLA4R50_eta;   //!
   TBranch        *b_FFTJetLA4R50_phi;   //!
   TBranch        *b_FFTJetLA4R50_m;   //!
   TBranch        *b_FFTJetLA4R50_splft;   //!
   TBranch        *b_FFTJetLA4R50_mglft;   //!
   TBranch        *b_FFTJetLA4R50_scale;   //!
   TBranch        *b_NumPhotons;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_e;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_Theta;   //!
   TBranch        *b_Photon_Vx;   //!
   TBranch        *b_Photon_Vy;   //!
   TBranch        *b_Photon_Vz;   //!
   TBranch        *b_Photon_SC_Et;   //!
   TBranch        *b_Photon_SC_E;   //!
   TBranch        *b_Photon_SC_Eta;   //!
   TBranch        *b_Photon_SC_Phi;   //!
   TBranch        *b_Photon_SC_Theta;   //!
   TBranch        *b_Photon_SC_x;   //!
   TBranch        *b_Photon_SC_y;   //!
   TBranch        *b_Photon_SC_z;   //!
   TBranch        *b_PFisocharged03;   //!
   TBranch        *b_PFisophoton03;   //!
   TBranch        *b_PFisoneutral03;   //!
   TBranch        *b_Photon_HoverE;   //!
   TBranch        *b_Photon_HoverE2011;   //!
   TBranch        *b_Photon_SigmaIetaIeta;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!
   TBranch        *b_Photon_passElecVeto;   //!
   TBranch        *b_Photon_Id2011;   //!
   TBranch        *b_Photon_Id2012;   //!

   tmpClass(TTree *tree=0);
   virtual ~tmpClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tmpClass_cxx
tmpClass::tmpClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uscms_data/d3/tlibeiro/VPlusGamma/condorfiles/fftjetTest/VplusGTree_1_fftjetTest.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uscms_data/d3/tlibeiro/VPlusGamma/condorfiles/fftjetTest/VplusGTree_1_fftjetTest.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/uscms_data/d3/tlibeiro/VPlusGamma/condorfiles/fftjetTest/VplusGTree_1_fftjetTest.root:/vplusganalyzer_fftjettest");
      dir->GetObject("VplusGTree",tree);

   }
   Init(tree);
}

tmpClass::~tmpClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tmpClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tmpClass::LoadTree(Long64_t entry)
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

void tmpClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("event_runNo", &event_runNo, &b_event_runNo);
   fChain->SetBranchAddress("event_evtNo", &event_evtNo, &b_event_evtNo);
   fChain->SetBranchAddress("event_lumi", &event_lumi, &b_event_lumi);
   fChain->SetBranchAddress("event_bunch", &event_bunch, &b_event_bunch);
   fChain->SetBranchAddress("event_HT", &event_HT, &b_event_HT);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_m", GenPart_m, &b_GenPart_m);
   fChain->SetBranchAddress("GenPart_id", GenPart_id, &b_GenPart_id);
   fChain->SetBranchAddress("GenPart_id_mother", GenPart_id_mother, &b_GenPart_id_mother);
   fChain->SetBranchAddress("FFTJetLA2R20_pt", FFTJetLA2R20_pt, &b_FFTJetLA2R20_pt);
   fChain->SetBranchAddress("FFTJetLA2R20_eta", FFTJetLA2R20_eta, &b_FFTJetLA2R20_eta);
   fChain->SetBranchAddress("FFTJetLA2R20_phi", FFTJetLA2R20_phi, &b_FFTJetLA2R20_phi);
   fChain->SetBranchAddress("FFTJetLA2R20_m", FFTJetLA2R20_m, &b_FFTJetLA2R20_m);
   fChain->SetBranchAddress("FFTJetLA2R20_splft", FFTJetLA2R20_splft, &b_FFTJetLA2R20_splft);
   fChain->SetBranchAddress("FFTJetLA2R20_mglft", FFTJetLA2R20_mglft, &b_FFTJetLA2R20_mglft);
   fChain->SetBranchAddress("FFTJetLA2R20_scale", FFTJetLA2R20_scale, &b_FFTJetLA2R20_scale);
   fChain->SetBranchAddress("FFTJetLA3R20_pt", FFTJetLA3R20_pt, &b_FFTJetLA3R20_pt);
   fChain->SetBranchAddress("FFTJetLA3R20_eta", FFTJetLA3R20_eta, &b_FFTJetLA3R20_eta);
   fChain->SetBranchAddress("FFTJetLA3R20_phi", FFTJetLA3R20_phi, &b_FFTJetLA3R20_phi);
   fChain->SetBranchAddress("FFTJetLA3R20_m", FFTJetLA3R20_m, &b_FFTJetLA3R20_m);
   fChain->SetBranchAddress("FFTJetLA3R20_splft", FFTJetLA3R20_splft, &b_FFTJetLA3R20_splft);
   fChain->SetBranchAddress("FFTJetLA3R20_mglft", FFTJetLA3R20_mglft, &b_FFTJetLA3R20_mglft);
   fChain->SetBranchAddress("FFTJetLA3R20_scale", FFTJetLA3R20_scale, &b_FFTJetLA3R20_scale);
   fChain->SetBranchAddress("FFTJetLA4R20_pt", FFTJetLA4R20_pt, &b_FFTJetLA4R20_pt);
   fChain->SetBranchAddress("FFTJetLA4R20_eta", FFTJetLA4R20_eta, &b_FFTJetLA4R20_eta);
   fChain->SetBranchAddress("FFTJetLA4R20_phi", FFTJetLA4R20_phi, &b_FFTJetLA4R20_phi);
   fChain->SetBranchAddress("FFTJetLA4R20_m", FFTJetLA4R20_m, &b_FFTJetLA4R20_m);
   fChain->SetBranchAddress("FFTJetLA4R20_splft", FFTJetLA4R20_splft, &b_FFTJetLA4R20_splft);
   fChain->SetBranchAddress("FFTJetLA4R20_mglft", FFTJetLA4R20_mglft, &b_FFTJetLA4R20_mglft);
   fChain->SetBranchAddress("FFTJetLA4R20_scale", FFTJetLA4R20_scale, &b_FFTJetLA4R20_scale);
   fChain->SetBranchAddress("FFTJetLA2R30_pt", FFTJetLA2R30_pt, &b_FFTJetLA2R30_pt);
   fChain->SetBranchAddress("FFTJetLA2R30_eta", FFTJetLA2R30_eta, &b_FFTJetLA2R30_eta);
   fChain->SetBranchAddress("FFTJetLA2R30_phi", FFTJetLA2R30_phi, &b_FFTJetLA2R30_phi);
   fChain->SetBranchAddress("FFTJetLA2R30_m", FFTJetLA2R30_m, &b_FFTJetLA2R30_m);
   fChain->SetBranchAddress("FFTJetLA2R30_splft", FFTJetLA2R30_splft, &b_FFTJetLA2R30_splft);
   fChain->SetBranchAddress("FFTJetLA2R30_mglft", FFTJetLA2R30_mglft, &b_FFTJetLA2R30_mglft);
   fChain->SetBranchAddress("FFTJetLA2R30_scale", FFTJetLA2R30_scale, &b_FFTJetLA2R30_scale);
   fChain->SetBranchAddress("FFTJetLA3R30_pt", FFTJetLA3R30_pt, &b_FFTJetLA3R30_pt);
   fChain->SetBranchAddress("FFTJetLA3R30_eta", FFTJetLA3R30_eta, &b_FFTJetLA3R30_eta);
   fChain->SetBranchAddress("FFTJetLA3R30_phi", FFTJetLA3R30_phi, &b_FFTJetLA3R30_phi);
   fChain->SetBranchAddress("FFTJetLA3R30_m", FFTJetLA3R30_m, &b_FFTJetLA3R30_m);
   fChain->SetBranchAddress("FFTJetLA3R30_splft", FFTJetLA3R30_splft, &b_FFTJetLA3R30_splft);
   fChain->SetBranchAddress("FFTJetLA3R30_mglft", FFTJetLA3R30_mglft, &b_FFTJetLA3R30_mglft);
   fChain->SetBranchAddress("FFTJetLA3R30_scale", FFTJetLA3R30_scale, &b_FFTJetLA3R30_scale);
   fChain->SetBranchAddress("FFTJetLA4R30_pt", FFTJetLA4R30_pt, &b_FFTJetLA4R30_pt);
   fChain->SetBranchAddress("FFTJetLA4R30_eta", FFTJetLA4R30_eta, &b_FFTJetLA4R30_eta);
   fChain->SetBranchAddress("FFTJetLA4R30_phi", FFTJetLA4R30_phi, &b_FFTJetLA4R30_phi);
   fChain->SetBranchAddress("FFTJetLA4R30_m", FFTJetLA4R30_m, &b_FFTJetLA4R30_m);
   fChain->SetBranchAddress("FFTJetLA4R30_splft", FFTJetLA4R30_splft, &b_FFTJetLA4R30_splft);
   fChain->SetBranchAddress("FFTJetLA4R30_mglft", FFTJetLA4R30_mglft, &b_FFTJetLA4R30_mglft);
   fChain->SetBranchAddress("FFTJetLA4R30_scale", FFTJetLA4R30_scale, &b_FFTJetLA4R30_scale);
   fChain->SetBranchAddress("FFTJetLA2R40_pt", FFTJetLA2R40_pt, &b_FFTJetLA2R40_pt);
   fChain->SetBranchAddress("FFTJetLA2R40_eta", FFTJetLA2R40_eta, &b_FFTJetLA2R40_eta);
   fChain->SetBranchAddress("FFTJetLA2R40_phi", FFTJetLA2R40_phi, &b_FFTJetLA2R40_phi);
   fChain->SetBranchAddress("FFTJetLA2R40_m", FFTJetLA2R40_m, &b_FFTJetLA2R40_m);
   fChain->SetBranchAddress("FFTJetLA2R40_splft", FFTJetLA2R40_splft, &b_FFTJetLA2R40_splft);
   fChain->SetBranchAddress("FFTJetLA2R40_mglft", FFTJetLA2R40_mglft, &b_FFTJetLA2R40_mglft);
   fChain->SetBranchAddress("FFTJetLA2R40_scale", FFTJetLA2R40_scale, &b_FFTJetLA2R40_scale);
   fChain->SetBranchAddress("FFTJetLA3R40_pt", FFTJetLA3R40_pt, &b_FFTJetLA3R40_pt);
   fChain->SetBranchAddress("FFTJetLA3R40_eta", FFTJetLA3R40_eta, &b_FFTJetLA3R40_eta);
   fChain->SetBranchAddress("FFTJetLA3R40_phi", FFTJetLA3R40_phi, &b_FFTJetLA3R40_phi);
   fChain->SetBranchAddress("FFTJetLA3R40_m", FFTJetLA3R40_m, &b_FFTJetLA3R40_m);
   fChain->SetBranchAddress("FFTJetLA3R40_splft", FFTJetLA3R40_splft, &b_FFTJetLA3R40_splft);
   fChain->SetBranchAddress("FFTJetLA3R40_mglft", FFTJetLA3R40_mglft, &b_FFTJetLA3R40_mglft);
   fChain->SetBranchAddress("FFTJetLA3R40_scale", FFTJetLA3R40_scale, &b_FFTJetLA3R40_scale);
   fChain->SetBranchAddress("FFTJetLA4R40_pt", FFTJetLA4R40_pt, &b_FFTJetLA4R40_pt);
   fChain->SetBranchAddress("FFTJetLA4R40_eta", FFTJetLA4R40_eta, &b_FFTJetLA4R40_eta);
   fChain->SetBranchAddress("FFTJetLA4R40_phi", FFTJetLA4R40_phi, &b_FFTJetLA4R40_phi);
   fChain->SetBranchAddress("FFTJetLA4R40_m", FFTJetLA4R40_m, &b_FFTJetLA4R40_m);
   fChain->SetBranchAddress("FFTJetLA4R40_splft", FFTJetLA4R40_splft, &b_FFTJetLA4R40_splft);
   fChain->SetBranchAddress("FFTJetLA4R40_mglft", FFTJetLA4R40_mglft, &b_FFTJetLA4R40_mglft);
   fChain->SetBranchAddress("FFTJetLA4R40_scale", FFTJetLA4R40_scale, &b_FFTJetLA4R40_scale);
   fChain->SetBranchAddress("FFTJetLA2R50_pt", FFTJetLA2R50_pt, &b_FFTJetLA2R50_pt);
   fChain->SetBranchAddress("FFTJetLA2R50_eta", FFTJetLA2R50_eta, &b_FFTJetLA2R50_eta);
   fChain->SetBranchAddress("FFTJetLA2R50_phi", FFTJetLA2R50_phi, &b_FFTJetLA2R50_phi);
   fChain->SetBranchAddress("FFTJetLA2R50_m", FFTJetLA2R50_m, &b_FFTJetLA2R50_m);
   fChain->SetBranchAddress("FFTJetLA2R50_splft", FFTJetLA2R50_splft, &b_FFTJetLA2R50_splft);
   fChain->SetBranchAddress("FFTJetLA2R50_mglft", FFTJetLA2R50_mglft, &b_FFTJetLA2R50_mglft);
   fChain->SetBranchAddress("FFTJetLA2R50_scale", FFTJetLA2R50_scale, &b_FFTJetLA2R50_scale);
   fChain->SetBranchAddress("FFTJetLA3R50_pt", FFTJetLA3R50_pt, &b_FFTJetLA3R50_pt);
   fChain->SetBranchAddress("FFTJetLA3R50_eta", FFTJetLA3R50_eta, &b_FFTJetLA3R50_eta);
   fChain->SetBranchAddress("FFTJetLA3R50_phi", FFTJetLA3R50_phi, &b_FFTJetLA3R50_phi);
   fChain->SetBranchAddress("FFTJetLA3R50_m", FFTJetLA3R50_m, &b_FFTJetLA3R50_m);
   fChain->SetBranchAddress("FFTJetLA3R50_splft", FFTJetLA3R50_splft, &b_FFTJetLA3R50_splft);
   fChain->SetBranchAddress("FFTJetLA3R50_mglft", FFTJetLA3R50_mglft, &b_FFTJetLA3R50_mglft);
   fChain->SetBranchAddress("FFTJetLA3R50_scale", FFTJetLA3R50_scale, &b_FFTJetLA3R50_scale);
   fChain->SetBranchAddress("FFTJetLA4R50_pt", FFTJetLA4R50_pt, &b_FFTJetLA4R50_pt);
   fChain->SetBranchAddress("FFTJetLA4R50_eta", FFTJetLA4R50_eta, &b_FFTJetLA4R50_eta);
   fChain->SetBranchAddress("FFTJetLA4R50_phi", FFTJetLA4R50_phi, &b_FFTJetLA4R50_phi);
   fChain->SetBranchAddress("FFTJetLA4R50_m", FFTJetLA4R50_m, &b_FFTJetLA4R50_m);
   fChain->SetBranchAddress("FFTJetLA4R50_splft", FFTJetLA4R50_splft, &b_FFTJetLA4R50_splft);
   fChain->SetBranchAddress("FFTJetLA4R50_mglft", FFTJetLA4R50_mglft, &b_FFTJetLA4R50_mglft);
   fChain->SetBranchAddress("FFTJetLA4R50_scale", FFTJetLA4R50_scale, &b_FFTJetLA4R50_scale);
   fChain->SetBranchAddress("NumPhotons", &NumPhotons, &b_NumPhotons);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_e", Photon_e, &b_Photon_e);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_Theta", Photon_Theta, &b_Photon_Theta);
   fChain->SetBranchAddress("Photon_Vx", Photon_Vx, &b_Photon_Vx);
   fChain->SetBranchAddress("Photon_Vy", Photon_Vy, &b_Photon_Vy);
   fChain->SetBranchAddress("Photon_Vz", Photon_Vz, &b_Photon_Vz);
   fChain->SetBranchAddress("Photon_SC_Et", Photon_SC_Et, &b_Photon_SC_Et);
   fChain->SetBranchAddress("Photon_SC_E", Photon_SC_E, &b_Photon_SC_E);
   fChain->SetBranchAddress("Photon_SC_Eta", Photon_SC_Eta, &b_Photon_SC_Eta);
   fChain->SetBranchAddress("Photon_SC_Phi", Photon_SC_Phi, &b_Photon_SC_Phi);
   fChain->SetBranchAddress("Photon_SC_Theta", Photon_SC_Theta, &b_Photon_SC_Theta);
   fChain->SetBranchAddress("Photon_SC_x", Photon_SC_x, &b_Photon_SC_x);
   fChain->SetBranchAddress("Photon_SC_y", Photon_SC_y, &b_Photon_SC_y);
   fChain->SetBranchAddress("Photon_SC_z", Photon_SC_z, &b_Photon_SC_z);
   fChain->SetBranchAddress("PFisocharged03", PFisocharged03, &b_PFisocharged03);
   fChain->SetBranchAddress("PFisophoton03", PFisophoton03, &b_PFisophoton03);
   fChain->SetBranchAddress("PFisoneutral03", PFisoneutral03, &b_PFisoneutral03);
   fChain->SetBranchAddress("Photon_HoverE", Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon_HoverE2011", Photon_HoverE2011, &b_Photon_HoverE2011);
   fChain->SetBranchAddress("Photon_SigmaIetaIeta", Photon_SigmaIetaIeta, &b_Photon_SigmaIetaIeta);
   fChain->SetBranchAddress("Photon_hasPixelSeed", Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   fChain->SetBranchAddress("Photon_passElecVeto", Photon_passElecVeto, &b_Photon_passElecVeto);
   fChain->SetBranchAddress("Photon_Id2011", Photon_Id2011, &b_Photon_Id2011);
   fChain->SetBranchAddress("Photon_Id2012", Photon_Id2012, &b_Photon_Id2012);
   Notify();
}

Bool_t tmpClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tmpClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tmpClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tmpClass_cxx
