//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 30 15:04:16 2013 by ROOT version 5.32/00
// from TTree VplusGTree/VplusGTree
// found on file: VplusGTree.root
//////////////////////////////////////////////////////////

#ifndef VplusGTreeAnalyzer_h
#define VplusGTreeAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <sstream>
#include "TLorentzVector.h"
#include <TDCacheFile.h>
using namespace std;
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class VplusGTreeAnalyzer {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain         *chain;   //!pointer to the analyzed TTree or TChain
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
   Float_t         GroomedJet_CA8_pt[4];
   Float_t         GroomedJet_CA8_eta[4];
   Float_t         GroomedJet_CA8_phi[4];
   Float_t         GroomedJet_CA8_m[4];
   Float_t         GroomedJet_CA8_e[4];
   Float_t         GroomedJet_CA8_ar[4];
   Float_t         GroomedJet_CA8_pt_uncorr[4];
   Float_t         GroomedJet_CA8_m_uncorr[4];
   Float_t         GroomedJet_CA8_pt_tr[4];
   Float_t         GroomedJet_CA8_eta_tr[4];
   Float_t         GroomedJet_CA8_phi_tr[4];
   Float_t         GroomedJet_CA8_m_tr[4];
   Float_t         GroomedJet_CA8_e_tr[4];
   Float_t         GroomedJet_CA8_ar_tr[4];
   Float_t         GroomedJet_CA8_pt_ft[4];
   Float_t         GroomedJet_CA8_eta_ft[4];
   Float_t         GroomedJet_CA8_phi_ft[4];
   Float_t         GroomedJet_CA8_m_ft[4];
   Float_t         GroomedJet_CA8_e_ft[4];
   Float_t         GroomedJet_CA8_ar_ft[4];
   Float_t         GroomedJet_CA8_pt_pr[4];
   Float_t         GroomedJet_CA8_eta_pr[4];
   Float_t         GroomedJet_CA8_phi_pr[4];
   Float_t         GroomedJet_CA8_m_pr[4];
   Float_t         GroomedJet_CA8_e_pr[4];
   Float_t         GroomedJet_CA8_ar_pr[4];
   Float_t         GroomedJet_CA8_tau1[4];
   Float_t         GroomedJet_CA8_tau2[4];
   Float_t         GroomedJet_CA8_tau3[4];
   Float_t         GroomedJet_CA8_tau4[4];
   Float_t         GroomedJet_CA8_tau2tau1[4];
   Float_t         GroomedJet_CA8_pt_ft_uncorr[4];
   Float_t         GroomedJet_CA8_m_ft_uncorr[4];
   Float_t         GroomedJet_CA8_pt_tr_uncorr[4];
   Float_t         GroomedJet_CA8_m_tr_uncorr[4];
   Float_t         GroomedJet_CA8_pt_pr_uncorr[4];
   Float_t         GroomedJet_CA8_m_pr_uncorr[4];
   Float_t         GroomedJet_CA8_pt_prsubj1[4];
   Float_t         GroomedJet_CA8_eta_prsubj1[4];
   Float_t         GroomedJet_CA8_phi_prsubj1[4];
   Float_t         GroomedJet_CA8_m_prsubj1[4];
   Float_t         GroomedJet_CA8_pt_prsubj2[4];
   Float_t         GroomedJet_CA8_eta_prsubj2[4];
   Float_t         GroomedJet_CA8_phi_prsubj2[4];
   Float_t         GroomedJet_CA8_m_prsubj2[4];
   Float_t         GenGroomedJet_CA8_pt[4];
   Float_t         GenGroomedJet_CA8_eta[4];
   Float_t         GenGroomedJet_CA8_phi[4];
   Float_t         GenGroomedJet_CA8_m[4];
   Float_t         GenGroomedJet_CA8_e[4];
   Float_t         GenGroomedJet_CA8_ar[4];
   Float_t         GenGroomedJet_CA8_pt_uncorr[4];
   Float_t         GenGroomedJet_CA8_m_uncorr[4];
   Float_t         GenGroomedJet_CA8_pt_tr[4];
   Float_t         GenGroomedJet_CA8_eta_tr[4];
   Float_t         GenGroomedJet_CA8_phi_tr[4];
   Float_t         GenGroomedJet_CA8_m_tr[4];
   Float_t         GenGroomedJet_CA8_e_tr[4];
   Float_t         GenGroomedJet_CA8_ar_tr[4];
   Float_t         GenGroomedJet_CA8_pt_ft[4];
   Float_t         GenGroomedJet_CA8_eta_ft[4];
   Float_t         GenGroomedJet_CA8_phi_ft[4];
   Float_t         GenGroomedJet_CA8_m_ft[4];
   Float_t         GenGroomedJet_CA8_e_ft[4];
   Float_t         GenGroomedJet_CA8_ar_ft[4];
   Float_t         GenGroomedJet_CA8_pt_pr[4];
   Float_t         GenGroomedJet_CA8_eta_pr[4];
   Float_t         GenGroomedJet_CA8_phi_pr[4];
   Float_t         GenGroomedJet_CA8_m_pr[4];
   Float_t         GenGroomedJet_CA8_e_pr[4];
   Float_t         GenGroomedJet_CA8_ar_pr[4];
   Float_t         GenGroomedJet_CA8_tau1[4];
   Float_t         GenGroomedJet_CA8_tau2[4];
   Float_t         GenGroomedJet_CA8_tau3[4];
   Float_t         GenGroomedJet_CA8_tau4[4];
   Float_t         GenGroomedJet_CA8_tau2tau1[4];
   Float_t         GenGroomedJet_CA8_pt_prsubj1[4];
   Float_t         GenGroomedJet_CA8_eta_prsubj1[4];
   Float_t         GenGroomedJet_CA8_phi_prsubj1[4];
   Float_t         GenGroomedJet_CA8_m_prsubj1[4];
   Float_t         GenGroomedJet_CA8_pt_prsubj2[4];
   Float_t         GenGroomedJet_CA8_eta_prsubj2[4];
   Float_t         GenGroomedJet_CA8_phi_prsubj2[4];
   Float_t         GenGroomedJet_CA8_m_prsubj2[4];
   Float_t         FFTJetR20_pt[6];
   Float_t         FFTJetR20_eta[6];
   Float_t         FFTJetR20_phi[6];
   Float_t         FFTJetR20_m[6];
   Float_t         FFTJetR20_splft[6];
   Float_t         FFTJetR20_mglft[6];
   Float_t         FFTJetR20_scale[6];
   Float_t         FFTJetR30_pt[6];
   Float_t         FFTJetR30_eta[6];
   Float_t         FFTJetR30_phi[6];
   Float_t         FFTJetR30_m[6];
   Float_t         FFTJetR30_splft[6];
   Float_t         FFTJetR30_mglft[6];
   Float_t         FFTJetR30_scale[6];
   Float_t         FFTJetR40_pt[6];
   Float_t         FFTJetR40_eta[6];
   Float_t         FFTJetR40_phi[6];
   Float_t         FFTJetR40_m[6];
   Float_t         FFTJetR40_splft[6];
   Float_t         FFTJetR40_mglft[6];
   Float_t         FFTJetR40_scale[6];
   Int_t           NumPhotons;
   Float_t         Photon_pt[6];   //[NumPhotons]
   Float_t         Photon_e[6];   //[NumPhotons]
   Float_t         Photon_eta[6];   //[NumPhotons]
   Float_t         Photon_phi[6];   //[NumPhotons]
   Float_t         Photon_Theta[6];   //[NumPhotons]
   Float_t         Photon_Vx[6];   //[NumPhotons]
   Float_t         Photon_Vy[6];   //[NumPhotons]
   Float_t         Photon_Vz[6];   //[NumPhotons]
   Float_t         Photon_SC_Et[6];   //[NumPhotons]
   Float_t         Photon_SC_E[6];   //[NumPhotons]
   Float_t         Photon_SC_Eta[6];   //[NumPhotons]
   Float_t         Photon_SC_Phi[6];   //[NumPhotons]
   Float_t         Photon_SC_Theta[6];   //[NumPhotons]
   Float_t         Photon_SC_x[6];   //[NumPhotons]
   Float_t         Photon_SC_y[6];   //[NumPhotons]
   Float_t         Photon_SC_z[6];   //[NumPhotons]
   Float_t         PFisocharged03[6];   //[NumPhotons]
   Float_t         PFisophoton03[6];   //[NumPhotons]
   Float_t         PFisoneutral03[6];   //[NumPhotons]
   Float_t         Photon_HoverE[6];   //[NumPhotons]
   Float_t         Photon_HoverE2011[6];   //[NumPhotons]
   Float_t         Photon_SigmaIetaIeta[6];   //[NumPhotons]
   Int_t           Photon_hasPixelSeed[6];   //[NumPhotons]
   Int_t           Photon_passElecVeto[6];   //[NumPhotons]
   Int_t           Photon_Id2011[6];   //[NumPhotons]
   Int_t           Photon_Id2012[6];   //[NumPhotons]

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
   TBranch        *b_GroomedJet_CA8_pt;   //!
   TBranch        *b_GroomedJet_CA8_eta;   //!
   TBranch        *b_GroomedJet_CA8_phi;   //!
   TBranch        *b_GroomedJet_CA8_m;   //!
   TBranch        *b_GroomedJet_CA8_e;   //!
   TBranch        *b_GroomedJet_CA8_ar;   //!
   TBranch        *b_GroomedJet_CA8_pt_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_m_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_pt_tr;   //!
   TBranch        *b_GroomedJet_CA8_eta_tr;   //!
   TBranch        *b_GroomedJet_CA8_phi_tr;   //!
   TBranch        *b_GroomedJet_CA8_m_tr;   //!
   TBranch        *b_GroomedJet_CA8_e_tr;   //!
   TBranch        *b_GroomedJet_CA8_ar_tr;   //!
   TBranch        *b_GroomedJet_CA8_pt_ft;   //!
   TBranch        *b_GroomedJet_CA8_eta_ft;   //!
   TBranch        *b_GroomedJet_CA8_phi_ft;   //!
   TBranch        *b_GroomedJet_CA8_m_ft;   //!
   TBranch        *b_GroomedJet_CA8_e_ft;   //!
   TBranch        *b_GroomedJet_CA8_ar_ft;   //!
   TBranch        *b_GroomedJet_CA8_pt_pr;   //!
   TBranch        *b_GroomedJet_CA8_eta_pr;   //!
   TBranch        *b_GroomedJet_CA8_phi_pr;   //!
   TBranch        *b_GroomedJet_CA8_m_pr;   //!
   TBranch        *b_GroomedJet_CA8_e_pr;   //!
   TBranch        *b_GroomedJet_CA8_ar_pr;   //!
   TBranch        *b_GroomedJet_CA8_tau1;   //!
   TBranch        *b_GroomedJet_CA8_tau2;   //!
   TBranch        *b_GroomedJet_CA8_tau3;   //!
   TBranch        *b_GroomedJet_CA8_tau4;   //!
   TBranch        *b_GroomedJet_CA8_tau2tau1;   //!
   TBranch        *b_GroomedJet_CA8_pt_ft_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_m_ft_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_pt_tr_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_m_tr_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_pt_pr_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_m_pr_uncorr;   //!
   TBranch        *b_GroomedJet_CA8_pt_prsubj1;   //!
   TBranch        *b_GroomedJet_CA8_eta_prsubj1;   //!
   TBranch        *b_GroomedJet_CA8_phi_prsubj1;   //!
   TBranch        *b_GroomedJet_CA8_m_prsubj1;   //!
   TBranch        *b_GroomedJet_CA8_pt_prsubj2;   //!
   TBranch        *b_GroomedJet_CA8_eta_prsubj2;   //!
   TBranch        *b_GroomedJet_CA8_phi_prsubj2;   //!
   TBranch        *b_GroomedJet_CA8_m_prsubj2;   //!
   TBranch        *b_GenGroomedJet_CA8_pt;   //!
   TBranch        *b_GenGroomedJet_CA8_eta;   //!
   TBranch        *b_GenGroomedJet_CA8_phi;   //!
   TBranch        *b_GenGroomedJet_CA8_m;   //!
   TBranch        *b_GenGroomedJet_CA8_e;   //!
   TBranch        *b_GenGroomedJet_CA8_ar;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_uncorr;   //!
   TBranch        *b_GenGroomedJet_CA8_m_uncorr;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_eta_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_phi_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_m_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_e_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_ar_tr;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_eta_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_phi_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_m_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_e_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_ar_ft;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_eta_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_phi_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_m_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_e_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_ar_pr;   //!
   TBranch        *b_GenGroomedJet_CA8_tau1;   //!
   TBranch        *b_GenGroomedJet_CA8_tau2;   //!
   TBranch        *b_GenGroomedJet_CA8_tau3;   //!
   TBranch        *b_GenGroomedJet_CA8_tau4;   //!
   TBranch        *b_GenGroomedJet_CA8_tau2tau1;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_prsubj1;   //!
   TBranch        *b_GenGroomedJet_CA8_eta_prsubj1;   //!
   TBranch        *b_GenGroomedJet_CA8_phi_prsubj1;   //!
   TBranch        *b_GenGroomedJet_CA8_m_prsubj1;   //!
   TBranch        *b_GenGroomedJet_CA8_pt_prsubj2;   //!
   TBranch        *b_GenGroomedJet_CA8_eta_prsubj2;   //!
   TBranch        *b_GenGroomedJet_CA8_phi_prsubj2;   //!
   TBranch        *b_GenGroomedJet_CA8_m_prsubj2;   //!
   TBranch        *b_FFTJetR20_pt;   //!
   TBranch        *b_FFTJetR20_eta;   //!
   TBranch        *b_FFTJetR20_phi;   //!
   TBranch        *b_FFTJetR20_m;   //!
   TBranch        *b_FFTJetR20_splft;   //!
   TBranch        *b_FFTJetR20_mglft;   //!
   TBranch        *b_FFTJetR20_scale;   //!
   TBranch        *b_FFTJetR30_pt;   //!
   TBranch        *b_FFTJetR30_eta;   //!
   TBranch        *b_FFTJetR30_phi;   //!
   TBranch        *b_FFTJetR30_m;   //!
   TBranch        *b_FFTJetR30_splft;   //!
   TBranch        *b_FFTJetR30_mglft;   //!
   TBranch        *b_FFTJetR30_scale;   //!
   TBranch        *b_FFTJetR40_pt;   //!
   TBranch        *b_FFTJetR40_eta;   //!
   TBranch        *b_FFTJetR40_phi;   //!
   TBranch        *b_FFTJetR40_m;   //!
   TBranch        *b_FFTJetR40_splft;   //!
   TBranch        *b_FFTJetR40_mglft;   //!
   TBranch        *b_FFTJetR40_scale;   //!
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

   VplusGTreeAnalyzer(const vector<string>& infile, const string& outfile, TTree *tree=0);
   virtual ~VplusGTreeAnalyzer();
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

#ifdef VplusGTreeAnalyzer_cxx
VplusGTreeAnalyzer::VplusGTreeAnalyzer(const vector<string>& infiles, const string& ofile, TTree *tree) : fChain(0),outfile(ofile) 
{
	TChain * chain = new TChain("vplusganalyzer/VplusGTree","");
	const unsigned numinput(infiles.size());
	for(unsigned i(0);i<numinput;++i)
	{
		string infile(infiles[i]);
		if(infile.find("pnfs")!=string::npos)
			infile = "dcache:"+infile;
		chain->Add(infile.c_str());
	}
	tree = chain;
	Init(tree);

	ostringstream ntup_vars; 
	vector<string> vars; 
	vars.push_back("run"); vars.push_back("event");   vars.push_back("lumi");
	vars.push_back("HT");
	vars.push_back("m_V_gr");
	vars.push_back("pt_V_gr");
	vars.push_back("pt_V_ungr");
	vars.push_back("eta_V_gr");
	vars.push_back("pt_G");
	vars.push_back("eta_G");
	vars.push_back("dr_VG");
	vars.push_back("dphi_VG");
	vars.push_back("m_VG");
	vars.push_back("m_V_fft");
	vars.push_back("pt_V_fft");
	vars.push_back("eta_V_fft");
	vars.push_back("maxdrjj_V_fft");
	vars.push_back("dr_VG_fft");
	vars.push_back("dphi_VG_fft");
	vars.push_back("m_VG_fft");
	vars.push_back("photonIso");
	vars.push_back("tau2tau1");
	const unsigned numvars(22);
	for(unsigned i(0);i<numvars;++i)
		if(!i) ntup_vars<<vars[i];
		else ntup_vars<<':'<<vars[i];

	ntup_VG = new TNtuple("ntup_VG","ntup_VG",ntup_vars.str().c_str());
}

VplusGTreeAnalyzer::~VplusGTreeAnalyzer()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t VplusGTreeAnalyzer::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t VplusGTreeAnalyzer::LoadTree(Long64_t entry)
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

void VplusGTreeAnalyzer::Init(TTree *tree)
{
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
	fChain->SetBranchAddress("GroomedJet_CA8_pt", GroomedJet_CA8_pt, &b_GroomedJet_CA8_pt);
	fChain->SetBranchAddress("GroomedJet_CA8_eta", GroomedJet_CA8_eta, &b_GroomedJet_CA8_eta);
	fChain->SetBranchAddress("GroomedJet_CA8_phi", GroomedJet_CA8_phi, &b_GroomedJet_CA8_phi);
	fChain->SetBranchAddress("GroomedJet_CA8_m", GroomedJet_CA8_m, &b_GroomedJet_CA8_m);
	fChain->SetBranchAddress("GroomedJet_CA8_e", GroomedJet_CA8_e, &b_GroomedJet_CA8_e);
	fChain->SetBranchAddress("GroomedJet_CA8_ar", GroomedJet_CA8_ar, &b_GroomedJet_CA8_ar);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_uncorr", GroomedJet_CA8_pt_uncorr, &b_GroomedJet_CA8_pt_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_uncorr", GroomedJet_CA8_m_uncorr, &b_GroomedJet_CA8_m_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_tr", GroomedJet_CA8_pt_tr, &b_GroomedJet_CA8_pt_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_tr", GroomedJet_CA8_eta_tr, &b_GroomedJet_CA8_eta_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_tr", GroomedJet_CA8_phi_tr, &b_GroomedJet_CA8_phi_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_tr", GroomedJet_CA8_m_tr, &b_GroomedJet_CA8_m_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_e_tr", GroomedJet_CA8_e_tr, &b_GroomedJet_CA8_e_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_ar_tr", GroomedJet_CA8_ar_tr, &b_GroomedJet_CA8_ar_tr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_ft", GroomedJet_CA8_pt_ft, &b_GroomedJet_CA8_pt_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_ft", GroomedJet_CA8_eta_ft, &b_GroomedJet_CA8_eta_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_ft", GroomedJet_CA8_phi_ft, &b_GroomedJet_CA8_phi_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_m_ft", GroomedJet_CA8_m_ft, &b_GroomedJet_CA8_m_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_e_ft", GroomedJet_CA8_e_ft, &b_GroomedJet_CA8_e_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_ar_ft", GroomedJet_CA8_ar_ft, &b_GroomedJet_CA8_ar_ft);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_pr", GroomedJet_CA8_pt_pr, &b_GroomedJet_CA8_pt_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_pr", GroomedJet_CA8_eta_pr, &b_GroomedJet_CA8_eta_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_pr", GroomedJet_CA8_phi_pr, &b_GroomedJet_CA8_phi_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_pr", GroomedJet_CA8_m_pr, &b_GroomedJet_CA8_m_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_e_pr", GroomedJet_CA8_e_pr, &b_GroomedJet_CA8_e_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_ar_pr", GroomedJet_CA8_ar_pr, &b_GroomedJet_CA8_ar_pr);
	fChain->SetBranchAddress("GroomedJet_CA8_tau1", GroomedJet_CA8_tau1, &b_GroomedJet_CA8_tau1);
	fChain->SetBranchAddress("GroomedJet_CA8_tau2", GroomedJet_CA8_tau2, &b_GroomedJet_CA8_tau2);
	fChain->SetBranchAddress("GroomedJet_CA8_tau3", GroomedJet_CA8_tau3, &b_GroomedJet_CA8_tau3);
	fChain->SetBranchAddress("GroomedJet_CA8_tau4", GroomedJet_CA8_tau4, &b_GroomedJet_CA8_tau4);
	fChain->SetBranchAddress("GroomedJet_CA8_tau2tau1", GroomedJet_CA8_tau2tau1, &b_GroomedJet_CA8_tau2tau1);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_ft_uncorr", GroomedJet_CA8_pt_ft_uncorr, &b_GroomedJet_CA8_pt_ft_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_ft_uncorr", GroomedJet_CA8_m_ft_uncorr, &b_GroomedJet_CA8_m_ft_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_tr_uncorr", GroomedJet_CA8_pt_tr_uncorr, &b_GroomedJet_CA8_pt_tr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_tr_uncorr", GroomedJet_CA8_m_tr_uncorr, &b_GroomedJet_CA8_m_tr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_pr_uncorr", GroomedJet_CA8_pt_pr_uncorr, &b_GroomedJet_CA8_pt_pr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_m_pr_uncorr", GroomedJet_CA8_m_pr_uncorr, &b_GroomedJet_CA8_m_pr_uncorr);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_prsubj1", GroomedJet_CA8_pt_prsubj1, &b_GroomedJet_CA8_pt_prsubj1);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_prsubj1", GroomedJet_CA8_eta_prsubj1, &b_GroomedJet_CA8_eta_prsubj1);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_prsubj1", GroomedJet_CA8_phi_prsubj1, &b_GroomedJet_CA8_phi_prsubj1);
	fChain->SetBranchAddress("GroomedJet_CA8_m_prsubj1", GroomedJet_CA8_m_prsubj1, &b_GroomedJet_CA8_m_prsubj1);
	fChain->SetBranchAddress("GroomedJet_CA8_pt_prsubj2", GroomedJet_CA8_pt_prsubj2, &b_GroomedJet_CA8_pt_prsubj2);
	fChain->SetBranchAddress("GroomedJet_CA8_eta_prsubj2", GroomedJet_CA8_eta_prsubj2, &b_GroomedJet_CA8_eta_prsubj2);
	fChain->SetBranchAddress("GroomedJet_CA8_phi_prsubj2", GroomedJet_CA8_phi_prsubj2, &b_GroomedJet_CA8_phi_prsubj2);
	fChain->SetBranchAddress("GroomedJet_CA8_m_prsubj2", GroomedJet_CA8_m_prsubj2, &b_GroomedJet_CA8_m_prsubj2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt", GenGroomedJet_CA8_pt, &b_GenGroomedJet_CA8_pt);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta", GenGroomedJet_CA8_eta, &b_GenGroomedJet_CA8_eta);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi", GenGroomedJet_CA8_phi, &b_GenGroomedJet_CA8_phi);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m", GenGroomedJet_CA8_m, &b_GenGroomedJet_CA8_m);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e", GenGroomedJet_CA8_e, &b_GenGroomedJet_CA8_e);
	fChain->SetBranchAddress("GenGroomedJet_CA8_ar", GenGroomedJet_CA8_ar, &b_GenGroomedJet_CA8_ar);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_uncorr", GenGroomedJet_CA8_pt_uncorr, &b_GenGroomedJet_CA8_pt_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_uncorr", GenGroomedJet_CA8_m_uncorr, &b_GenGroomedJet_CA8_m_uncorr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_tr", GenGroomedJet_CA8_pt_tr, &b_GenGroomedJet_CA8_pt_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_tr", GenGroomedJet_CA8_eta_tr, &b_GenGroomedJet_CA8_eta_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_tr", GenGroomedJet_CA8_phi_tr, &b_GenGroomedJet_CA8_phi_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_tr", GenGroomedJet_CA8_m_tr, &b_GenGroomedJet_CA8_m_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_tr", GenGroomedJet_CA8_e_tr, &b_GenGroomedJet_CA8_e_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_ar_tr", GenGroomedJet_CA8_ar_tr, &b_GenGroomedJet_CA8_ar_tr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_ft", GenGroomedJet_CA8_pt_ft, &b_GenGroomedJet_CA8_pt_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_ft", GenGroomedJet_CA8_eta_ft, &b_GenGroomedJet_CA8_eta_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_ft", GenGroomedJet_CA8_phi_ft, &b_GenGroomedJet_CA8_phi_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_ft", GenGroomedJet_CA8_m_ft, &b_GenGroomedJet_CA8_m_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_ft", GenGroomedJet_CA8_e_ft, &b_GenGroomedJet_CA8_e_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_ar_ft", GenGroomedJet_CA8_ar_ft, &b_GenGroomedJet_CA8_ar_ft);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_pr", GenGroomedJet_CA8_pt_pr, &b_GenGroomedJet_CA8_pt_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_pr", GenGroomedJet_CA8_eta_pr, &b_GenGroomedJet_CA8_eta_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_pr", GenGroomedJet_CA8_phi_pr, &b_GenGroomedJet_CA8_phi_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_pr", GenGroomedJet_CA8_m_pr, &b_GenGroomedJet_CA8_m_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_e_pr", GenGroomedJet_CA8_e_pr, &b_GenGroomedJet_CA8_e_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_ar_pr", GenGroomedJet_CA8_ar_pr, &b_GenGroomedJet_CA8_ar_pr);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau1", GenGroomedJet_CA8_tau1, &b_GenGroomedJet_CA8_tau1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau2", GenGroomedJet_CA8_tau2, &b_GenGroomedJet_CA8_tau2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau3", GenGroomedJet_CA8_tau3, &b_GenGroomedJet_CA8_tau3);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau4", GenGroomedJet_CA8_tau4, &b_GenGroomedJet_CA8_tau4);
	fChain->SetBranchAddress("GenGroomedJet_CA8_tau2tau1", GenGroomedJet_CA8_tau2tau1, &b_GenGroomedJet_CA8_tau2tau1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_prsubj1", GenGroomedJet_CA8_pt_prsubj1, &b_GenGroomedJet_CA8_pt_prsubj1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_prsubj1", GenGroomedJet_CA8_eta_prsubj1, &b_GenGroomedJet_CA8_eta_prsubj1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_prsubj1", GenGroomedJet_CA8_phi_prsubj1, &b_GenGroomedJet_CA8_phi_prsubj1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_prsubj1", GenGroomedJet_CA8_m_prsubj1, &b_GenGroomedJet_CA8_m_prsubj1);
	fChain->SetBranchAddress("GenGroomedJet_CA8_pt_prsubj2", GenGroomedJet_CA8_pt_prsubj2, &b_GenGroomedJet_CA8_pt_prsubj2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_eta_prsubj2", GenGroomedJet_CA8_eta_prsubj2, &b_GenGroomedJet_CA8_eta_prsubj2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_phi_prsubj2", GenGroomedJet_CA8_phi_prsubj2, &b_GenGroomedJet_CA8_phi_prsubj2);
	fChain->SetBranchAddress("GenGroomedJet_CA8_m_prsubj2", GenGroomedJet_CA8_m_prsubj2, &b_GenGroomedJet_CA8_m_prsubj2);
	fChain->SetBranchAddress("FFTJetR20_pt", FFTJetR20_pt, &b_FFTJetR20_pt);
	fChain->SetBranchAddress("FFTJetR20_eta", FFTJetR20_eta, &b_FFTJetR20_eta);
	fChain->SetBranchAddress("FFTJetR20_phi", FFTJetR20_phi, &b_FFTJetR20_phi);
	fChain->SetBranchAddress("FFTJetR20_m", FFTJetR20_m, &b_FFTJetR20_m);
	fChain->SetBranchAddress("FFTJetR20_splft", FFTJetR20_splft, &b_FFTJetR20_splft);
	fChain->SetBranchAddress("FFTJetR20_mglft", FFTJetR20_mglft, &b_FFTJetR20_mglft);
	fChain->SetBranchAddress("FFTJetR20_scale", FFTJetR20_scale, &b_FFTJetR20_scale);
	fChain->SetBranchAddress("FFTJetR30_pt", FFTJetR30_pt, &b_FFTJetR30_pt);
	fChain->SetBranchAddress("FFTJetR30_eta", FFTJetR30_eta, &b_FFTJetR30_eta);
	fChain->SetBranchAddress("FFTJetR30_phi", FFTJetR30_phi, &b_FFTJetR30_phi);
	fChain->SetBranchAddress("FFTJetR30_m", FFTJetR30_m, &b_FFTJetR30_m);
	fChain->SetBranchAddress("FFTJetR30_splft", FFTJetR30_splft, &b_FFTJetR30_splft);
	fChain->SetBranchAddress("FFTJetR30_mglft", FFTJetR30_mglft, &b_FFTJetR30_mglft);
	fChain->SetBranchAddress("FFTJetR30_scale", FFTJetR30_scale, &b_FFTJetR30_scale);
	fChain->SetBranchAddress("FFTJetR40_pt", FFTJetR40_pt, &b_FFTJetR40_pt);
	fChain->SetBranchAddress("FFTJetR40_eta", FFTJetR40_eta, &b_FFTJetR40_eta);
	fChain->SetBranchAddress("FFTJetR40_phi", FFTJetR40_phi, &b_FFTJetR40_phi);
	fChain->SetBranchAddress("FFTJetR40_m", FFTJetR40_m, &b_FFTJetR40_m);
	fChain->SetBranchAddress("FFTJetR40_splft", FFTJetR40_splft, &b_FFTJetR40_splft);
	fChain->SetBranchAddress("FFTJetR40_mglft", FFTJetR40_mglft, &b_FFTJetR40_mglft);
	fChain->SetBranchAddress("FFTJetR40_scale", FFTJetR40_scale, &b_FFTJetR40_scale);
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

Bool_t VplusGTreeAnalyzer::Notify()
{
	return kTRUE;
}

void VplusGTreeAnalyzer::Show(Long64_t entry)
{
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t VplusGTreeAnalyzer::Cut(Long64_t entry)
{
	return 1;
}
//helper classes and functions
class GenParticle {
	public:
		double pt;
		double eta;
		double phi;
		double mass;
		double id; 
		double mother_id;
		//constructor
		GenParticle(double p, double e,double ph , double m,
				double i, double mid)
			:pt(p), eta(e), phi(ph), mass(m), id(i), 
			mother_id(mid) 
	{};
};
class Jet {
	public : 
		double pt;
		double eta;
		double phi;
		double mass;
		double mglft;
		double sptlft;
		double scale;
		double photonIso;
		//constructors
		Jet(double p, double e, double ph, double m) 
			:mglft(0), sptlft(0), scale(0),photonIso(-10)
		{
			pt = p;
			eta = e;
			phi = ph;
			mass = m; 
		};  
		Jet(double p, double e, double ph, double m,
				double slft, double mlft) 
			:scale(0),photonIso(-10)
		{
			pt = p;    sptlft = slft;
			eta = e;   mglft = mlft;
			phi = ph;
			mass = m; 
		};  

		Jet(double* jetvec) {
			pt = jetvec[0];  phi = jetvec[2];
			eta = jetvec[1]; mass = jetvec[3];
		};

		Jet()
			: pt(1e-8),eta(1e-8),phi(1e-8),mass(1e-8),
			mglft(0), sptlft(0), scale(0),photonIso(-10)
	{ };
		~Jet(){ };

		Jet(const TLorentzVector& vec ){
			if(vec.Pt()>0)
			{
				pt = vec.Pt();eta = vec.Eta();phi = vec.Phi();mass = vec.M(); 
			}
			else 
			{
				pt = 1e-8;eta = 1e-8;phi = 1e-8;mass = 1e-8; 
			}
		};
		inline TLorentzVector  getLorentzVector(){
			TLorentzVector vec;
			vec.SetPtEtaPhiM(pt,eta,phi,mass);
			return (vec);
		};
}; 

//class deltaR
class DeltaRDistance
{
	public:
		double operator()(const Jet& jet1, const Jet& jet2) const
		{
			const double deltaEta = jet1.eta - jet2.eta;
			double deltaPhi = jet1.phi - jet2.phi;
			if (deltaPhi < -M_PI)
				deltaPhi += 2.0*M_PI;
			if (deltaPhi > M_PI)
				deltaPhi -= 2.0*M_PI;
			const double distance = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			return (distance);
		}
		double operator()(const GenParticle& parton, const Jet& jet) const
		{
			const double deltaEta = parton.eta - jet.eta;
			double deltaPhi = parton.phi - jet.phi;
			if (deltaPhi < -M_PI)
				deltaPhi += 2.0*M_PI;
			if (deltaPhi > M_PI)
				deltaPhi -= 2.0*M_PI;
			const double distance = sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
			return (distance);
		}
};
double DeltaPhiDistance (double phi1, double phi2)
{
	double deltaPhi = phi1-phi2;
	if (deltaPhi < -M_PI)
		deltaPhi += 2.0*M_PI;
	if (deltaPhi > M_PI)
		deltaPhi -= 2.0*M_PI;

	return deltaPhi;
};
//
void getVecFromJets(const vector<GenParticle>& genparticles, vector<Jet>& jets, 
		const vector<int>&  matchVector,const DeltaRDistance& distance,
		TLorentzVector* outputVec, double* maxDrgj, double* maxDrjj)
{
	bool verbose(false);
	const unsigned numDaus(genparticles.size());
	for(unsigned i(0);i<numDaus;++i)
	{
		GenParticle gp(genparticles[i]);
		double mtch(matchVector[i]);
		if(mtch>-1) 
		{
			double tmpdr(distance(gp,jets[mtch]));
			if(tmpdr>(*maxDrgj)) (*maxDrgj)=tmpdr;
			(*outputVec) += jets[mtch].getLorentzVector();
			//find the farthest W jet matched to a genparton
			for(unsigned j(0);j<numDaus;++j)
			{
				double mtchjj (matchVector[j]);
				if(mtchjj>-1) 
				{
					double tmpdrjj(jets[mtchjj].getLorentzVector().DeltaR
							(jets[mtch].getLorentzVector()));
					if(tmpdrjj>(*maxDrjj)) (*maxDrjj)=tmpdrjj;
				}
			}

			if(verbose) {
				cout<<"------------------\nnum Daugther:"<<i<<endl;
				cout<<"Gen Particle pt eta phi mass\n"
					<<genparticles[i].pt<<' '
					<<genparticles[i].eta<<' '
					<<genparticles[i].phi<<' '
					<<genparticles[i].mass<<' '
					<<outputVec->M()<<' '
					<<endl;
				cout<<"Matched Jet pt eta phi mass\n"
					<<jets[mtch].pt<<' '
					<<jets[mtch].eta<<' '
					<<jets[mtch].phi<<' '
					<<jets[mtch].mass<<' '
					<<outputVec->M()<<' '
					<<endl; }
		}
	}
};
//
template <class T1>
void findBestJetsWithInvMass (const double ptCut, const double maxdrCut,
		const double invariantMass,
		const std::vector<T1>& inputJets, 
		T1* jet1, T1* jet2, double* mtmp ,double* drjj,
		const bool verbose,	int* jetcheck = NULL
		)
{
	const unsigned numjets = inputJets.size();
	TLorentzVector mass_vec;
	TLorentzVector daus_vec[numjets];
	for(unsigned i(0);i<numjets;++i)
		if(inputJets[i].pt>1e-8)
			daus_vec[i].SetPtEtaPhiM(inputJets[i].pt,
					inputJets[i].eta,
					inputJets[i].phi,
					inputJets[i].mass);
		else daus_vec[i].SetPtEtaPhiM(1e-8,1e-8,1e-8,1e-8);

	//	for(unsigned i(0);i<numjets;++i)
	//		if(daus_vec[i].Pt()<=1e-8)
	//			cout<<daus_vec[i].Pt()<<endl;


	for(unsigned i(0);i<numjets;++i)
		for(unsigned j(0);j<numjets;++j)
		{
			if(i!=j && daus_vec[i].Pt()>ptCut && daus_vec[j].Pt()>ptCut)
			{
				mass_vec = daus_vec[i]+daus_vec[j];
				double mass = mass_vec.M();
				double dr = daus_vec[j].DeltaR(daus_vec[i]);
				if(verbose)
					cout<<"TotalJets "<<numjets<<" Mass "<<mass<<" dr "<<dr
						<<" jet 1 "<<daus_vec[i].Pt()<<" eta "<<daus_vec[i].Eta()
						<<" jet 2 "<<daus_vec[j].Pt()<<" eta "<<daus_vec[j].Eta()
						//					<<" fabs(mass-invariantMass)<*mtmp "<<*mtmp
						//					<<" 0.<dr && dr<maxdrCut "<<(0.<dr && dr<maxdrCut)
						<<endl;

				if(fabs(mass-invariantMass)<*mtmp && 0.<dr && dr<maxdrCut)
				{
					*mtmp = fabs(mass-invariantMass);
					*drjj = dr;
					*jet1 = inputJets[i];
					*jet2 = inputJets[j];
					if(verbose)
						cout<<"\nPicked jet, total jets "<<numjets<<", mtmp "<<' '<<*mtmp
							<<" mass "<<mass
							<<endl;
					if(jetcheck) {
						jetcheck[0]=i;jetcheck[1]=j;jetcheck[2]=numjets;//to check jet selection
					}
				}
			}
		}
}
//
void analyzeJetSelection(const int* jetcheck,
		vector<unsigned>& analyzeVector
		)
{
	if(jetcheck[2]==2)
		analyzeVector[0]+=1;
	if(jetcheck[2]==3)
	{
		if(jetcheck[0]==0&&jetcheck[1]==1)
			analyzeVector[1]+=1;
		if(jetcheck[0]==0&&jetcheck[1]==2)
			analyzeVector[2]+=1;
		if(jetcheck[0]==1&&jetcheck[1]==2)
			analyzeVector[3]+=1;
	}//numjets3 fla3
	if(jetcheck[2]==4)
	{
		if(jetcheck[0]==0&&jetcheck[1]==1)
			analyzeVector[4]+=1;

		if(jetcheck[0]==0&&jetcheck[2]==2)
			analyzeVector[5]+=1;

		if(jetcheck[0]==0&&jetcheck[3]==4)
			analyzeVector[6]+=1;

		if(jetcheck[1]==1&&jetcheck[2]==2)
			analyzeVector[7]+=1;

		if(jetcheck[1]==2&&jetcheck[3]==3)
			analyzeVector[8]+=1;

		if(jetcheck[2]==3&&jetcheck[3]==4)
			analyzeVector[9]+=1;
	}
}

#endif // #ifdef VplusGTreeAnalyzer_cxx
