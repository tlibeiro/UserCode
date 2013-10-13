#ifndef VPLUSGANALYZER_GROOMEDJETFILLER
#define VPLUSGANALYZER_GROOMEDJETFILLER
#include <fastjet/ClusterSequence.hh>
    //#include <fastjet/ActiveAreaSpec.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TLorentzVector.h"
#include "TopQuarkAnalysis/VplusGAnalyzer/src/NjettinessPlugin.hh"
#include "TopQuarkAnalysis/VplusGAnalyzer/src/Nsubjettiness.hh"
#include "TVector3.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include <iostream>
///filler classes for groomed jet 
class GroomedJetFiller {
	public: 
		static const unsigned maxJets=4;
		GroomedJetFiller(const std::string& name,
				TTree* t,const std::string& JEC_GlobalTag,const bool runMC);
		edm::Handle<std::vector<pat::Jet> > jetCollection;
		const std::string name;
    const bool runOnMC;
		TTree* tree;
		void fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag,	
										const edm::InputTag& srcPrimaryVertex, const std::string& jetsForRho);
		  float tau1   [maxJets];
			float tau2   [maxJets];
			float tau3   [maxJets];
			float tau4   [maxJets];
			float tau2tau1[maxJets];

			float jetpt [maxJets];      float jetpt_uncorr[maxJets];
			float jetm  [maxJets];      float jetm_uncorr[maxJets];
			float jeteta[maxJets];     
			float jetphi[maxJets];    
			float jete  [maxJets];
      float jetar [maxJets];

			float jetpt_ft  [maxJets];  float jetpt_pr  [maxJets];
			float jetphi_ft [maxJets];  float jetphi_pr [maxJets];
			float jeteta_ft [maxJets];  float jeteta_pr [maxJets];
			float jetm_ft   [maxJets];  float jetm_pr   [maxJets];
			float jete_ft   [maxJets];  float jete_pr   [maxJets];
			float jetar_ft  [maxJets];  float jetar_pr  [maxJets];

			float jetpt_tr  [maxJets];   float jetpt_pr_uncorr   [maxJets];
			float jetphi_tr [maxJets];   float jetm_pr_uncorr    [maxJets];
			float jeteta_tr [maxJets];   float jetpt_tr_uncorr   [maxJets];
			float jetm_tr   [maxJets];   float jetm_tr_uncorr    [maxJets];
			float jete_tr   [maxJets];   float jetpt_ft_uncorr   [maxJets];
			float jetar_tr  [maxJets];   float jetm_ft_uncorr    [maxJets];

			float prsubjet1_pt[maxJets];    float prsubjet2_pt[maxJets];
			float prsubjet1_eta[maxJets];   float prsubjet2_eta[maxJets];
			float prsubjet1_phi[maxJets];   float prsubjet2_phi[maxJets];
			float prsubjet1_m[maxJets];     float prsubjet2_m[maxJets];

			void SetBranches();
			void SetBranch(float* x,const char* name,unsigned maxJetsOrPhotons);
			void SetBranch(int* x  ,const char* name,unsigned maxJetsOrPhotons);
			bool isGen;
			double rhoVal_;
			unsigned nPV_;
			FactorizedJetCorrector* jec_;
			JetCorrectionUncertainty* jecUnc_;
			const TLorentzVector getCorrectedJet(fastjet::PseudoJet& jet);
			double getJEC(double curJetEta,double curJetPt,double curJetE,double curJetArea);
};
//definitions
GroomedJetFiller::GroomedJetFiller(const std::string& n,
		TTree* t,const std::string& JEC_GlobalTag,const bool runMC)
:name(n),runOnMC(runMC)
{
	assert(t);
	tree = t;
	isGen = (name.find("Gen")!=std::string::npos);
	for(unsigned i(0);i<maxJets;++i)
	{
		jetpt[i] =-1;     jetpt_uncorr[i]=-1; 
		jetm[i]  =-1;     jetm_uncorr[i]=-1; 
		jeteta[i]=-1;
		jetphi[i]=-1;
		jetar[i] =-1;
		jete[i]  =-1;

		jetpt_tr[i] =-1;  jetpt_pr_uncorr[i]=-1;
		jetphi_tr[i]=-1;  jetm_pr_uncorr[i]=-1;
		jeteta_tr[i]=-1;  jetpt_ft_uncorr[i]=-1;
		jetm_tr[i]  =-1;  jetm_ft_uncorr[i]=-1;  
		jetar_tr[i]=-1;   jetpt_tr_uncorr[i]=-1;  
		jete_tr[i]=-1;    jetm_tr_uncorr[i]=-1; 

		jetpt_ft[i] =-1;  jetpt_pr[i] =-1;
		jetphi_ft[i]=-1;  jetphi_pr[i]=-1;
		jeteta_ft[i]=-1;  jeteta_pr[i]=-1;
		jetm_ft[i]  =-1;  jetm_pr[i]  =-1;
		jete_ft[i]  =-1;  jete_pr[i]  =-1;
		jetar_ft[i]  =-1; jetar_pr[i]  =-1;

		tau1[i]=-1;  
		tau2[i]=-1;
		tau3[i]=-1;
		tau4[i]=-1;
		tau2tau1[i]=-1;

		prsubjet1_pt[i]=-1;   prsubjet2_pt[i]=-1;
		prsubjet1_eta[i]=-1;  prsubjet2_eta[i]=-1;
		prsubjet1_phi[i]=-1;  prsubjet2_phi[i]=-1;
		prsubjet1_m[i]=-1;    prsubjet2_m[i]=-1;
	}
	SetBranches();
	// ---- setting up the jec on-the-fly from text files...
	const std::string& fDir = JEC_GlobalTag;
	std::vector< JetCorrectorParameters > jecPars;
	std::vector< std::string > jecStr;
	const bool applyJECToGroomedJets(true);

	if(applyJECToGroomedJets) {
		jecStr.push_back( fDir + "_L1FastJet_AK7PFchs.txt" );
		jecStr.push_back( fDir + "_L2Relative_AK7PFchs.txt" );
		jecStr.push_back( fDir + "_L3Absolute_AK7PFchs.txt" );
		if (!runOnMC)
			jecStr.push_back( fDir + "_L2L3Residual_AK7PFchs.txt" );

		for (unsigned int i = 0; i < jecStr.size(); ++i){
			JetCorrectorParameters* ijec = new JetCorrectorParameters( jecStr[i] );
			jecPars.push_back( *ijec );
		}
		jec_ = new FactorizedJetCorrector(jecPars);
		jecUnc_ = new JetCorrectionUncertainty( fDir + "_Uncertainty_AK7PFchs.txt" );
	}
	//end corrections
}
//
void GroomedJetFiller::SetBranches()
{
	//set branches
	SetBranch(jetpt    ,(name+"_pt").c_str(),maxJets);
	SetBranch(jeteta   ,(name+"_eta").c_str(),maxJets);
	SetBranch(jetphi   ,(name+"_phi").c_str(),maxJets);
	SetBranch(jetm     ,(name+"_m").c_str(),maxJets);
	SetBranch(jete     ,(name+"_e").c_str(),maxJets);
	SetBranch(jetar    ,(name+"_ar").c_str(),maxJets);
	SetBranch(jetpt_uncorr,(name+"_pt_uncorr").c_str(),maxJets);
	SetBranch(jetm_uncorr ,(name+"_m_uncorr").c_str(),maxJets);
	//
	SetBranch(jetpt_tr ,(name+"_pt_tr").c_str(),maxJets);
	SetBranch(jeteta_tr,(name+"_eta_tr").c_str(),maxJets);
	SetBranch(jetphi_tr,(name+"_phi_tr").c_str(),maxJets);
	SetBranch(jetm_tr  ,(name+"_m_tr").c_str(),maxJets);
	SetBranch(jete_tr  ,(name+"_e_tr").c_str(),maxJets);
	SetBranch(jetar_tr  ,(name+"_ar_tr").c_str(),maxJets);
	//
	SetBranch(jetpt_ft ,(name+"_pt_ft").c_str(),maxJets);
	SetBranch(jeteta_ft,(name+"_eta_ft").c_str(),maxJets);
	SetBranch(jetphi_ft,(name+"_phi_ft").c_str(),maxJets);
	SetBranch(jetm_ft  ,(name+"_m_ft").c_str(),maxJets);
	SetBranch(jete_ft  ,(name+"_e_ft").c_str(),maxJets);
	SetBranch(jetar_ft  ,(name+"_ar_ft").c_str(),maxJets);
	//
	SetBranch(jetpt_pr ,(name+"_pt_pr").c_str(),maxJets);
	SetBranch(jeteta_pr,(name+"_eta_pr").c_str(),maxJets);
	SetBranch(jetphi_pr,(name+"_phi_pr").c_str(),maxJets);
	SetBranch(jetm_pr  ,(name+"_m_pr").c_str(),maxJets);
	SetBranch(jete_pr  ,(name+"_e_pr").c_str(),maxJets);
	SetBranch(jetar_pr  ,(name+"_ar_pr").c_str(),maxJets);
	//
	SetBranch(tau1     ,(name+"_tau1").c_str(),maxJets);
	SetBranch(tau2     ,(name+"_tau2").c_str(),maxJets);
	SetBranch(tau3     ,(name+"_tau3").c_str(),maxJets);
	SetBranch(tau4     ,(name+"_tau4").c_str(),maxJets);
	SetBranch(tau2tau1 ,(name+"_tau2tau1").c_str(),maxJets);
	//
	if(!isGen) {
		SetBranch(jetpt_ft_uncorr ,(name+"_pt_ft_uncorr").c_str(),maxJets);
		SetBranch(jetm_ft_uncorr  ,(name+"_m_ft_uncorr").c_str(),maxJets);
		SetBranch(jetpt_tr_uncorr ,(name+"_pt_tr_uncorr").c_str(),maxJets);
		SetBranch(jetm_tr_uncorr  ,(name+"_m_tr_uncorr").c_str(),maxJets);
		SetBranch(jetpt_pr_uncorr ,(name+"_pt_pr_uncorr").c_str(),maxJets);
		SetBranch(jetm_pr_uncorr  ,(name+"_m_pr_uncorr").c_str(),maxJets);
	}
	//
	SetBranch(prsubjet1_pt ,(name+"_pt_prsubj1").c_str(),maxJets);
	SetBranch(prsubjet1_eta,(name+"_eta_prsubj1").c_str(),maxJets);
	SetBranch(prsubjet1_phi,(name+"_phi_prsubj1").c_str(),maxJets);
	SetBranch(prsubjet1_m  ,(name+"_m_prsubj1").c_str(),maxJets);
	SetBranch(prsubjet2_pt ,(name+"_pt_prsubj2").c_str(),maxJets);
	SetBranch(prsubjet2_eta,(name+"_eta_prsubj2").c_str(),maxJets);
	SetBranch(prsubjet2_phi,(name+"_phi_prsubj2").c_str(),maxJets);
	SetBranch(prsubjet2_m  ,(name+"_m_prsubj2").c_str(),maxJets);
}
//
void GroomedJetFiller::fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag, 
		const edm::InputTag& srcPrimaryVertex, const std::string& jetsForRho) 
{
	//    edm::Handle<vector<reco::GenJet> > genJets;
	edm::Handle<reco::GenParticleRefVector> genParticles;
	reco::GenParticleRefVector genParticlesforJets;
	if (isGen) {
		iEvent.getByLabel(inputTag, genParticles);
		genParticlesforJets = *genParticles;
	}
	else
		iEvent.getByLabel(inputTag,jetCollection); 

	// ------ get rho --------
	rhoVal_ = -99.;
	edm::Handle<double> rho;
	const edm::InputTag eventrho(jetsForRho, "rho");
	iEvent.getByLabel(eventrho,rho);
	rhoVal_ = *rho;
	// ------ get nPV: primary/secondary vertices------
	nPV_ = 0.;
	double nPVval = 0;
	edm::Handle <edm::View<reco::Vertex> > recVtxs;
	iEvent.getByLabel(srcPrimaryVertex,recVtxs);
	for(unsigned int ind=0;ind<recVtxs->size();ind++){
		if (!((*recVtxs)[ind].isFake()) && ((*recVtxs)[ind].ndof()>=4)
				&& (fabs((*recVtxs)[ind].z())<=24.0) &&
				((*recVtxs)[ind].position().Rho()<=2.0) ) {
			nPVval += 1;
		}
	}
	nPV_ = nPVval;
	//get pf candidates for jets   
	std::vector<fastjet::PseudoJet> FJparticles;
	if(isGen) 
	{
		const unsigned numJetCandidates(genParticles->size());
		for(unsigned i(0);i<numJetCandidates;++i)
		{
			const reco::GenParticle& P = *(genParticlesforJets[i]);
			FJparticles.push_back(fastjet::PseudoJet(P.px(),
						P.py(),P.pz(),P.energy()) );
		}
	}
	else {
		const unsigned numFatJets(jetCollection->size());
		for(unsigned j(0);j<numFatJets;++j) { 
			std::vector <reco::PFCandidatePtr> jetCandidates  (jetCollection->at(j).getPFConstituents());
			const unsigned numJetCandidates(jetCandidates.size());
			for(unsigned i(0);i<numJetCandidates;++i)
			{
				FJparticles.push_back(fastjet::PseudoJet(jetCandidates.at(i)->px(),
							jetCandidates.at(i)->py(),
							jetCandidates.at(i)->pz(),
							jetCandidates.at(i)->energy() ) );
			}//make candidates 
		}//cycle over jets 
	}
	//// //test genjets 
	// edm::Handle<vector<reco::GenJet> > gjtest;
	//  iEvent.getByLabel("ca8GenJetsNoNu",gjtest);
	//  const unsigned numgj = gjtest->size();
	//  for(unsigned i(0);i<numgj;++i)
	//   {
	//		const reco::GenJet& gps = gjtest->at(i);
	//    cout<<"gps "<<gps.getGenConstituent(0)->pt()<<endl;
	//	 }
	// do re-clustering
	fastjet::JetDefinition jetDef(fastjet::cambridge_algorithm, 0.8);
	jetDef.set_jet_algorithm( fastjet::cambridge_algorithm );

	int activeAreaRepeats = 1;
	double ghostArea = 0.01;
	double ghostEtaMax = 5.0;
	// fastjet::ActiveAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
	fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
	fjActiveArea.set_fj2_placement(true);
	fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area_explicit_ghosts, fjActiveArea );
	fastjet::ClusterSequenceArea thisClustering(FJparticles, jetDef, fjAreaDefinition);

	std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering.inclusive_jets(50.0));
	fastjet::ClusterSequence thisClustering_basic(FJparticles, jetDef);
	std::vector<fastjet::PseudoJet> out_jets_basic = sorted_by_pt(thisClustering_basic.inclusive_jets(50.0));
	// define groomers
	fastjet::Filter trimmer(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), 
				fastjet::SelectorPtFractionMin(0.03)));
	fastjet::Filter filter(fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm, 0.3), 
				fastjet::SelectorNHardest(3)));
	fastjet::Pruner pruner(fastjet::cambridge_algorithm, 0.1, 0.5);

	std::vector<fastjet::Transformer const *> transformers;
	transformers.push_back(&trimmer);
	transformers.push_back(&filter);
	transformers.push_back(&pruner);
	// Defining Nsubjettiness parameters
	// power for angular dependence, 
	//e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
	double beta = 1; 
	double R0 = 0.8; 	// Characteristic jet radius for normalization
	double Rcut = 0.8; 	// maximum R particles can be from axis to be included in jet

	// s t a r t   l o o p   o n   j e t s
	const unsigned numJets(out_jets.size());
	for(unsigned i(0);i<maxJets && i<numJets;++i)
	{
		// save ungroomed jet
		jetm_uncorr[i] = out_jets.at(i).m();
		jetpt_uncorr[i] = out_jets.at(i).pt();
		const TLorentzVector& jet_corr = getCorrectedJet(out_jets.at(i));
		jetm[i] = jet_corr.M();
		jetpt[i] = jet_corr.Pt();
		jeteta[i] = jet_corr.Eta();
		jetphi[i] = jet_corr.Phi();
		jete[i]   = jet_corr.Energy();
		jetar[i] = out_jets.at(i).area();
		// pruning, trimming, filtering  -------------
		int transctr = 0;
		for(std::vector<fastjet::Transformer const *>::const_iterator
				itransf = transformers.begin(), itransfEnd = transformers.end();
				itransf != itransfEnd; ++itransf)
		{
			fastjet::PseudoJet transformedJet = out_jets.at(i);
			transformedJet = (**itransf)(transformedJet);
			fastjet::PseudoJet transformedJet_basic = out_jets_basic.at(i);
			transformedJet_basic = (**itransf)(transformedJet_basic);

			if(transctr == 0) {
				jetm_tr_uncorr[i] = transformedJet.m();
				jetpt_tr_uncorr[i] = transformedJet.pt();
				const TLorentzVector& jet_tr_corr = getCorrectedJet(transformedJet);
				jetm_tr[i] = jet_tr_corr.M();
				jetpt_tr[i] = jet_tr_corr.Pt();
				jeteta_tr[i] = jet_tr_corr.Eta();
				jetphi_tr[i] = jet_tr_corr.Phi();
				jete_tr[i]   = jet_tr_corr.Energy();
				jetar_tr[i] = transformedJet.area();
			}
			else if(transctr == 1) {
				jetm_ft_uncorr[i] = transformedJet.m();
				jetpt_ft_uncorr[i] = transformedJet.pt();
				const TLorentzVector& jet_ft_corr = getCorrectedJet(transformedJet);
				jetm_ft[i] = jet_ft_corr.M();
				jetpt_ft[i] = jet_ft_corr.Pt();
				jeteta_ft[i] = jet_ft_corr.Eta();
				jetphi_ft[i] = jet_ft_corr.Phi();
				jete_ft[i]   = jet_ft_corr.Energy();
				jetar_ft[i] = transformedJet.area();
			}
			else if(transctr == 2) {
				jetm_pr_uncorr[i] = transformedJet.m();
				jetpt_pr_uncorr[i] = transformedJet.pt();
				const TLorentzVector& jet_pr_corr = getCorrectedJet(transformedJet);
				jetm_pr[i] = jet_pr_corr.M();
				jetpt_pr[i] = jet_pr_corr.Pt();
				jeteta_pr[i] = jet_pr_corr.Eta();
				jetphi_pr[i] = jet_pr_corr.Phi();
				jete_pr[i]   = jet_pr_corr.Energy();
				jetar_pr[i] = transformedJet.area();
				//decompose into requested number of subjets:
				if (transformedJet_basic.constituents().size() > 1){
					int nsubjetstokeep = 2;
					std::vector<fastjet::PseudoJet> subjets = transformedJet_basic.associated_cluster_sequence()->exclusive_subjets(transformedJet_basic,nsubjetstokeep);
					TLorentzVector sj1( subjets.at(0).px(),subjets.at(0).py(),subjets.at(0).pz(),subjets.at(0).e());
					TLorentzVector sj2( subjets.at(1).px(),subjets.at(1).py(),subjets.at(1).pz(),subjets.at(1).e());
					prsubjet1_pt[i] = subjets.at(0).pt();  prsubjet1_eta[i] = subjets.at(0).eta(); 
					prsubjet2_pt[i] = subjets.at(1).pt();  prsubjet2_eta[i] = subjets.at(1).eta(); 
					prsubjet1_phi[i] = subjets.at(0).phi(); prsubjet1_m[i] = subjets.at(0).m();
					prsubjet2_phi[i] = subjets.at(1).phi(); prsubjet2_m[i] = subjets.at(1).m();
				}
			}
			else
				std::cout << "error in number of transformers" << std::endl;
			++transctr;
		}
		fastjet::Nsubjettiness nSub1KT(1, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		tau1[i] = nSub1KT(out_jets.at(i));
		fastjet::Nsubjettiness nSub2KT(2, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		tau2[i] = nSub2KT(out_jets.at(i));
		fastjet::Nsubjettiness nSub3KT(3, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		tau3[i] = nSub3KT(out_jets.at(i));
		fastjet::Nsubjettiness nSub4KT(4, Njettiness::onepass_kt_axes, beta, R0, Rcut);
		tau4[i] = nSub4KT(out_jets.at(i));
		tau2tau1[i] = tau2[i]/tau1[i];
	}//jet loop
}
/////////Helper functions //////////////
void GroomedJetFiller::SetBranch(float* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/F";
	tree->Branch( name, x, oss.str().c_str() );
}
//
void GroomedJetFiller::SetBranch(int* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/I";
	tree->Branch( name, x, oss.str().c_str() );
}
//
const TLorentzVector GroomedJetFiller::getCorrectedJet(fastjet::PseudoJet& jet) {

	double jecVal = 1.0;

	if(!isGen)
		jecVal = getJEC( jet.eta(), jet.pt(), jet.e(), jet.area());

	return(TLorentzVector (jet.px() * jecVal,
				jet.py() * jecVal,
				jet.pz() * jecVal,
				jet.e() * jecVal));

}
//
double GroomedJetFiller::getJEC(double curJetEta,
		double curJetPt,
		double curJetE,
		double curJetArea){

	jec_->setJetEta( curJetEta );
	jec_->setJetPt ( curJetPt );
	jec_->setJetE  ( curJetE );
	jec_->setJetA  ( curJetArea );
	jec_->setRho   ( rhoVal_ );
	jec_->setNPV   ( nPV_ );
	double corr = jec_->getCorrection();
	return corr;
}
#endif
