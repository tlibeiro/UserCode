#ifndef VPLUSGANALYZER_FFTJETFILLER_
#define VPLUSGANALYZER_FFTJETFILLER_

#include <algorithm>
#include <vector>
#include <iostream>
#include "RecoJets/FFTJetProducers/interface/FFTJetParameterParser.h"
// functions which manipulate storable trees
#include "RecoJets/FFTJetAlgorithms/interface/clusteringTreeConverters.h"
#include "peakLifetime.hh"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
template<class Ptr>
void checkConfig(const Ptr& ptr, const char* message);
///filler classes for groomed jet 
class FFTJetFiller {
	public: 
		static const unsigned maxJets=6;
		FFTJetFiller(const std::string& name,
				TTree* t);
		edm::Handle<reco::FFTPFJetCollection> jetCollection;
		const std::string name;
		TTree* tree;
		void fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag,	
				const edm::ParameterSet& iConfig);
		void SetBranches();
		void SetBranch(float* x, const char* name,unsigned maxJetsOrPhotons);
		void SetBranch(int* x  , const char* name,unsigned maxJetsOrPhotons);

		bool isGen;

		float jetpt[maxJets];   float sp_lft[maxJets];
		float jeteta[maxJets];  float mg_lft[maxJets];
		float jetphi[maxJets];  float pk_scl[maxJets];
		float jetmass[maxJets];
};


#endif
