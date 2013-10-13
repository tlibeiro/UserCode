#include <algorithm>
#include <vector>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//fftheaders
#include "DataFormats/JetReco/interface/FFTPFJetCollection.h"
#include "RecoJets/FFTJetProducers/interface/FFTJetParameterParser.h"
#include "RecoJets/FFTJetAlgorithms/interface/clusteringTreeConverters.h"
#include "peakLifetime.hh"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTree.h"
#include "FFTJetFiller.h"
using namespace std;
//definitions
FFTJetFiller::FFTJetFiller(const std::string& n,
		TTree* t)
:name(n)
{
	assert(t);
	tree = t;
	isGen = (name.find("Gen")!=std::string::npos);
	for(unsigned i(0);i<maxJets;++i)
	{
		jetpt[i]=-1;   jeteta[i]=-1;
		jetphi[i]=-1;  jetmass[i]=-1;
		sp_lft[i]=-1;  mg_lft[i]=-1;
		pk_scl[i]=-1;
	}
	SetBranches();
}
//
void FFTJetFiller::SetBranches()
{
	SetBranch(jetpt,  (name+"_pt").c_str(),maxJets);
	SetBranch(jeteta, (name+"_eta").c_str(),maxJets);
	SetBranch(jetphi, (name+"_phi").c_str(),maxJets);
	SetBranch(jetmass,(name+"_m").c_str(),maxJets);
	SetBranch(sp_lft, (name+"_splft").c_str(),maxJets);
	SetBranch(mg_lft, (name+"_mglft").c_str(),maxJets);
	SetBranch(pk_scl, (name+"_scale").c_str(),maxJets);
}
//
void FFTJetFiller::fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag,
		const edm::ParameterSet& iConfig)
{
	//import tree helpers
	const edm::InputTag& treeTag(iConfig.getParameter<edm::InputTag>("treeLabel"));
	double completeEventScale = iConfig.getParameter<double>("completeEventScale");
	const edm::ParameterSet& initialscales_ (iConfig.getParameter<edm::ParameterSet>("InitialScales"));
	auto_ptr<std::vector<double> >  iniScales = fftjetcms::fftjet_ScaleSet_parser(initialscales_);
	checkConfig(iniScales, "invalid set of scales");
	std::sort(iniScales->begin(), iniScales->end(), std::greater<double>());
	//import clustering tree
	typedef reco::PattRecoTree<float,reco::PattRecoPeak<float> > StoredTree;
	SparseTree sparseTree;
	edm::Handle<StoredTree> input;
	iEvent.getByLabel(treeTag,input);
	fftjetcms::sparsePeakTreeFromStorable(*input, iniScales.get(),completeEventScale,&sparseTree);
	sparseTree.sortNodes();
	const double minScale(sparseTree.minScale()),
				maxScale(sparseTree.maxScale());
	iEvent.getByLabel(inputTag,jetCollection);  
	//
	const unsigned numJets(jetCollection->size());
	for(unsigned i(0);i<numJets && i<maxJets;++i)
	{
		reco::FFTAnyJet<reco::PFJet> jet(jetCollection->at(i));
		jetpt[i]=jet.pt(); 
		jeteta[i]=jet.eta(); 
		jetphi[i]=jet.phi(); 
		jetmass[i]=jet.mass(); 
		pk_scl[i]=jet.getFFTSpecific().f_precluster().scale();
		sp_lft[i]=peakSplitTime(sparseTree, jet.getFFTSpecific().f_code(),minScale);
		mg_lft[i]=peakMergeTime(sparseTree, jet.getFFTSpecific().f_code(),maxScale);
	}
}
/////// Helper functions ////////////////////////////////
//////////////////////////////////////////////////////////////////
	template<class Ptr>
void checkConfig(const Ptr& ptr, const char* message)
{
	if (ptr.get() == NULL)
		throw cms::Exception("FFTJetBadConfig") << message << std::endl;
}
//
void FFTJetFiller::SetBranch(float* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/F";
	tree->Branch( name, x, oss.str().c_str() );
}
//
void FFTJetFiller::SetBranch(int* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/I";
	tree->Branch( name, x, oss.str().c_str() );
}
