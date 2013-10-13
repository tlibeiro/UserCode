// -*- C++ -*-
//
// Package:    VplusGAnalyzer_fftjetTest
// Class:      VplusGAnalyzer_fftjetTest
// 
/**\class VplusGAnalyzer_fftjetTest VplusGAnalyzer_fftjetTest.cc TopQuarkAnalysis/VplusGAnalyzer_fftjetTest/src/VplusGAnalyzer_fftjetTest.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Terence Libeiro
//         Created:  Mon Jul 29 13:47:06 CDT 2013
// $Id$
//
// system include files
#include <memory>
#include <algorithm>
#include <vector>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/FFTPFJetCollection.h"
#include "RecoJets/FFTJetProducers/interface/FFTJetParameterParser.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
// root headers
#include "TTree.h"
#include "FFTJetFiller.h"
#include "PhotonFiller.h"
using namespace std;
using namespace reco;
//
// class declaration
//
class VplusGAnalyzer_fftjetTest : public edm::EDAnalyzer {
   public:
      explicit VplusGAnalyzer_fftjetTest(const edm::ParameterSet&);
      ~VplusGAnalyzer_fftjetTest();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
			virtual void beginRun(edm::Run const&, edm::EventSetup const&);
			virtual void endRun(edm::Run const&, edm::EventSetup const&);
			virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
			virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
			// ----------member data ---------------------------
			const edm::InputTag srcPhoton;
			vector<edm::InputTag> srcFFTJetLA2; 
			vector<edm::InputTag> srcFFTJetLA3; 
			vector<edm::InputTag> srcFFTJetLA4; 
      vector<double> rValues;
			const std::string jetsForRho;
			const edm::InputTag jetsForHT;
			const bool runOnMC;
      const edm::ParameterSet iC;
			TTree* outTree;
			int run  ; 
			int event;
			int lumi ;
			int bunch;
      float HT;
      const static unsigned numRValues =4;
      const static unsigned maxGenParticles = 10;
      float genP_pt[maxGenParticles];
      float genP_eta[maxGenParticles];
      float genP_phi[maxGenParticles];
      float genP_m[maxGenParticles];
      float genP_id[maxGenParticles];
      float genP_motherid[maxGenParticles];
			auto_ptr<PhotonFiller<vector<reco::Photon> > >  photonfiller;
			auto_ptr<FFTJetFiller>fftjetfillerla2[numRValues];
			auto_ptr<FFTJetFiller>fftjetfillerla3[numRValues];
			auto_ptr<FFTJetFiller>fftjetfillerla4[numRValues];
};
VplusGAnalyzer_fftjetTest::VplusGAnalyzer_fftjetTest(const edm::ParameterSet& iConfig)
	: 
		srcPhoton(iConfig.getParameter<edm::InputTag>("srcPhoton")),
		srcFFTJetLA2(iConfig.getParameter<vector<edm::InputTag> >("srcFFTJetLA2")),
		srcFFTJetLA3(iConfig.getParameter<vector<edm::InputTag> >("srcFFTJetLA3")),
		srcFFTJetLA4(iConfig.getParameter<vector<edm::InputTag> >("srcFFTJetLA4")),
    rValues(iConfig.getParameter<vector<double> >("rValues")),
		jetsForRho(iConfig.getParameter<string>("jetsForRho")),
		jetsForHT(iConfig.getParameter<edm::InputTag>("jetsForHT")),
		runOnMC(iConfig.getParameter<bool>("runOnMC")),
		iC(iConfig),
		outTree(0)
{
	assert(numRValues==srcFFTJetLA2.size());
	assert(numRValues==srcFFTJetLA3.size());
	assert(numRValues==srcFFTJetLA4.size());
	assert(numRValues==rValues.size());
}
VplusGAnalyzer_fftjetTest::~VplusGAnalyzer_fftjetTest()
{
}
//
// member functions
//
// ------------ method called for each event  ------------
	void
VplusGAnalyzer_fftjetTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	//Write the tree
	// write event information: run, event, bunch crossing, ....
	run   = iEvent.id().run();
	event = iEvent.id().event();
	lumi  = iEvent.luminosityBlock();
	bunch = iEvent.bunchCrossing();

	for(unsigned i(0);i<maxGenParticles;++i)
	{
		genP_pt[i]=-1;    genP_eta[i]=-1;
		genP_phi[i]=-1;   genP_m[i]=-1;
		genP_id[i]=-10000;    genP_motherid[i]=-1;
	}
	Handle<vector<reco::GenParticle> > genParticlesForHT;
	iEvent.getByLabel(jetsForHT,genParticlesForHT);
	const unsigned numJetsHT (genParticlesForHT->size());
	const reco::Candidate *V=NULL;
	unsigned numParticlesSaved(0);
	HT=0;
	for(unsigned i(0);i<numJetsHT;++i) 
	{
		V = &((*genParticlesForHT)[i]);
		const int status(V->status()); 
		const int pdgid (V->pdgId()); 
		const int num_moth (V->numberOfMothers());
		const int mothid (num_moth?V->mother()->pdgId():-10000); 
		if(status==3 &&                    
				(abs(pdgid)<=6 || abs(pdgid)==21) //find gluons and quarks for HT
			)
			HT+=V->pt();
		if (status==3 &&                                             //stable Matrix element  part
				(abs(mothid) ==24 || abs(mothid)==23 || abs(pdgid)==22) && // photon, or W,Z daughters
				(numParticlesSaved<maxGenParticles)																		
			 )
		{	
			genP_pt[numParticlesSaved]=V->pt();
			genP_phi[numParticlesSaved]=V->phi();
			genP_eta[numParticlesSaved]=V->eta();
			genP_m[numParticlesSaved]=V->mass();
			genP_id[numParticlesSaved]=pdgid;
			genP_motherid[numParticlesSaved]=mothid;
			++numParticlesSaved;
		}
	}

	for(unsigned i(0);i<numRValues;++i) 
	{
   const InputTag LA2Label(srcFFTJetLA2[i]), 
									LA3Label(srcFFTJetLA3[i]), 
									LA4Label(srcFFTJetLA4[i]);
     
	(fftjetfillerla2[i].get())->fillBranch(iEvent,LA2Label,iC);
	(fftjetfillerla3[i].get())->fillBranch(iEvent,LA3Label,iC);
	(fftjetfillerla4[i].get())->fillBranch(iEvent,LA4Label,iC);
	}
	(photonfiller.get())->fillBranch(iEvent,srcPhoton,jetsForRho);
	outTree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
	void 
VplusGAnalyzer_fftjetTest::beginJob()
{
	edm::Service<TFileService> fs;
	ostringstream maxGenP; maxGenP<<maxGenParticles;
	outTree = fs->make<TTree>("VplusGTree","VplusGTree");
	outTree->Branch("event_runNo",  &run,   "event_runNo/I");
	outTree->Branch("event_evtNo",  &event, "event_evtNo/I");
	outTree->Branch("event_lumi",   &lumi,  "event_lumi/I");
	outTree->Branch("event_bunch",  &bunch, "event_bunch/I");
	outTree->Branch("event_HT",     &HT,    "event_HT/F");
	outTree->Branch("GenPart_pt",   genP_pt,("GenPart_pt["  +maxGenP.str()+"]/F").c_str());
	outTree->Branch("GenPart_eta",  genP_eta,("GenPart_eta["+maxGenP.str()+"]/F").c_str());
	outTree->Branch("GenPart_phi",  genP_phi,("GenPart_phi["+maxGenP.str()+"]/F").c_str());
	outTree->Branch("GenPart_m"  ,  genP_m,("GenPart_m["    +maxGenP.str()+"]/F").c_str());
	outTree->Branch("GenPart_id" ,  genP_id,("GenPart_id["  +maxGenP.str()+"]/F").c_str());
	outTree->Branch("GenPart_id_mother",genP_motherid,("GenPart_id_mother["+maxGenP.str()+"]/F").c_str());
	//	//tree fillers
	for(unsigned i(0);i<numRValues;++i)
   {
    ostringstream ossla2, ossla3, ossla4;
		ossla2<<"FFTJetLA2R"<<rValues[i]*100;
		ossla3<<"FFTJetLA3R"<<rValues[i]*100;
		ossla4<<"FFTJetLA4R"<<rValues[i]*100;
		fftjetfillerla2[i] = auto_ptr<FFTJetFiller>(new FFTJetFiller(ossla2.str(),outTree));
	  fftjetfillerla3[i] = auto_ptr<FFTJetFiller>(new FFTJetFiller(ossla3.str(),outTree));
	  fftjetfillerla4[i] = auto_ptr<FFTJetFiller>(new FFTJetFiller(ossla4.str(),outTree));
	}
	photonfiller = auto_ptr<PhotonFiller<vector<reco::Photon> > >(new PhotonFiller<vector<reco::Photon> >("Photon",outTree,iC));
}
// ------------ method called once each job just after ending the event loop  ------------
	void 
VplusGAnalyzer_fftjetTest::endJob() 
{
}
// ------------ method called when starting to processes a run  ------------
	void 
VplusGAnalyzer_fftjetTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a run  ------------
	void 
VplusGAnalyzer_fftjetTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
	void 
VplusGAnalyzer_fftjetTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a luminosity block  ------------
	void 
VplusGAnalyzer_fftjetTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VplusGAnalyzer_fftjetTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(VplusGAnalyzer_fftjetTest);
