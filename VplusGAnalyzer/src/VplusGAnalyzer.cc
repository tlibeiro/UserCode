// -*- C++ -*-
//
// Package:    VplusGAnalyzer
// Class:      VplusGAnalyzer
// 
/**\class VplusGAnalyzer VplusGAnalyzer.cc TopQuarkAnalysis/VplusGAnalyzer/src/VplusGAnalyzer.cc
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
#include "VplusGAnalyzer.h"
#include "PhotonFiller.h"
#include "GroomedJetFiller.h"
#include "FFTJetFiller.h"
using namespace std;
using namespace reco;
//
// class declaration
//
class VplusGAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VplusGAnalyzer(const edm::ParameterSet&);
      ~VplusGAnalyzer();
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
			const edm::InputTag srcPhoton, 
						srcFatJet,
						srcGenFatJet,
						srcFFTJetR20,srcFFTJetR30,srcFFTJetR40; 
			const edm::InputTag srcPrimaryVertex;
			const std::string jetsForRho;
			const std::string JEC_GlobalTag;
			const edm::InputTag jetsForHT;
			const bool runOnMC;
      const edm::ParameterSet iC;
			TTree* outTree;
			double jetpt;
			double photonpt;
			int run  ; 
			int event;
			int lumi ;
			int bunch;
      float HT;
      const static unsigned maxGenParticles = 10;
      float genP_pt[maxGenParticles];
      float genP_eta[maxGenParticles];
      float genP_phi[maxGenParticles];
      float genP_m[maxGenParticles];
      float genP_id[maxGenParticles];
      float genP_motherid[maxGenParticles];
			auto_ptr<GroomedJetFiller> fatjetfiller;
			auto_ptr<GroomedJetFiller> genfatjetfiller;
			auto_ptr<PhotonFiller<vector<reco::Photon> > >  photonfiller;
			auto_ptr<FFTJetFiller>fftjetr20filler;
			auto_ptr<FFTJetFiller>fftjetr30filler;
			auto_ptr<FFTJetFiller>fftjetr40filler;
};
VplusGAnalyzer::VplusGAnalyzer(const edm::ParameterSet& iConfig)
	: 
		srcPhoton(iConfig.getParameter<edm::InputTag>("srcPhoton")),
		srcFatJet(iConfig.getParameter<edm::InputTag>("srcFatJet")),
		srcGenFatJet(iConfig.getParameter<edm::InputTag>("srcGenFatJet")),
		srcFFTJetR20(iConfig.getParameter<edm::InputTag>("srcFFTJetR20")),
		srcFFTJetR30(iConfig.getParameter<edm::InputTag>("srcFFTJetR30")),
		srcFFTJetR40(iConfig.getParameter<edm::InputTag>("srcFFTJetR40")),
		srcPrimaryVertex(iConfig.getParameter<edm::InputTag>("srcPrimaryVertex")),
		jetsForRho(iConfig.getParameter<string>("jetsForRho")),
		JEC_GlobalTag(iConfig.getParameter<string>("JEC_GlobalTag")),
		jetsForHT(iConfig.getParameter<edm::InputTag>("jetsForHT")),
		runOnMC(iConfig.getParameter<bool>("runOnMC")),
		iC(iConfig),
		outTree(0)
{
}
VplusGAnalyzer::~VplusGAnalyzer()
{
}
//
// member functions
//
// ------------ method called for each event  ------------
	void
VplusGAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
//	Handle<vector<reco::Photon> > photons;
//	Handle<vector<pat::Jet> > fatJets;
//  iEvent.getByLabel(srcPhoton,photons);
//	iEvent.getByLabel(srcFatJet,fatJets);
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
	if(runOnMC) 
	{
		iEvent.getByLabel(jetsForHT,genParticlesForHT);
		const unsigned numJetsHT (genParticlesForHT->size());
		const reco::Candidate *V=NULL;
		unsigned numParticlesSaved(0);
		HT=0;
		for(unsigned i(0);i<numJetsHT;++i) 
		{
			V = &((*genParticlesForHT)[i]);
			//    cout<<"Particle "<<i
			//				<<" id "<<V->pdgId()
			//				<<" status "<<V->status()
			//				<<" daus "<<V->numberOfDaughters()
			//				<<" mother: "<<(V->numberOfMothers()?V->mother()->pdgId():-1)
			//				<<endl;
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
		} // cycle one gen particles
	} // runonmc condition

	(fatjetfiller.get())->fillBranch(iEvent,srcFatJet,srcPrimaryVertex,jetsForRho);
	if(runOnMC)
		(genfatjetfiller.get())->fillBranch(iEvent,srcGenFatJet,srcPrimaryVertex,jetsForRho);
	//	//		(genfatjetfiller.get())->fillBranch(iEvent,srcGenFatJet,srcPrimaryVertex,jetsForRho,JEC_GlobalTag,runOnMC);
	(fftjetr20filler.get())->fillBranch(iEvent,srcFFTJetR20,iC);
	(fftjetr30filler.get())->fillBranch(iEvent,srcFFTJetR30,iC);
	(fftjetr40filler.get())->fillBranch(iEvent,srcFFTJetR40,iC);
	(photonfiller.get())->fillBranch(iEvent,srcPhoton,jetsForRho);
	outTree->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
	void 
VplusGAnalyzer::beginJob()
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
	//tree fillers
	fatjetfiller    = auto_ptr<GroomedJetFiller>(new GroomedJetFiller("GroomedJet_CA8",outTree,JEC_GlobalTag,runOnMC));
	genfatjetfiller = auto_ptr<GroomedJetFiller>(new GroomedJetFiller("GenGroomedJet_CA8",outTree,JEC_GlobalTag,runOnMC));
	fftjetr20filler = auto_ptr<FFTJetFiller>(new FFTJetFiller("FFTJetR20",outTree));
	fftjetr30filler = auto_ptr<FFTJetFiller>(new FFTJetFiller("FFTJetR30",outTree));
	fftjetr40filler = auto_ptr<FFTJetFiller>(new FFTJetFiller("FFTJetR40",outTree));
	photonfiller = auto_ptr<PhotonFiller<vector<reco::Photon> > >(new PhotonFiller<vector<reco::Photon> >("Photon",outTree,iC));
}
// ------------ method called once each job just after ending the event loop  ------------
	void 
VplusGAnalyzer::endJob() 
{
}
// ------------ method called when starting to processes a run  ------------
	void 
VplusGAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a run  ------------
	void 
VplusGAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
// ------------ method called when starting to processes a luminosity block  ------------
	void 
VplusGAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method called when ending the processing of a luminosity block  ------------
	void 
VplusGAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VplusGAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(VplusGAnalyzer);
