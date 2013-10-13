#define VplusGTreeAnalyzer_fftjetTest_cxx
#include "VplusGTreeAnalyzer_fftjetTest.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <cassert>
#include <iostream>
#include "matchOneToOne.hh"
using namespace std;

void BranchOn(TTree* fChain, double rValue, unsigned laValue);
void VplusGTreeAnalyzer_fftjetTest::Loop()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;

	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("event_runNo", 0); fChain->SetBranchStatus("event_lumi" , 0);
	fChain->SetBranchStatus("event_evtNo", 0); fChain->SetBranchStatus("event_HT"   , 0);
	fChain->SetBranchStatus("Photon_pt",1);
	fChain->SetBranchStatus("Photon_eta",1);
	fChain->SetBranchStatus("Photon_phi",1);
	fChain->SetBranchStatus("Photon_Id2012",1);
	//turn braches on 
	double rVals[]  = {0.2,0.3,0.4,0.5};
	unsigned laVals[] = {2,3,4};
	const unsigned numR(sizeof(rVals)/sizeof(rVals[0])), 
				numLA(sizeof(laVals)/sizeof(laVals[0]));
	for(unsigned i(0);i<numR;++i) 
		for(unsigned j(0);j<numLA;++j)
			BranchOn(fChain,rVals[i],laVals[j]); 

	fChain->SetBranchStatus("GenPart_pt",     1);
	fChain->SetBranchStatus("GenPart_eta",    1);
	fChain->SetBranchStatus("GenPart_phi",    1);
	fChain->SetBranchStatus("GenPart_m",      1);
	fChain->SetBranchStatus("GenPart_id",     1);
	fChain->SetBranchStatus("GenPart_id_mother",1);

	//../Set analysis parameters
	const unsigned totalJetsOrPhoton(4);
	const unsigned totalFFTJets(6);
	const unsigned totalGenParticles(10);
	const double ptCut(20), drCut(1.2),
				invariantWMass(80.3);
	vector<unsigned> analyzejetsel(10,0);
	vector<Jet> photons;
	Jet fftJetsla2[numR][totalFFTJets];
	Jet fftJetsla3[numR][totalFFTJets];
	Jet fftJetsla4[numR][totalFFTJets];
	DeltaRDistance distance;
	//Run over events
	//	  nentries = 10000;
  cout<<"Running over "<<nentries<<endl;
	for (Long64_t jentry=0; jentry<nentries;++jentry) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		//
		photons.clear();
		for(unsigned i(0);i<totalJetsOrPhoton;++i)
		{
			Jet photon(Photon_pt[i],
					Photon_eta[i],
					Photon_phi[i],
					0);
			photon.photonIso = Photon_Id2012[i];
			photons.push_back(photon);
		}
		vector<GenParticle> genparticles;
		vector<GenParticle> genphotons;
		for(unsigned i(0);i<totalFFTJets;++i) 
		{
			//la 2 
			fftJetsla2[0][i] = Jet(FFTJetLA2R20_pt[i],FFTJetLA2R20_eta[i],
					FFTJetLA2R20_phi[i],FFTJetLA2R20_m[i]);
			fftJetsla2[1][i] = Jet(FFTJetLA2R30_pt[i],FFTJetLA2R30_eta[i],
					FFTJetLA2R30_phi[i],FFTJetLA2R30_m[i]);
			fftJetsla2[2][i] = Jet(FFTJetLA2R40_pt[i],FFTJetLA2R40_eta[i],
					FFTJetLA2R40_phi[i],FFTJetLA2R40_m[i]);
			fftJetsla2[3][i] = Jet(FFTJetLA2R40_pt[i],FFTJetLA2R40_eta[i],
					FFTJetLA2R40_phi[i],FFTJetLA2R40_m[i]);
			//la 3             
			fftJetsla3[0][i] = Jet(FFTJetLA3R20_pt[i],FFTJetLA3R20_eta[i],
					FFTJetLA3R20_phi[i],FFTJetLA3R20_m[i]);
			fftJetsla3[1][i] = Jet(FFTJetLA3R30_pt[i],FFTJetLA3R30_eta[i],
					FFTJetLA3R30_phi[i],FFTJetLA3R30_m[i]);
			fftJetsla3[2][i] = Jet(FFTJetLA3R40_pt[i],FFTJetLA3R40_eta[i],
					FFTJetLA3R40_phi[i],FFTJetLA3R40_m[i]);
      fftJetsla3[3][i] = Jet(FFTJetLA3R40_pt[i],FFTJetLA3R40_eta[i],
					FFTJetLA3R40_phi[i],FFTJetLA3R40_m[i]);
			//la 4             
			fftJetsla4[0][i] = Jet(FFTJetLA4R20_pt[i],FFTJetLA4R20_eta[i],
					FFTJetLA4R20_phi[i],FFTJetLA4R20_m[i]);
			fftJetsla4[1][i] = Jet(FFTJetLA4R30_pt[i],FFTJetLA4R30_eta[i],
					FFTJetLA4R30_phi[i],FFTJetLA4R30_m[i]);
			fftJetsla4[2][i] = Jet(FFTJetLA4R40_pt[i],FFTJetLA4R40_eta[i],
					FFTJetLA4R40_phi[i],FFTJetLA4R40_m[i]);
 	    fftJetsla4[3][i] = Jet(FFTJetLA4R40_pt[i],FFTJetLA4R40_eta[i],
					FFTJetLA4R40_phi[i],FFTJetLA4R40_m[i]);
		}
		//GenParticle processing
		for(unsigned i(0);i<totalGenParticles;++i)
		{
			GenParticle gp(GenPart_pt[i],GenPart_eta[i],GenPart_phi[i],
					GenPart_m[i],GenPart_id[i],GenPart_id_mother[i]);
			// cout<<"particle:"<<gp.id<<" mohter:"<<gp.mother_id<<endl;
			if(abs(gp.mother_id)==24) 
				genparticles.push_back(gp);
			if(abs(gp.id)==22)
				genphotons.push_back(gp);
			//				cout<<"particle "<<i<<" id: "<<gp.id
			//				<<" mother id "<<gp.mother_id
			//				<<endl;
		}
		assert(genphotons.size()==1);
		const unsigned numWDaughters(genparticles.size());
		//get generated W 
		TLorentzVector vecGen;
		for(unsigned i(0);i<numWDaughters;++i)
		{
			GenParticle gp(genparticles[i]);
			TLorentzVector tmpvec;
			tmpvec.SetPtEtaPhiM(gp.pt,gp.eta,gp.phi,gp.mass);
			vecGen += tmpvec;
		} 
		//ntup fill vec
		vector<float> ntup_vec_radscan;
		for(unsigned i(0);i<numR;++i) {

			vector<Jet> fftVecla2,fftVecla3,fftVecla4;
			for(unsigned j(0);j<totalFFTJets;++j) {
				fftVecla2.push_back(fftJetsla2[i][j]);
				fftVecla3.push_back(fftJetsla3[i][j]);
				fftVecla4.push_back(fftJetsla4[i][j]);
			}
			//match jets to particles 
			vector<int> matchedJetsla2, matchedJetsla3,matchedJetsla4;
			matchOneToOne(genparticles,fftVecla2,DeltaRDistance(),&matchedJetsla2,1000.0);
			matchOneToOne(genparticles,fftVecla3,DeltaRDistance(),&matchedJetsla3,1000.0);
			matchOneToOne(genparticles,fftVecla4,DeltaRDistance(),&matchedJetsla4,1000.0);
			/// fill ntuples variables 
			double maxdr20(-1), maxdr30(-1), maxdr40(-1);
			double maxdrjj20(-1), maxdrjj30(-1), maxdrjj40(-1);
			bool allM20, allM30, allM40;
			TLorentzVector vec20,vec30,vec40;
			getVecFromJets(genparticles,fftVecla2,matchedJetsla2,distance,&vec20,&maxdr20,&maxdrjj20,&allM20);
			getVecFromJets(genparticles,fftVecla3,matchedJetsla3,distance,&vec30,&maxdr30,&maxdrjj30,&allM30);
			getVecFromJets(genparticles,fftVecla4,matchedJetsla4,distance,&vec40,&maxdr40,&maxdrjj40,&allM40);
			//start push back to vec
			ntup_vec_radscan.push_back(vec20.M());
			ntup_vec_radscan.push_back(maxdr20);
			ntup_vec_radscan.push_back(maxdrjj20);
			ntup_vec_radscan.push_back(allM20);
			ntup_vec_radscan.push_back(vec30.M());
			ntup_vec_radscan.push_back(maxdr30);
			ntup_vec_radscan.push_back(maxdrjj30);
			ntup_vec_radscan.push_back(allM30);
			ntup_vec_radscan.push_back(vec40.M());
			ntup_vec_radscan.push_back(maxdr40);
			ntup_vec_radscan.push_back(maxdrjj40);
			ntup_vec_radscan.push_back(allM40);
		}// ntupe fill vec mehtod 
		ntup_vec_radscan.push_back(numWDaughters?vecGen.DeltaR(genphotons[0].getLorentzVector()):-1);
		ntup_vec_radscan.push_back(vecGen.M());
		ntup_vec_radscan.push_back(vecGen.Pt());
		ntup_vec_radscan.push_back(numWDaughters);
		ntup_vec_radscan.push_back(photons[0].photonIso);
		assert(ntup_vec_radscan.size()==53);
		ntup_RadScan->Fill(&ntup_vec_radscan[0]);
		///
		//print out the event num
		if(!(jentry%20000))
			cout<<"Processing Event "<<jentry<<endl;
	}
	cout<<"Done\n";
	f_out = new TFile(outfile.c_str(),"RECREATE");
	f_out->cd();
	ntup_RadScan->Write();
	f_out->Write();
	f_out->Close();

}


void BranchOn(TTree* fChain, double rValue, unsigned laValue)
{
	ostringstream rval, laval, jettag;
	rval<<rValue*100;
	laval<<laValue;
	jettag<<"FFTJet"<<"LA"<<laval.str()<<"R"<<rval.str();

	fChain->SetBranchStatus((jettag.str()+"_pt").c_str(),   1); 
	fChain->SetBranchStatus((jettag.str()+"_eta").c_str(),  1); 
	fChain->SetBranchStatus((jettag.str()+"_phi").c_str(),  1); 
	fChain->SetBranchStatus((jettag.str()+"_m").c_str(),    1); 
	fChain->SetBranchStatus((jettag.str()+"_splft").c_str(),0); 
	fChain->SetBranchStatus((jettag.str()+"_mglft").c_str(),0); 
	fChain->SetBranchStatus((jettag.str()+"_scale").c_str(),0);
};
