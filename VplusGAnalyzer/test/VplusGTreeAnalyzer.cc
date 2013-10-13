#define VplusGTreeAnalyzer_cxx
#include "VplusGTreeAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <string>
#include <cassert>
#include <iostream>
#include "matchOneToOne.hh"
using namespace std;

void VplusGTreeAnalyzer::Loop()
{
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntries();
	Long64_t nbytes = 0, nb = 0;

	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("event_runNo", 1); fChain->SetBranchStatus("event_lumi" , 1);
	fChain->SetBranchStatus("event_evtNo", 1); fChain->SetBranchStatus("event_HT"   , 1);
	fChain->SetBranchStatus("GroomedJet_CA8_pt_pr" ,1);  fChain->SetBranchStatus("Photon_pt",1);
	fChain->SetBranchStatus("GroomedJet_CA8_eta_pr",1);  fChain->SetBranchStatus("Photon_eta",1);
	fChain->SetBranchStatus("GroomedJet_CA8_phi_pr",1);  fChain->SetBranchStatus("Photon_phi",1);
	fChain->SetBranchStatus("GroomedJet_CA8_m_pr"  ,1);  fChain->SetBranchStatus("Photon_Id2012",1);
	fChain->SetBranchStatus("GroomedJet_CA8_pt" ,1);     fChain->SetBranchStatus("GroomedJet_CA8_tau2tau1",1);
	fChain->SetBranchStatus("GroomedJet_CA8_eta",0);
	fChain->SetBranchStatus("GroomedJet_CA8_phi",0);
	fChain->SetBranchStatus("GroomedJet_CA8_m"  ,0);

	fChain->SetBranchStatus("FFTJetR20_pt",   1);   fChain->SetBranchStatus("FFTJetR30_pt",   1);
	fChain->SetBranchStatus("FFTJetR20_eta",  1);   fChain->SetBranchStatus("FFTJetR30_eta",  1);
	fChain->SetBranchStatus("FFTJetR20_phi",  1);   fChain->SetBranchStatus("FFTJetR30_phi",  1);
	fChain->SetBranchStatus("FFTJetR20_m",    1);   fChain->SetBranchStatus("FFTJetR30_m",    1);
	fChain->SetBranchStatus("FFTJetR20_splft",1);   fChain->SetBranchStatus("FFTJetR30_splft",1);
	fChain->SetBranchStatus("FFTJetR20_mglft",1);   fChain->SetBranchStatus("FFTJetR30_mglft",1);
	fChain->SetBranchStatus("FFTJetR20_scale",0);   fChain->SetBranchStatus("FFTJetR30_scale",0);

	fChain->SetBranchStatus("FFTJetR40_pt",   1);   fChain->SetBranchStatus("GenPart_pt",     1);
	fChain->SetBranchStatus("FFTJetR40_eta",  1);   fChain->SetBranchStatus("GenPart_eta",    1);
	fChain->SetBranchStatus("FFTJetR40_phi",  1);   fChain->SetBranchStatus("GenPart_phi",    1);
	fChain->SetBranchStatus("FFTJetR40_m",    1);   fChain->SetBranchStatus("GenPart_m",      1);
	fChain->SetBranchStatus("FFTJetR40_splft",1);   fChain->SetBranchStatus("GenPart_id",     1);
	fChain->SetBranchStatus("FFTJetR40_mglft",1);   fChain->SetBranchStatus("GenPart_id_mother",1);
	fChain->SetBranchStatus("FFTJetR40_scale",0);


  //../Set analysis parameters
	const unsigned totalJetsOrPhoton(4);
	const unsigned totalFFTJets(6);
  const unsigned totalGenParticles(5);
  const double ptCut(20), drCut(1.2),
							 invariantWMass(80.3);
  vector<unsigned> analyzejetsel(10,0);
	vector<float> ntup_vg_vec;
	vector<Jet> grJets, photons,fftJets;
  DeltaRDistance distance;
//Run over events
//	  nentries = 10000;
  cout<<"Processing Events "<<nentries<<endl;
	for (Long64_t jentry=0; jentry<nentries;++jentry) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		//
		grJets.clear();photons.clear();
		for(unsigned i(0);i<totalJetsOrPhoton;++i)
		{
			grJets.push_back(Jet(GroomedJet_CA8_pt_pr[i],
						GroomedJet_CA8_eta_pr[i],
						GroomedJet_CA8_phi_pr[i],
						GroomedJet_CA8_m_pr[i]
						));
      Jet photon(Photon_pt[i],
            Photon_eta[i],
            Photon_phi[i],
            0);
      photon.photonIso = Photon_Id2012[i];
			photons.push_back(photon);
		}
		fftJets.clear();
   vector<Jet> fftJets20,fftJets30,fftJets40;
   vector<GenParticle> genparticles;
   vector<GenParticle> genphotons;
	 for(unsigned i(0);i<totalFFTJets;++i) 
	 {
		 fftJets.push_back(Jet(FFTJetR30_pt[i],
					 FFTJetR30_eta[i],
					 FFTJetR30_phi[i],
					 FFTJetR30_m[i],				
					 FFTJetR30_splft[i],
					 FFTJetR30_mglft[i]
					 ));
		 fftJets20.push_back(Jet(FFTJetR20_pt[i],
					 FFTJetR20_eta[i],
					 FFTJetR20_phi[i],
					 FFTJetR20_m[i],				
					 FFTJetR20_splft[i],
					 FFTJetR20_mglft[i]
					 ));
		 fftJets30.push_back(Jet(FFTJetR30_pt[i],
					 FFTJetR30_eta[i],
					 FFTJetR30_phi[i],
					 FFTJetR30_m[i],				
					 FFTJetR30_splft[i],
					 FFTJetR30_mglft[i]
					 ));
		 fftJets40.push_back(Jet(FFTJetR40_pt[i],
					 FFTJetR40_eta[i],
					 FFTJetR40_phi[i],
					 FFTJetR40_m[i],				
					 FFTJetR40_splft[i],
					 FFTJetR40_mglft[i]
					 ));
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
	 }
	 //match jets to particles
	 vector<int> matchedJets20, matchedJets30,matchedJets40;
	 matchOneToOne(genparticles,fftJets20,DeltaRDistance(),&matchedJets20,1000.0);
	 matchOneToOne(genparticles,fftJets30,DeltaRDistance(),&matchedJets30,1000.0);
	 matchOneToOne(genparticles,fftJets40,DeltaRDistance(),&matchedJets40,1000.0);
	 /// fill ntuples variables 
   const unsigned numWDaughters(genparticles.size());
   double maxdr20(-1), maxdr30(-1), maxdr40(-1);
   double maxdrjj20(-1), maxdrjj30(-1), maxdrjj40(-1);
   TLorentzVector vec20,vec30,vec40,vecGen;
   getVecFromJets(genparticles,fftJets20,matchedJets20,distance,&vec20,&maxdr20,&maxdrjj20);
   getVecFromJets(genparticles,fftJets30,matchedJets30,distance,&vec30,&maxdr30,&maxdrjj30);
   getVecFromJets(genparticles,fftJets40,matchedJets40,distance,&vec40,&maxdr40,&maxdrjj40);
   for(unsigned i(0);i<numWDaughters;++i)
		{
      GenParticle gp(genparticles[i]);
      TLorentzVector tmpvec;
      tmpvec.SetPtEtaPhiM(gp.pt,gp.eta,gp.phi,gp.mass);
      vecGen += tmpvec;
		}
      ///
   ////Analyze fftjet collections
   Jet jet1, jet2;
   double mtmp(DBL_MAX), drjj(DBL_MAX);
   const bool verbose(false);
   int jetcheck[3]={-1,-1,-1};  //check the selection of jets, [0]=jet1 index,[1]=jet2 index,[2]=numjets(2,3,4)
   findBestJetsWithInvMass (ptCut,drCut,invariantWMass,fftJets20,&jet1,&jet2,&mtmp ,&drjj,verbose,jetcheck);
   findBestJetsWithInvMass (ptCut,drCut,invariantWMass,fftJets30,&jet1,&jet2,&mtmp ,&drjj,verbose,jetcheck);
   findBestJetsWithInvMass (ptCut,drCut,invariantWMass,fftJets40,&jet1,&jet2,&mtmp ,&drjj,verbose,jetcheck);
   if(jetcheck[0]>-1)
     analyzeJetSelection(jetcheck,analyzejetsel);
   //cout<<"jetcheck "<<jetcheck[0]<<":"
   //                 <<jetcheck[1]<<":"
   //                 <<jetcheck[2]<<":"
   //<<endl;

   TLorentzVector jfft(jet1.getLorentzVector()+jet2.getLorentzVector());
	 Jet V_fft(jfft.Pt(),jfft.Eta(),jfft.Phi(),jfft.M());
   
	 const double mass_V_fft(V_fft.mass),	
				 pt_V_fft  (V_fft.pt),
				 eta_V_fft (V_fft.eta),	
				 dr_VG_fft (V_fft.getLorentzVector().DeltaR(photons[0].getLorentzVector())),
				 dphi_VG_fft (abs(DeltaPhiDistance(V_fft.phi,photons[0].phi))),
				 mass_VG_fft ((V_fft.getLorentzVector()+
							 photons[0].getLorentzVector()).M() );
	 const double mass_V (grJets[0].mass),
				 pt_V   (grJets[0].pt),
				 eta_V  (grJets[0].eta),
				 pt_G   (photons[0].pt),
				 eta_G  (photons[0].eta),
				 dr_VG  (grJets[0].getLorentzVector().DeltaR(photons[0].getLorentzVector())),
				 dphi_VG  (abs(DeltaPhiDistance(grJets[0].phi,photons[0].phi))),
				 photonIso(photons[0].photonIso), 
				 tau2tau1 (GroomedJet_CA8_tau2tau1[0]);
	 const double mass_VG((grJets[0].getLorentzVector()+
				 photons[0].getLorentzVector()).M()
			 );

	 ntup_vg_vec.clear();
	 ntup_vg_vec.push_back(event_runNo);
	 ntup_vg_vec.push_back(event_evtNo);
	 ntup_vg_vec.push_back(event_lumi);
	 ntup_vg_vec.push_back(event_HT);
	 ntup_vg_vec.push_back(mass_V );
	 ntup_vg_vec.push_back(pt_V   );
	 ntup_vg_vec.push_back(GroomedJet_CA8_pt[0]);
	 ntup_vg_vec.push_back(eta_V  );
	 ntup_vg_vec.push_back(pt_G   );
	 ntup_vg_vec.push_back(eta_G  );
	 ntup_vg_vec.push_back(dr_VG  );
	 ntup_vg_vec.push_back(dphi_VG);
	 ntup_vg_vec.push_back(mass_VG);
	 ntup_vg_vec.push_back(mass_V_fft);
	 ntup_vg_vec.push_back(pt_V_fft);
	 ntup_vg_vec.push_back(eta_V_fft);
	 ntup_vg_vec.push_back(drjj);
	 ntup_vg_vec.push_back(dr_VG_fft);
	 ntup_vg_vec.push_back(dphi_VG_fft);
	 ntup_vg_vec.push_back(mass_VG_fft);
	 ntup_vg_vec.push_back(photonIso);
	 ntup_vg_vec.push_back(tau2tau1);
	 assert(ntup_vg_vec.size()==22);
	 ntup_VG->Fill(&ntup_vg_vec[0]);
	 //print out the event num
	 if(!(jentry%200000))
		 cout<<"Processing Event "<<(float)jentry<<endl;
	}
	cout<<"Done\n";
  f_out = new TFile(outfile.c_str(),"RECREATE");
	f_out->cd();
	ntup_VG->Write();
	f_out->Write();
	f_out->Close();

//analyze jet selection
	cout<<"\nTotal Entries "<<nentries
		<<"\nLA2 selected jet1 jet2  "<<analyzejetsel[0] 
		<<"\nLA3 selected jet1 jet2  "<<analyzejetsel[1] 
		<<"\nLA3 selected jet1 jet3  "<<analyzejetsel[2] 
		<<"\nLA3 selected jet2 jet3  "<<analyzejetsel[3] 
		<<"\nLA4 selected jet1 jet2  "<<analyzejetsel[4] 
		<<"\nLA4 selected jet1 jet3  "<<analyzejetsel[5] 
		<<"\nLA4 selected jet1 jet4  "<<analyzejetsel[6] 
		<<"\nLA4 selected jet2 jet3  "<<analyzejetsel[7] 
		<<"\nLA4 selected jet2 jet4  "<<analyzejetsel[8] 
		<<"\nLA4 selected jet3 jet4  "<<analyzejetsel[9] 
		<<endl;
	}
