#define VerifyLHEEventClass_cxx
#include "VerifyLHEEventClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TLorentzVector.h"
#include <cassert>

void VerifyLHEEventClass::Loop()
{
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	vector<float> ntup_vg_vec;

	for (Long64_t jentry=0; jentry<nentries;++jentry) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;
    TLorentzVector wvec, gvec;
    wvec.SetPtEtaPhiM(Wpt,Weta,Wphi,Wmass);	
    gvec.SetPtEtaPhiM(photonpt,photoneta,photonphi,0);
  
		ntup_vg_vec.clear();
		ntup_vg_vec.push_back(Wmass );
		ntup_vg_vec.push_back(Wpt   );
		ntup_vg_vec.push_back(Wpt   );
		ntup_vg_vec.push_back(Weta  );
		ntup_vg_vec.push_back(photonpt);
		ntup_vg_vec.push_back(photoneta);
		ntup_vg_vec.push_back(photonWdr);
		ntup_vg_vec.push_back(photonWdphi);
		ntup_vg_vec.push_back((wvec+gvec).M());
		assert(ntup_vg_vec.size()==9);
		ntup_VG->Fill(&ntup_vg_vec[0]);
		//print out the event num
		if(!(jentry%200000))
			cout<<"Processing Event "<<jentry<<endl;
	}
	cout<<"Done\n";
	f_out = new TFile(outfile.c_str(),"RECREATE");
	f_out->cd();
	ntup_VG->Write();
	f_out->Write();
	f_out->Close();

}
