#ifndef VPLUSGANALYZER_HH
#define VPLUSGANALYZER_HH
#include <algorithm>
#include <vector>
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
using namespace std;

template<class T>
class TreeFiller {
	public:
		static const unsigned maxJetsOrPhotons=4;
		TreeFiller(const std::string& name,
				TTree* t);
		void fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag);
	private:
		TreeFiller(){ };
		/// Helper function for main constructor
		void SetBranch( float* x, const char* name,unsigned maxJetsOrPhotons);
		void SetBranch( int* x, const char* name, unsigned maxJetsOrPhotons);
		void SetBranchSingle( float* x, const std::string& name);
		void SetBranchSingle( int* x, const std::string& name);
    void SetBranches();
		TTree* tree;
		edm::Handle<T> jetsOrPhotons;
		float jetpt  [maxJetsOrPhotons];
		float jetphi [maxJetsOrPhotons];
		float jeteta [maxJetsOrPhotons];
		float jetm   [maxJetsOrPhotons];
    const std::string name;
    bool isGen;
};
	template <class T>
TreeFiller<T>::TreeFiller(const std::string& n,
		TTree* t)
:name(n)
{
	assert(t);
	tree = t;
  isGen = (name.find("Gen")!=std::string::npos);
	for(unsigned i(0);i<maxJetsOrPhotons;++i)
	{
		jetpt[i] =-1;
		jetphi[i]=-1;
		jeteta[i]=-1;
		jetm[i]  =-1;
	}
  SetBranches();
}
//
template<class T>
void TreeFiller<T>::SetBranches()
{
  //set branches 
	SetBranch(tree,jetpt  ,(name+"_pt").c_str(),maxJetsOrPhotons);
	SetBranch(tree,jeteta ,(name+"_eta").c_str(),maxJetsOrPhotons);
	SetBranch(tree,jetphi ,(name+"_phi").c_str(),maxJetsOrPhotons);
	SetBranch(tree,jetm   ,(name+"_m").c_str(),maxJetsOrPhotons);
}
	template<class T>
void TreeFiller<T>::fillBranch(const edm::Event& iEvent,const edm::InputTag& inputTag) 
{
	iEvent.getByLabel(inputTag,jetsOrPhotons); 
	const unsigned numJetsOrPhotons(jetsOrPhotons->size());
	for(unsigned i(0);i<maxJetsOrPhotons;++i)
		if(i<numJetsOrPhotons)
		{
			jetpt[i] =jetsOrPhotons->at(i).pt();
			jeteta[i]=jetsOrPhotons->at(i).eta();
			jetphi[i]=jetsOrPhotons->at(i).phi();
			jetm[i]  =jetsOrPhotons->at(i).mass();
		}
		else 
		{
			jetpt[i] = -1.;
			jeteta[i]= -1.;
			jetphi[i]= -1.;
			jetm[i]  = -1.;
		}
}

/////// Helper for above function ////////////////////////////////
//////////////////////////////////////////////////////////////////
template<class T>
void TreeFiller<T>::SetBranchSingle(float* x, const std::string& name)
{
	tree->Branch( name.c_str(), x, ( name+"/F").c_str() );
}
//
template<class T>
void TreeFiller<T>::SetBranchSingle(int* x, const std::string& name)
{
	tree->Branch( name.c_str(), x, ( name+"/I").c_str() );
}
//
template<class T>
void TreeFiller<T>::SetBranch(float* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/F";
	tree->Branch( name, x, oss.str().c_str() );
}
//
template<class T>
void TreeFiller<T>::SetBranch(int* x, const char* name,unsigned maxJetsOrPhotons)
{
	std::ostringstream oss;
	oss<<name<<"["<<maxJetsOrPhotons<<"]"<<"/I";
	tree->Branch( name, x, oss.str().c_str() );
}
#endif
