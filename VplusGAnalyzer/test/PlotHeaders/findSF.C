#include "samples.h"
#include <iostream> 
#include <iomanip> 
#include <TCanvas.h>
#include <TRegexp.h>
#include <TLegend.h>
#include <THStack.h>
#include <TLatex.h>
#include "setTDRStyle.C"
#include <TROOT.h>
#include "customize.h"
#include <cmath>
#include <sstream>
#include <fstream>
using namespace std;
//namespace std { using namespace std::tr1; }

void findSF()
{
	setTDRStyle();
	gROOT->ForceStyle();
	initOptions();
	// gStyle->SetOptFit(1111);
	// gStyle->SetOptStat(0);

	std::vector<Sample> s = samples();
	const unsigned numSamples(s.size());
	Sample data(1,1,"fake data","../VJetAnalyzer_data_contreg.root",0,true,10000);
	const double mHigh(105.),mLow(65.);

	for(size_t i=0;i< s.size();++i) if(s[i].data) {data=s[i];break;}
	data.file()->ls(); 
	for(size_t i=0;i< s.size();++i) {s[i].dump(data.lumi()); }

	std::vector<Variable> variables;
	variables.push_back(Variable("m_W_la","m_{W}",25,10,200));
	TTree* ntup[numSamples];
	for(unsigned j=0;j<numSamples ;++j)
		ntup[j]=(TTree*)(s[j].file()->Get("ntup_W"));
	const unsigned numVars(variables.size());
	//get mva cut values from the h600 samples
	Sample h600(1,1,"SignalMC","../VJetAnalyzer_h600_contreg.root",0,false);
	TTree* ntup_h600 = (TTree*)(h600.file()->Get("quantiles_ntup"));
	float qvalue(0),mvacut(0);  
	ntup_h600->SetBranchAddress("qvalue",&qvalue);
	ntup_h600->SetBranchAddress("cutvalue",&mvacut);
	//
	TH1F* hd; TH1F* h;
	TCanvas *c = new TCanvas();
	//loop over different quantile values 
	const unsigned numQuantiles (ntup_h600->GetEntries()); 
	for(unsigned i(0);i<numQuantiles-2;++i)
	{ 
		ntup_h600->GetEntry(i);
		TString n=variables[0].name;
		double events_precut_mc(0), events_postcut_mc(0);
		double events_precut_data(0), events_postcut_data(0);
		Options o;
    c->Clear();
		//cuts on ntuples
    ostringstream cutstr1,cutstr2; 
    cutstr1<<"m_W_la<"<<mHigh<<"&&m_W_la>"<<mLow;//masswindow
    cutstr2<<cutstr1.str()<<"&&mva3Value_fft>"<<mvacut;//mvacut
		//histo title
		ostringstream oss;
		oss<<variables[0].name<<"data"<<"q"<<qvalue*100;
     
		hd = new TH1F(oss.str().c_str(),oss.str().c_str(),
				variables[0].nbins,
				variables[0].xmin,
				variables[0].xmax);
		ostringstream oss2;
		oss2<<variables[0].name<<">>"<<oss.str();
		for(unsigned sam(0);sam<numSamples;++sam) 	
			if(s[sam].data)
		{
				ntup[sam]->Draw(oss2.str().c_str(),cutstr1.str().c_str());
				events_precut_data+=hd->Integral();
				ntup[sam]->Draw(oss2.str().c_str(),cutstr2.str().c_str());
				events_postcut_data+=hd->Integral();
		}
		hd->SetMarkerStyle(20);
		hd->Draw("E1");
		ostringstream ylabel; 
		ylabel<<"Events"<<" / "<<(float)hd->GetBinWidth(1)<<" GeV";

		hd->SetYTitle(ylabel.str().c_str()); 
		hd->SetXTitle(variables[0].xlabel.c_str()); 
		THStack * sta = new THStack("sta",hd->GetTitle());
		TLegend * l = new TLegend(o.legendx1,o.legendy1,o.legendx2,o.legendy2); //0.7,0.1,0.9,0.6);

		l->AddEntry(hd, "Data","LP");
		l->SetFillColor(kWhite);
		for(size_t j=0;j<numSamples ;++j) 
			if(!s[j].data) 
			{
				ostringstream oss3;
				oss3<<variables[0].name<<s[j].name<<"q"<<qvalue*100;
				h = new TH1F(oss3.str().c_str(),oss3.str().c_str(),
						variables[0].nbins,
						variables[0].xmin,
						variables[0].xmax);

				ostringstream oss4;
				oss4<<variables[0].name<<">>"<<oss3.str();
			  //mass window cut 
				ntup[j]->Draw(oss4.str().c_str(),cutstr1.str().c_str());
				h->Scale(s[j].scale(data.lumi()));
				events_precut_mc+=h->Integral();
				//mva cut
				ntup[j]->Draw(oss4.str().c_str(),cutstr2.str().c_str());
				h->Scale(s[j].scale(data.lumi()));
				events_postcut_mc+=h->Integral();

				h ->SetLineColor(s[j].color);
				h->SetFillColor(s[j].color);
				l->AddEntry(h,s[j].name.c_str(),"F");

				sta->Add(h);
			}

		hd->Draw("E1");
		sta->Draw("same");
		hd->Draw("E1same");
		cout<<"Quantile "<<qvalue<<" mvaCut "<<mvacut
		<<"\n SF Data "<<events_postcut_data/events_precut_data
		<<"\t SF MC   "<<events_postcut_mc/events_precut_mc
		<<"\t Data/MC SF  "<<(events_postcut_data/events_precut_data)/(events_postcut_mc/events_precut_mc)
		<<endl;
    cout<<" Data pre/post cut Events "<<events_postcut_data<<'/'<<events_precut_data
				<<"\n MC pre/post cut          "<<events_postcut_mc<<'/'<<events_precut_mc
				<<endl;

		float maxY1 = max (sta->GetMaximum(),hd->GetBinContent(hd->GetMaximumBin()));
		double maxX1 = hd->GetXaxis()->GetXmax()+hd->GetXaxis()->GetXmin();
		hd->GetYaxis()->SetRangeUser(0,maxY1*1.15);
		hd->SetTitle("");
		l->Draw();
	}// quantiles loop
}//findSF
