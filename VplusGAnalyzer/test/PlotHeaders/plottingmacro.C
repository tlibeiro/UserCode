#include "samples.h"
#include <iostream> 
#include <iomanip> 
#include <TCanvas.h>
#include <TRegexp.h>
#include <TLegend.h>
#include <THStack.h>
#include <TF1.h>
#include <TLatex.h>
#include "setTDRStyle.C"
#include <TROOT.h>
#include "customize.h"
#include <cmath>
#include <sstream>
#include <fstream>
using namespace std;

void plottingmacro()
{
	setTDRStyle();
	gROOT->ForceStyle();
	// gStyle->SetOptFit(1111);
	// gStyle->SetOptStat(0);

	std::vector<Sample> s = samples();
  const unsigned numSamples(s.size());
	Sample data(1,1,"fake data","../VplusGTreeAnalyzer_GJets200.root",0,true,10000);
  const string output_prefix("");
	double totalEv_data(0), totalEv_mc(0);

	for(size_t i=0;i< s.size();++i) if(s[i].data) {data=s[i];break;}
	data.file()->ls(); 
	for(size_t i=0;i< s.size();++i) {s[i].dump(data.lumi()); }

//legend variables
  legendVars legvar("legvar",0.57,0.65,0.82,0.80);
  legendVars tauleg("tauleg",0.27,0.65,0.52,0.80);
  //for test 
	tauleg=legvar;

	std::vector<Variable> variables;
//	variables.push_back(Variable("HT","HT",100,100,1200));
	variables.push_back(Variable("m_V_gr","pruned jet mass [GeV/c^{2}]",30,0,200));
	variables.push_back(Variable("pt_V_gr","pruned jet Pt [GeV/c]",30,0,800));
	variables.push_back(Variable("pt_V_ungr","jet Pt [GeV/c]",30,0,800));
	variables.push_back(Variable("eta_V_gr","pruned jet #eta",20,-5.5,5.5));
	variables.push_back(Variable("pt_G","#gamma Pt [GeV/c]",20,0,800));
	variables.push_back(Variable("dr_VG","dR^{#gamma-Jet}",20,0,6));
	variables.push_back(Variable("dphi_VG","dPhi^{#gamma-Jet}",20,0,3.2));
	variables.push_back(Variable("m_VG","mass^{#gamma-Jet}",30,0,1500));
	variables.push_back(Variable("tau2tau1","#tau_{2}/#tau_{1}",20,0,1.1));

	TTree* ntup[numSamples];
	for(unsigned j=0;j<numSamples ;++j)
		ntup[j]=(TTree*)(s[j].file()->Get("ntup_VG"));
	const unsigned numVars(variables.size());

	//	TCanvas *c = new TCanvas();
	for(size_t i = 0 ; i<numVars ; ++i) 
	{
		TString n=variables[i].name;
		//		c->Clear();
		TCanvas *c = new TCanvas();
		c->SetTitle(variables[i].name.c_str());
		ostringstream oss;
		oss<<variables[i].name<<"data";
		TH1F* hd = new TH1F(oss.str().c_str(),oss.str().c_str(),
				variables[i].nbins,
				variables[i].xmin,
				variables[i].xmax);
		ostringstream oss2;
		oss2<<variables[i].name<<">>"<<oss.str();
		ostringstream cutstr; cutstr<<"1";
      //quality cuts, not for HT plot
		if(variables[i].name.find("HT")==string::npos)
      cutstr
   			<<"&&photonIso==1"
        <<"&&pt_G>250"
				<<"&&tau2tau1<0.5"
			;
			//dr, pt  cuts 
			if(variables[i].name.find("fft")==string::npos)
				cutstr
				<<"&&dr_VG>3"
        <<"&&dphi_VG>2.9"
				<<"&&pt_V_ungr>260"
				<<"&&m_VG>600"
				;
     else if (variables[i].name.find("fft")!=string::npos)
				cutstr
				<<"&&dr_VG_fft>1"
				<<"&&pt_V_fft>200"
				;

		for(unsigned sam(0);sam<numSamples;++sam) 	
			if(s[sam].data)
				ntup[sam]->Draw(oss2.str().c_str(),cutstr.str().c_str());
		hd->SetMinimum(1e-3);
		hd->SetMarkerStyle(21);
		hd->Draw("E1");
		ostringstream ylabel; 
		ylabel<<"Events"<<" / "<<(float)hd->GetBinWidth(1);
		if(variables[i].name.find("pt")!=string::npos ||
				variables[i].name.find("m_")!=string::npos)
			ylabel<<" GeV";

		hd->SetYTitle(ylabel.str().c_str()); 
		hd->SetXTitle(variables[i].xlabel.c_str()); 
		THStack * sta = new THStack("sta",hd->GetTitle());
		TLegend * l;
		if(variables[i].name.find("tau")!=string::npos) 
				 l=new TLegend(tauleg.xmin,tauleg.ymin,
											 tauleg.xmax,tauleg.ymax);
		else l=new TLegend(0.6,0.7,0.85,0.85);
		l->SetTextFont(42);
		l->SetTextSize(0.035);
		//declare th1f for stack 
		ostringstream osshsta;
		osshsta<<"hstack"<<variables[i].name;
		TH1F*	hstack = new TH1F(osshsta.str().c_str(),
				osshsta.str().c_str(),
				variables[i].nbins,
				variables[i].xmin,
				variables[i].xmax);  
	  l->AddEntry(hd, "Data","LPE1");
		l->SetFillColor(kWhite);
		for(size_t j=0;j<numSamples ;++j) 
			if(!s[j].data) 
			{
				ostringstream oss3;
				oss3<<variables[i].name<<s[j].name;
				TH1F* h = new TH1F(oss3.str().c_str(),oss3.str().c_str(),
						variables[i].nbins,
						variables[i].xmin,
						variables[i].xmax);

				ostringstream oss4;
				oss4<<variables[i].name<<">>"<<oss3.str();
				ntup[j]->Draw(oss4.str().c_str(),cutstr.str().c_str());
		    h->SetMinimum(1e-3);
				h->Scale(s[j].scale(data.lumi()));
				h ->SetLineColor(s[j].color);
				h->SetFillColor(s[j].color);
				h->SetTitle("");
				l->AddEntry(h,s[j].name.c_str(),"F");

				sta->Add(h);
				hstack->Add(h);
				//	      h->Draw("same");
				//for test
				//			totalEv_mc+=h->Integral();
			}

		TCanvas *c2 = new TCanvas();
		hd->Draw("P");
		sta->Draw("same");
		hd->Draw("Psame");
		//for test
		//		totalEv_data+=hd->Integral();
		//		cout<<variables[i].name<<"\nTotal Data/MC "<<totalEv_data
		//				<<"/"<<totalEv_mc
		//				<<endl;
		//		totalEv_data=0; totalEv_mc=0;

		float maxY1  = max (sta->GetMaximum(),hd->GetBinContent(hd->GetMaximumBin()));
		float maxX1 = hd->GetBinLowEdge(hd->GetNbinsX());
		float minX1 = hd->GetBinLowEdge(1);
		hd->GetYaxis()->SetRangeUser(0,maxY1*1.15);
		hd->SetTitle("");

		ostringstream lumioss,ksoss;
		lumioss<<"#sqrt{s}=8TeV, L="<<data.lumi()/1000<<"/fb^{-1}";
		ksoss<<"KS="<<hd->KolmogorovTest(hstack);
		TLatex latex;
		latex.SetTextSize(0.04);
		latex.SetTextFont(42);
		latex.DrawLatex(maxX1*0.1+minX1,maxY1*1.05,lumioss.str().c_str());
		l->Draw();

		string outfile  = "Plots/"+variables[i].name+output_prefix+".png";
		string outfile2 = "Plots/"+variables[i].name+output_prefix+".eps";
		c2->SaveAs(outfile.c_str());
		c2->SaveAs(outfile2.c_str());
	}//variables names loop
}//plottingmacro
