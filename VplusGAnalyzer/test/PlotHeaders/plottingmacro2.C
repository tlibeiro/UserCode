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

void plottingmacro2()
{
  setTDRStyle();
  gROOT->ForceStyle();
  initOptions();

std::vector<Sample> s = samples();
 const unsigned numSamples(s.size());
  Sample data(1,1,"fake data","/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput/VplusGTreeAnalyzer_GJets200.root",0,true,10000);

	for(size_t i=0;i<numSamples;++i) if(s[i].data) {data=s[i];break;}
	data.file()->ls(); 
	for(size_t i=0;i<numSamples;++i) s[i].dump(data.lumi());

//legend variables
  legendVars legvar("legvar",0.57,0.65,0.82,0.80);
  legendVars tauleg("tauleg",0.27,0.65,0.52,0.80);
  //for test 
	tauleg=legvar;

  const string tag("");
	std::vector<Variable> variables;
////for test 
//	variables.push_back(Variable("Wpt","W Pt [GeV/c]",50,0,1000));
//	variables.push_back(Variable("photonpt","#gamma Pt [GeV/c]",50,0,1000));
//
	variables.push_back(Variable("m_V_gr","pruned jet mass [GeV/c^{2}]",30,0,200));
//	variables.push_back(Variable("pt_V_gr","pruned jet Pt [GeV/c]",30,0,800));
//	variables.push_back(Variable("pt_V_ungr","jet Pt [GeV/c]",30,0,800));
//	variables.push_back(Variable("eta_V_gr","pruned jet #eta",20,-5.5,5.5));
//	variables.push_back(Variable("pt_G","#gamma Pt [GeV/c]",20,0,800));
//	variables.push_back(Variable("dr_VG","dR^{#gamma-Jet}",20,0,6));
//	variables.push_back(Variable("dphi_VG","dPhi^{#gamma-Jet}",20,0,3.2));
//	variables.push_back(Variable("m_VG","mass^{#gamma-Jet}",30,0,1500));
//	variables.push_back(Variable("tau2tau1","#tau_{2}/#tau_{1}",20,0,1.1));
//	variables.push_back(Variable("m_V_fft","fftjet mass [GeV/c^{2}]",50,0,200));
//	variables.push_back(Variable("pt_V_fft","fftjet Pt [GeV/c]",100,0,500));
//	variables.push_back(Variable("eta_V_fft","fftjet #eta",20,-5.5,5.5));
//	variables.push_back(Variable("dr_VG_fft","dR^{#gamma-fftJet}",20,0,5));
//	variables.push_back(Variable("dphi_VG_fft","dPhi^{#gamma-fftJet}",20,0,3.2));
//	variables.push_back(Variable("m_VG_fft","mass^{#gamma-fftJet}",100,0,1000));

	const unsigned numVars(variables.size());
	TCanvas *c;
	TH1F* histo_arr[numVars][numSamples];
	TTree* ntup[numSamples];
	double maxY[numVars];
	TLegend* leg[numVars];
	for(unsigned j=0;j<numSamples ;++j)
		ntup[j]=(TTree*)(s[j].file()->Get("ntup_VG"));
 //for test
//		ntup[j]=(TTree*)(s[j].file()->Get("verifylheevents/verifylhe"));
	for(size_t i = 0 ; i<numVars ; ++i) 
	{
		TString n = variables[i].name;
		string n2 = variables[i].name;
		Options o;
		TLegend * l;
		if(variables[i].name.find("tau")!=string::npos) 
				 l=new TLegend(tauleg.xmin,tauleg.ymin,
											 tauleg.xmax,tauleg.ymax);
		else l=new TLegend(0.6,0.7,0.85,0.85);

		leg[i] = l;
		leg[i]->SetFillColor(0);
		leg[i]->SetTextFont(42);
		leg[i]->SetTextSize(0.03);
		maxY[i]=0;    

		for(size_t j=0;j<numSamples ;++j) 
		{ 
			cout<<"Making Sample "<<s[j].name<<" Variable "<<variables[i].name<<endl;
			ostringstream oss;
			oss<<variables[i].name<<s[j].name;
			TH1F* h = new TH1F(oss.str().c_str(),oss.str().c_str(),
					variables[i].nbins,
					variables[i].xmin,
					variables[i].xmax);

			ostringstream oss2;
			oss2<<variables[i].name<<">>"<<oss.str();
			ostringstream cutstr; cutstr<<"1";
			//quality cuts, not for HT plot
			if(variables[i].name.find("HT")==string::npos)
			{  
				cutstr
					<<"&&photonIso==1"
					<<"&&pt_G>250"
//					<<"&&tau2tau1<0.5"
					;
				//dr, pt  cuts 
				if(variables[i].name.find("fft")==string::npos)
					cutstr
						<<"&&dr_VG>3"
						<<"&&dphi_VG>2.9"
//						<<"&&pt_V_ungr>260"
//						<<"&&m_VG>600"
						;
				else if (variables[i].name.find("fft")!=string::npos)
					cutstr
						<<"&&dr_VG_fft>1"
						<<"&&pt_V_fft>200"
						;			
			}

			ntup[j]->Draw(oss2.str().c_str(),cutstr.str().c_str());
			h->Scale(s[j].scale(data.lumi()));
			//for shape ccomparision only 
//			if(n2.find("HT")==string::npos)
//				h->Scale(1./h->Integral());
			h->SetLineColor(s[j].color);
			h->SetLineWidth(2);
			ostringstream ylabel;
			ylabel<<o.yaxis<<" / "<<h->GetBinWidth(1)<<" ";
			h->GetYaxis()->SetTitle(ylabel.str().c_str());
			h->GetXaxis()->SetTitle(variables[i].xlabel.c_str());
			h->SetTitle("");
			//check for maximum Y
			double maxy=h->GetBinContent(h->GetMaximumBin());
			if(maxy>maxY[i]) 
				maxY[i] = maxy;       
			histo_arr[i][j] = h ;

			if(!s[j].data) 
				leg[i]->AddEntry(h,(s[j].name).c_str(),"L");
		} //samples loop
	}//varibales names loop

	for(size_t i = 0 ; i<numVars ; ++i) 
	{
		TCanvas* c = new TCanvas();
		for(size_t j=0;j<numSamples;++j)  
		{
			if(variables[i].name.find("HT")!=string::npos)
				histo_arr[i][j]->GetYaxis()->SetRangeUser(0,maxY[i]*100);
			else 
				histo_arr[i][j]->GetYaxis()->SetRangeUser(0,maxY[i]*1.2);
			histo_arr[i][j]->SetMinimum(1e-3);
			if(!j && (!s[j].data)) 
				histo_arr[i][j]->Draw();
			else if(!s[j].data)
				histo_arr[i][j]->Draw("same");
			if(variables[i].name.find("HT")!=string::npos
				||	variables[i].name.find("pt_G")!=string::npos
				||	variables[i].name.find("m_V")!=string::npos
				)
				c->SetLogy(1);
			leg[i]->Draw();
		}
		string outfile    = "Plots/Dists/"+variables[i].name+tag+"_dists.png";
		string outfile2   = "Plots/Dists/"+variables[i].name+tag+"_dists.eps";
		c->SaveAs(outfile.c_str());
		c->SaveAs(outfile2.c_str());
	}//


}//plottingmacro
