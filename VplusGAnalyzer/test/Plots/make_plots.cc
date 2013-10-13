#include <TH2.h>
#include "fftplots_procedures.hh"
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>

void ScaleToInt(TH1D* h) {
h->Scale(1.0/h->Integral());
}
void DrawHist1dRes(TH1D* h_fft,string xtitle,float* etaLim)
{
	//find max
	float maxY1 = h_fft->GetMaximum();
	float maxX1 = h_fft->GetXaxis()->GetXmax();
  //
	h_fft->SetStats(kFALSE); 
	h_fft->SetTitle("");
	h_fft->GetXaxis()->SetTitle(xtitle.c_str());
	h_fft->SetLineColor(kBlue);
	h_fft->SetMaximum(maxY1*1.15);
	h_fft->Draw(); 
	//write stats 
	ostringstream statstr;
	TLatex stats;
	stats.SetTextAlign(21);
	stats.SetTextSize(0.03);
	stats.SetTextFont(42);
	statstr<<"M_{FFT} "<<getMedian(*h_fft);
	stats.DrawLatex(maxX1*0.75,maxY1*0.45,statstr.str().c_str());statstr.str("");
	statstr<<"#sigma_{FFT} "<<getRMS(*h_fft);
	stats.DrawLatex(maxX1*0.75,maxY1*0.37,statstr.str().c_str());statstr.str("");
	statstr<<etaLim[0]<<"<#left|#eta#right|<"<<etaLim[1];
	stats.DrawLatex(maxX1*0.5,maxY1*1.05,statstr.str().c_str());statstr.str("");
}
void DrawHist1dRes(TH1D* h_fft,TH1D* h_ak5,string xtitle,float* etaLim)
{
	//find max
	float maxY1 = max (h_fft->GetMaximum(),h_ak5->GetMaximum());
	float maxX1 = h_fft->GetXaxis()->GetXmax()+h_ak5->GetXaxis()->GetXmin();
  //
	h_fft->SetStats(kFALSE); 
	h_ak5->SetStats(kFALSE); 
	h_fft->SetTitle("");
	h_fft->GetXaxis()->SetTitle(xtitle.c_str());
	h_fft->SetLineColor(kBlue);
	h_fft->SetMaximum(maxY1*1.15);
	h_ak5->SetLineColor(kRed);
	h_fft->Draw(); h_ak5->Draw("same");
	//write stats 
	ostringstream statstr;
	TLatex stats;
	stats.SetTextAlign(21);
	stats.SetTextSize(0.03);
	stats.SetTextFont(42);
	statstr<<"M_{FFT} "<<getMedian(*h_fft);
	stats.DrawLatex(maxX1*0.75,maxY1*0.45,statstr.str().c_str());statstr.str("");
	statstr<<"#sigma_{FFT} "<<getRMS(*h_fft);
	stats.DrawLatex(maxX1*0.75,maxY1*0.37,statstr.str().c_str());statstr.str("");
	statstr<<"M_{ak5}  "<<getMedian(*h_ak5);
	stats.DrawLatex(maxX1*0.75,maxY1*0.30,statstr.str().c_str());statstr.str("");
	statstr<<"#sigma_{ak5}  "<<getRMS(*h_ak5);
	stats.DrawLatex(maxX1*0.75,maxY1*0.23,statstr.str().c_str());statstr.str("");
	statstr<<etaLim[0]<<"<#left|#eta#right|<"<<etaLim[1];
	stats.DrawLatex(maxX1*0.5,maxY1*1.05,statstr.str().c_str());statstr.str("");
	//draw legend
	TLegend* leg = new TLegend(0.65,0.70,0.90,0.89); // xy coordinates
	leg->AddEntry(h_fft ,"FFT","L");
	leg->AddEntry(h_ak5 ,"AK5","L");
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.03);
	leg->Draw("same");
}


void DrawHist1d(TH1D* sg,TH1D* bk,string xtitle)
{
  sg->SetStats(kFALSE); 
  sg->SetTitle("");
  sg->GetXaxis()->SetTitle(xtitle.c_str());
  sg->SetLineColor(kRed);
  bk->SetLineColor(kBlack);
  sg->Draw(); bk->Draw("same");
}
void DrawHist1d(TH1D* h,string xtitle)
{
  h->SetStats(kFALSE); 
  h->SetTitle("");
  h->GetXaxis()->SetTitle(xtitle.c_str());
  h->SetLineColor(kBlack);
  h->Draw();
}


void DrawHist2d(TH2D* sg,TH2D* bk,string xtitle, 
								string ytitle)
{
  sg->SetStats(kFALSE); bk->SetStats(kFALSE); 
  sg->SetTitle(""); 
  sg->GetXaxis()->SetTitle(xtitle.c_str());
  sg->GetYaxis()->SetTitle(ytitle.c_str());
  sg->SetFillColor(kRed);    sg->SetMarkerColor(kRed);
  bk->SetFillColor(kBlack);  bk->SetMarkerColor(kBlack);
  sg->SetMarkerSize(0.5); sg->SetMarkerStyle(20);
  bk->SetMarkerSize(0.5); bk->SetMarkerStyle(20);
  sg->Draw("box"); bk->Draw("box same");
}
void DrawHist2d(TH2D* sg,string xtitle, 
								string ytitle)
{
  sg->SetStats(kFALSE);
	sg->SetTitle(""); 
	sg->GetXaxis()->SetTitle(xtitle.c_str());
	sg->GetYaxis()->SetTitle(ytitle.c_str());
	sg->SetFillColor(kRed);    sg->SetMarkerColor(kRed);
	sg->SetMarkerSize(0.05); sg->SetMarkerStyle(20);
	sg->GetYaxis()->SetTitleOffset(1.25);
	sg->Draw("colz"); 
}

int make_test_plots() {

	TFile *f_sg       = new TFile("../test.root","read");
	TNtuple* ntup_VG =  (TNtuple*)f_sg->Get("ntup_VG");
	//declare hists
	TH1D* h_mV  = new TH1D("h_mV","h_mV",  50,0,150);  
	TH1D* h_ptV = new TH1D("h_ptV","h_ptV",100,0,300);  
	TH1D* h_ptG = new TH1D("h_ptG","h_ptG",100,0,300);  
	TH1D* h_dphiVG = new TH1D("h_dphiVG","h_dphiVG",100,-7,7);  
	//draw hists
	ntup_VG->Draw("m_V_gr>>h_mV"); 
	ntup_VG->Draw("pt_V_gr>>h_ptV"); 
	ntup_VG->Draw("pt_G>>h_ptG"); 
	ntup_VG->Draw("dphi_VG>>h_dphiVG"); 
	//scale hists
	ScaleToInt(h_mV); ScaleToInt(h_ptV);
	ScaleToInt(h_ptG); 
	//
	//Draw Hist
	TCanvas* c1 = new TCanvas("c1","c1",0,0,600,600);
	DrawHist1d(h_mV,"m^{jet} [Gev/c^{2}]");
	c1->SaveAs("mass_V.png");
	//
  TCanvas* c2 = new TCanvas("c2","c2",0,0,600,600);
	DrawHist1d(h_ptV,"Jet Pt[Gev/c]");
	c2->SaveAs("pt_V.png");
	//
  TCanvas* c3 = new TCanvas("c3","c3",0,0,600,600);
	DrawHist1d(h_ptG,"Photon Pt[Gev/c]");
	c3->SaveAs("pt_G.png");
	//
  TCanvas* c4 = new TCanvas("c4","c4",0,0,600,600);
	DrawHist1d(h_dphiVG,"dPhi Photon-Jet");
	c4->SaveAs("dphi_VG.png");
}


int make_radvsres_plots() {

	const double ptCut(20);
	const double drCut(0.2);
	//radius value s
	double radius[] = {20,30,40};
	//eta bins 
	double etaBins[] = {0.0,3.3};
	const unsigned numRad(sizeof(radius)/sizeof(radius[0]));
	const unsigned numEta ( sizeof(etaBins)/sizeof(etaBins[0])) ;
	//input file
	const string dir("../");
 	//declare and draw  hists
		TH1D* hW1_fft[numEta-1][numRad];
		TCanvas* cw[numEta-1][numRad];
		TCanvas* ct[numEta-1][numRad];
	TFile* f1 = new TFile((dir+"/test.root").c_str(),"read");
  //resolution arrays
	double wres[numEta-1][numRad], tres[numEta-1][numRad];
	double werr[numEta-1][numRad], terr[numEta-1][numRad],
				 rerr[numEta-1][numRad];
  const double c(0.94247451);
	for(unsigned nrad(0);nrad<numRad;++nrad) {
    ostringstream inputNtuple; 
		inputNtuple<<"ntup_RadScan";
		TNtuple* ntup1 = (TNtuple*)f1->Get(inputNtuple.str().c_str());
			for(unsigned i(0);i<numEta-1;++i)
		{
			ostringstream oss_f,oss_a;
			ostringstream cut_f,cut_a,draw_f,draw_a;
			oss_f<<"hWf_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]<<"_rad_"<<radius[nrad]; // histos 
			draw_f<<"Wmass"<<radius[nrad]<<"/WmassGen>>"<<oss_f.str(); //draw strings
			cut_f<<"matchdr"<<radius[nrad]<<"<"<<drCut
						<<"&&maxdrVG"<<radius[nrad]<<">1.0"
						<<"&&matchJet1>-1&&matchJet2>-1&&numWDaus==2"
						<<"&&photonIso&&WptGen>200"		
			;  // cut string 
			//histos 
			hW1_fft[i][nrad] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0,5);
			ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
      //get resolutions
      //set width_w1f [ expr {($p84-$p15)*0.5/$median} ]
			wres[i][nrad] = getRMS(*hW1_fft[i][nrad])/getMedian(*hW1_fft[i][nrad]);
      werr[i][nrad] = wres[i][nrad]*c/hW1_fft[i][nrad]->GetEntries();
      rerr[i][nrad] = 0;
		}
		//Draw histograms
		for(unsigned i(0);i<numEta-1;++i)
		{
			const float etaLim[2] ={etaBins[i],etaBins[i+1]};
			//W
			ostringstream cn;
			cn<<"cw"<<i<<"_rad_"<<radius[nrad];
			cw[i][nrad] = new TCanvas(cn.str().c_str(),cn.str().c_str(),600,600);
			DrawHist1dRes(hW1_fft[i][nrad],"m^{W}(jj)",etaLim);
			//save histos
			cn.str(""); 
			cn<<"Wres_eta_"<<etaBins[i]<<"_rad_"<<radius[nrad];
			cw[i][nrad]->SaveAs((cn.str()+".png").c_str());
		}
	}//radius loop
	TCanvas*cres = new TCanvas("cres","cres",10,10,600,600);
	TGraphErrors* gr_W = new TGraphErrors(numRad,radius,wres[0],rerr[0],werr[0]);
  //set options
  gr_W->SetLineColor(1); gr_W->SetMarkerColor(2); gr_W->SetMarkerStyle(21);
  //draw multigraph
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Resolution Graph");
	mg->Add(gr_W);
	mg->Draw("APZC");
  //save the resolution plot
	cres->SaveAs("RadiusvsResolution.png");
}

int make_radvsresvslajets_plots() {

	const double ptCut(20);
	const double drCut(0.2);
	//radius value s
	double radius[] = {20,30,40,50};
  double laVals[] = {2,3,4};
  const unsigned numLA(sizeof(laVals)/sizeof(laVals[0]));
	const unsigned numRad(sizeof(radius)/sizeof(radius[0]));
	//input file
	const string dir("../");
 	//declare and draw  hists
		TH1D* hW1_fft[numLA][numRad];
		TCanvas* cw[numLA][numRad];
		TCanvas* ct[numLA][numRad];
//	TFile* f1 = new TFile((dir+"/VplusGTreeAnalyzer_fftjetTest.root").c_str(),"read");
	TFile* f1 = new TFile((dir+"/VplusGTreeAnalyzer_fftjetTest_100Karis.root").c_str(),"read");
  //resolution arrays
	double wres[numLA][numRad];
	double werr[numLA][numRad],
				 rerr[numLA][numRad];
  const double c(0.94247451);
	for(unsigned nrad(0);nrad<numRad;++nrad) {
    ostringstream inputNtuple; 
		inputNtuple<<"ntup_RadScan";
		TNtuple* ntup1 = (TNtuple*)f1->Get(inputNtuple.str().c_str());
			for(unsigned i(0);i<numLA;++i)
		{
			ostringstream oss_f,oss_a;
			ostringstream cut_f,cut_a,draw_f,draw_a;
			oss_f<<"hWf_la_"<<laVals[i]<<"_rad_"<<radius[nrad]; // histos 
			draw_f<<"WmassLA"<<laVals[i]<<"R"<<radius[nrad]<<"/WmassGen>>"<<oss_f.str(); //draw strings
			cut_f<<"matchdrLA"<<laVals[i]<<"R"<<radius[nrad]<<"<"<<drCut
						<<"&&maxgendrVG>3.0"
						<<"&&allMatchedLA"<<laVals[i]<<"R"<<radius[nrad]
						<<"&&numWDaus==2"
						<<"&&photonIso&&WptGen>200"		
			;  // cut string 
			//histos 
			hW1_fft[i][nrad] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0,5);
			ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
      //get resolutions
			wres[i][nrad] = getRMS(*hW1_fft[i][nrad])/getMedian(*hW1_fft[i][nrad]);
      werr[i][nrad] = wres[i][nrad]*c/hW1_fft[i][nrad]->GetEntries();
      rerr[i][nrad] = 0;
      
      cout<<"Computed la "<<laVals[i]<<" rad "<<radius[nrad]<<endl;
      cout<<"Entries "<<hW1_fft[i][nrad]->GetEntries()<<endl;
      cout<<"W res "<< wres[i][nrad]<<endl;
		}
		//Draw histograms
		for(unsigned i(0);i<numLA;++i)
		{
			const float etaLim[2] ={99,99};
			//W
			ostringstream cn;
			cn<<"cw_la"<<laVals[i]<<"_rad_"<<radius[nrad];
			cw[i][nrad] = new TCanvas(cn.str().c_str(),cn.str().c_str(),600,600);
			DrawHist1dRes(hW1_fft[i][nrad],"m^{W}(jj)",etaLim);
			//save histos
			cn.str(""); 
			cn<<"Wres_la_"<<laVals[i]<<"_rad_"<<radius[nrad];
			cw[i][nrad]->SaveAs((cn.str()+"_la.png").c_str());
		}
	}//radius loop
  TCanvas*      cres[numLA];
  TGraphErrors* gr_W[numLA];
  TMultiGraph*  mg[numLA];
	//Draw TGraphs of resolutionsvsradius for each la n jets
	for(unsigned i(0);i<numLA;++i) 
	{
    ostringstream ocanvas; ocanvas<<"cres_la"<<laVals[i];
		cres[i] = new TCanvas(ocanvas.str().c_str(),ocanvas.str().c_str(),10,10,600,600);
		gr_W[i] = new TGraphErrors(numRad,radius,wres[i],rerr[i],werr[i]);
		//set options
		gr_W[i]->SetLineColor(1); gr_W[i]->SetMarkerColor(2); gr_W[i]->SetMarkerStyle(21);
    gr_W[i]->GetYaxis()->SetRangeUser(0.12,0.16);
		//draw multigraph
		mg[i] = new TMultiGraph();
    ostringstream lajets; lajets<<laVals[i];
		mg[i]->SetTitle(("Resolution Graph LA Jets="+lajets.str()).c_str());
		mg[i]->Add(gr_W[i]);
		mg[i]->Draw("APZC");
		//save the resolution plot
		cres[i]->SaveAs(("RadiusvsResolution_la"+lajets.str()+".png").c_str());
	}
	//Draw 2d color plot 
	TCanvas* ccol = new TCanvas("ccol","ccol",10,10,600,600);
  double bwla(0.5*(laVals[numLA-1]-laVals[0])/(numLA-1)),
				 bwrd(0.5*(radius[numRad-1]-radius[0])/(numRad-1));
	TH2F* hcol=new TH2F("hcol","hcol",numLA,laVals[0]-bwla,laVals[numLA-1]+bwla,
																		numRad,radius[0]-bwrd,radius[numRad-1]+bwrd);
  	for(unsigned nrad(0);nrad<numRad;++nrad) 
		for(unsigned i(0);i<numLA;++i)
			hcol->Fill(laVals[i],radius[nrad],wres[i][nrad]);
  hcol->SetStats(0);
  hcol->SetTitle("");
  hcol->SetXTitle("Requested Jets"); hcol->SetYTitle("Cone Radius"); 
	hcol->Draw("colz");
  ccol->SaveAs("RadiusvsResolutionCombined.png");
}		
