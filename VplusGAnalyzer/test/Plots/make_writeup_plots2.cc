#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include "fftplots_procedures.hh"

void ScaleToInt(TH1D* h) {
	h->Scale(1.0/h->Integral());
}

void DrawHist1d(TH1D* h_fft,TH1D* h_ak5,string xtitle)
{
//find max
	float maxY1 = max (h_fft->GetMaximum(),h_ak5->GetMaximum());
	float maxX1 = h_fft->GetXaxis()->GetXmax()+h_ak5->GetXaxis()->GetXmin();

	h_fft->SetStats(kFALSE); 
	h_ak5->SetStats(kFALSE); 
	h_fft->SetTitle("");
	h_fft->GetXaxis()->SetTitle(xtitle.c_str());
	h_fft->SetLineColor(kBlue);
	h_fft->SetMaximum(maxY1*1.15);
	h_ak5->SetLineColor(kRed);
	h_fft->Draw(); h_ak5->Draw("same");
}
void DrawHist1dPU(vector<TH1F*>& hists,string& xtitle,string& ytitle)
{
	const unsigned numhists(hists.size());
	float maxY1(0), maxX1(0);

	for(unsigned i(0);i<numhists;++i) {
		maxY1 = max (hists[i]->GetMaximum(),maxY1);
		hists[i]->SetStats(kFALSE); 
		hists[i]->Rebin(2); 
		hists[i]->SetTitle("");
		hists[i]->GetXaxis()->SetTitle(xtitle.c_str());
		hists[i]->GetYaxis()->SetTitle(ytitle.c_str());
		hists[i]->SetLineColor(kAzure-9+3*i);
		if(!i) hists[i]->Draw("");
		else   hists[i]->Draw("same");
	}//hists loop
	hists[0]->SetMaximum(maxY1*2);
	//draw legend
	TLegend* leg = new TLegend(0.65,0.62,0.89,0.89); // xy coordinates
	for(unsigned i(0);i<numhists;++i) {
		ostringstream legstr;
		legstr<<"nPV "<<i;
		leg->AddEntry(hists[i],legstr.str().c_str(),"L");
	}
	leg->SetFillColor(0);  leg->SetTextFont(42);
	leg->SetBorderSize(0); leg->SetTextSize(0.03);
	leg->Draw("same");
}

void DrawHist1dPU(TH1F* hist,string& xtitle,string& ytitle)
{
	hist->SetStats(kFALSE); 
	hist->Rebin(2); 
	hist->SetTitle("");
	hist->GetXaxis()->SetTitle(xtitle.c_str());
	hist->GetYaxis()->SetTitle(ytitle.c_str());
	hist->SetLineColor(kBlack);
	hist->GetYaxis()->SetTitleOffset(1.2);
	hist->Draw("");
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
	sg->SetFillColor(kRed);    sg->SetMarkerColor(kGray+3);
	sg->SetMarkerSize(0.05);   sg->SetMarkerStyle(20);
	sg->GetYaxis()->SetTitleOffset(1.5);
	sg->Draw(""); 
}
void getShadedAreaFromTF1(const TF1*flow, const TF1*fhigh,TGraph* gr)
{
	const int num = gr->GetN()/2;
	double xmax,xmin;
	flow->GetRange(xmin,xmax);
	const double binWidth = (xmax-xmin)/num;
	for(unsigned i(0);i<num;++i)
	{	
		gr->SetPoint(i,binWidth*i,fhigh->Eval(binWidth*i));
		gr->SetPoint(num+i,xmax-binWidth*i,flow->Eval(xmax-binWidth*i));
	}
  TColor* col = gROOT->GetColor(kRed);
  TColor* mycolor = new TColor(1700,col->GetRed(),col->GetGreen(),
																	 col->GetBlue(),"",0.5);
  gr->SetFillStyle(1001);
  gr->SetFillColor(1700);
	gr->Draw("f");
}

int make_pu_plots() {
	//input files
	const string dir("/uscms_data/d3/tlibeiro/PileupCHS/CMSSW_5_3_8/src/Pile_study_Igor/PileupCorrections/PFchs/MC");
	TFile* f1 = new TFile((dir+"/pileup_chs_calibration_pf_PUS7.root").c_str(),"read");
	TFile* f2 = new TFile((dir+"/pfclean_2_mc_fit_0.4_True_bw_0.03.root").c_str(),"read");
	//get ntuple, hists 
	const unsigned nPV(46);
	TH1F* meanPU[nPV];
	for(unsigned i(0);i<nPV;++i)
	{
		ostringstream oss;
		oss<<"Eta_Dependence_of_Mean_Pileup__nPV___"<<i;
		meanPU[i]=(TH1F*)f1->Get(oss.str().c_str());
	}
	TH1F* etDensity = (TH1F*)f1->Get("Et_Scale_Factors");
	TH1F* etFactors = (TH1F*)f1->Get("Et_Flattening_Factors");
	TTree* predResp = (TTree*)f2->Get("Predictor-response_ntuple");
	TTree* regress  = (TTree*)f2->Get("Regression");
	TH2F*  estimate1= new TH2F("estimate1","estimate1",100,0,30,100,0,0.06);
	TH2F*  estimate2= new TH2F("estimate2","estimate2",100,0,30,100,0,0.06);
	TH2F*  estimate3= new TH2F("estimate3","estimate3",100,0,30,100,0,0.06);
	TH2F*  estimate4= new TH2F("estimate4","estimate4",100,0,30,100,0,0.06);
	//draw hists
	const unsigned drawPU(6);
	vector<TH1F*> meanPUDraw;
	for(unsigned i(0);i<drawPU;++i)
		meanPUDraw.push_back(meanPU[i]);
  //eta vs mean ET, nPV
	TCanvas* c1 = new TCanvas("c1","c1",600,600);
	DrawHist1dPU(meanPUDraw,"#eta","<E_{T}>");
  c1->SaveAs("eta_vs_et.eps");
	// eta vs rho 
	TCanvas* c2 = new TCanvas("c2","c2",600,600);
	DrawHist1dPU(etDensity,"#eta","E_{T} Density");
  c2->SaveAs("eta_vs_rho.eps");
	// eta vs scale factors
	TCanvas* c3 = new TCanvas("c3","c3",600,600);
	DrawHist1dPU(etFactors,"#eta","Scale Factor");
  c3->SaveAs("eta_vs_scf.eps");
	///regression plot
	TCanvas* c4 = new TCanvas("c4","c4",600,600);
	predResp->Draw("estimate:gridEtDensity>>estimate1","1","P");
	regress->Draw("q1_0:gridEtDensity>>estimate2","1","");
	regress->Draw("q1_1:gridEtDensity>>estimate3","1","");
	regress->Draw("q1_2:gridEtDensity>>estimate4","1","");
	TF1* fx1 = new TF1("fx1","[0]*x+[1]",0,30);  
	TF1* fx2 = new TF1("fx2","[0]*x+[1]",0,30);  
	TF1* fx3 = new TF1("fx3","[0]*x+[1]",0,30);  
	estimate2->Fit("fx1");
	estimate3->Fit("fx2");
	estimate4->Fit("fx3");
	TGraph* gr = new TGraph(estimate2->GetXaxis()->GetNbins()*2);
	TCanvas* c5 = new TCanvas("c5","c5",600,600);
	DrawHist2d((TH2D*)estimate1,"#rho","estimate");
	getShadedAreaFromTF1(fx1,fx3,gr);
  fx2->SetLineColor(kBlack);
  fx2->Draw("same");
	c5->SaveAs("rho_vs_puestimate.eps");
}

int make_radvsres_plots() {

	const double ptCut(30);
	const double drCut(0.25);
	//radius value s
	double radius[] = {20,30,40,50,60,70,80};
	//eta bins 
	double etaBins[] = {0.0,1.3};
	const unsigned numRad(sizeof(radius)/sizeof(radius[0]));
	const unsigned numEta ( sizeof(etaBins)/sizeof(etaBins[0])) ;
	//input file
	const string dir("/uscms_data/d3/tlibeiro/TuneFFTJet/TTJets_MassiveBinDECAY_TuneZ2star8TeV53V7Cv1_MultiRad");
 	//declare and draw  hists
		TH1D* hW1_fft[numEta-1][numRad];
		TH1D* ht1_fft[numEta-1][numRad];
		TCanvas* cw[numEta-1][numRad];
		TCanvas* ct[numEta-1][numRad];
	TFile* f1 = new TFile((dir+"/Analyzer3_merged.root").c_str(),"read");
  //resolution arrays
	double wres[numEta-1][numRad], tres[numEta-1][numRad];
	double werr[numEta-1][numRad], terr[numEta-1][numRad],
				 rerr[numEta-1][numRad];
  const double c(0.94247451);
	for(unsigned nrad(0);nrad<numRad;++nrad) {
    ostringstream inputNtuple; 
		inputNtuple<<"Analyzer3R"<<radius[nrad]<<"/Wt_ntup";
		TNtuple* ntup1 = (TNtuple*)f1->Get(inputNtuple.str().c_str());
			for(unsigned i(0);i<numEta-1;++i)
		{
			ostringstream oss_f,oss_a;
			ostringstream cut_f,cut_a,draw_f,draw_a;
			oss_f<<"hWf_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]<<"_rad_"<<radius[nrad]; // histos 
			draw_f<<"W1_f/W1_p>>"<<oss_f.str(); //draw strings
			cut_f<<"maxdrfj<"<<drCut<<"&&w1dau1_pt_fj>"<<ptCut
				<<"&&w1dau2_pt_fj>"<<ptCut
				<<"&&fabs(w1dau1_eta_fj)>"<<etaBins[i]
				<<"&&fabs(w1dau2_eta_fj)<"<<etaBins[i+1];  // cut string 
			//histos 
			hW1_fft[i][nrad] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0,5);
			ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
			//top 
			oss_f.str(""); draw_f.str(""); 
			//Add b cuts
			cut_f<<"&&b1_pt_fj>"<<ptCut
				<<"&&fabs(b1_eta_fj)>"<<etaBins[i]
				<<"&&fabs(b1_eta_fj)<"<<etaBins[i+1];
			oss_f<<"htf_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]<<"_rad_"<<radius[nrad]; //histos 
			draw_f<<"t1_f/t1_p>>"<<oss_f.str(); //draw strings
			//histos
			ht1_fft[i][nrad] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0,5);
			ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
      //get resolutions
      //set width_w1f [ expr {($p84-$p15)*0.5/$median} ]
			wres[i][nrad] = getRMS(*hW1_fft[i][nrad])/getMedian(*hW1_fft[i][nrad]);
			tres[i][nrad] = getRMS(*ht1_fft[i][nrad])/getMedian(*ht1_fft[i][nrad]);
      werr[i][nrad] = wres[i][nrad]*c/hW1_fft[i][nrad]->GetEntries();
      terr[i][nrad] = tres[i][nrad]*c/ht1_fft[i][nrad]->GetEntries();
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
			//top
			cn.str(""); cn<<"ct"<<i<<"_rad_"<<radius[nrad];
			ct[i][nrad] = new TCanvas(cn.str().c_str(),cn.str().c_str(),600,600);
			DrawHist1dRes(ht1_fft[i][nrad],"m^{t}(jjb)",etaLim);
			//save histos
			cn.str(""); 
			cn<<"Wres_eta_"<<etaBins[i]<<"_rad_"<<radius[nrad];
			cw[i][nrad]->SaveAs((cn.str()+".png").c_str());
			cn.str("");
			cn<<"tres_eta_"<<etaBins[i]<<"_rad_"<<radius[nrad];
			ct[i][nrad]->SaveAs((cn.str()+".png").c_str());
		}
	}//radius loop

	TCanvas*cres = new TCanvas("cres","cres",10,10,600,600);
	TGraphErrors* gr_W = new TGraphErrors(numRad,radius,wres[0],rerr[0],werr[0]);
	TGraphErrors* gr_t = new TGraphErrors(numRad,radius,tres[0],rerr[0],terr[0]);
  //set options
  gr_W->SetLineColor(1); gr_W->SetMarkerColor(2); gr_W->SetMarkerStyle(21);
  gr_t->SetLineColor(1); gr_t->SetMarkerColor(4); gr_t->SetMarkerStyle(21);
  //draw multigraph
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Resolution Graph");
	mg->Add(gr_W);
	mg->Add(gr_t);
	mg->Draw("APZC");
  //save the resolution plot
	cres->SaveAs("RadiusvsResolution.png");
}

int make_fjak_plots() {

	const double ptCut(30);
	const double drCut(0.25);
	//input file
	const string dir("/uscms_data/d3/tlibeiro/TuneFFTJet/TTJets_MassiveBinDECAY_TuneZ2star_8TeV_53_V7C_v1_TestForBoost");
	TFile* f1 = new TFile((dir+"/Analyzer3_merged.root").c_str(),"read");
	TNtuple* ntup1 = (TTree*)f1->Get("analyzer3_test/Wt_ntup");
	//eta bins 
	double etaBins[] = {0.0,1.3,3.0,5.0};
	const unsigned numEta ( sizeof(etaBins)/sizeof(etaBins[0])) ;
	//declare and draw  hists
	TH1D* hW1_fft[numEta-1];
	TH1D* ht1_fft[numEta-1];
	TH1D* hW1_ak5[numEta-1];
	TH1D* ht1_ak5[numEta-1];
	TCanvas* cw[numEta-1];
	TCanvas* ct[numEta-1];
	for(unsigned i(0);i<numEta-1;++i)
	{
		ostringstream oss_f,oss_a;
		ostringstream cut_f,cut_a,draw_f,draw_a;
		oss_f<<"hWf_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]; // histos 
		oss_a<<"hWa_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]; // histos
		draw_f<<"W1_f>>"<<oss_f.str(); //draw strings
		draw_a<<"W1_a>>"<<oss_a.str(); //draw strings 
		cut_f<<"maxdrfj<"<<drCut<<"&&w1dau1_pt_fj>"<<ptCut
			<<"&&w1dau2_pt_fj>"<<ptCut
			<<"&&fabs(w1dau1_eta_fj)>"<<etaBins[i]
			<<"&&fabs(w1dau2_eta_fj)<"<<etaBins[i+1];  // cut string 
		cut_a<<"maxdrak<"<<drCut<<"&&w1dau1_pt_ak>"<<ptCut
			<<"&&w1dau2_pt_ak>"<<ptCut
			<<"&&fabs(w1dau1_eta_ak)>"<<etaBins[i]
			<<"&&fabs(w1dau2_eta_ak)<"<<etaBins[i+1];  // cut string 
		//histos 
		hW1_fft[i] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0, 200);
		hW1_ak5[i] = new TH1D(oss_a.str().c_str(),oss_a.str().c_str(),100, 0, 200);
		ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
		ntup1->Draw(draw_a.str().c_str(),cut_a.str().c_str()); 
		//top 
		oss_f.str("");  oss_a.str("");
		draw_f.str(""); draw_a.str("");
		//Add b cuts
		cut_f<<"&&b1_pt_fj>"<<ptCut
			<<"&&fabs(b1_eta_fj)>"<<etaBins[i]
			<<"&&fabs(b1_eta_fj)<"<<etaBins[i+1];
		cut_a<<"&&b1_pt_ak>"<<ptCut
			<<"&&fabs(b1_eta_ak)>"<<etaBins[i]
			<<"&&fabs(b1_eta_ak)<"<<etaBins[i+1];
		oss_f<<"htf_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]; //histos 
		oss_a<<"hta_eta_"<<etaBins[i]<<"_"<<etaBins[i+1]; //histos 
		draw_f<<"t1_f>>"<<oss_f.str(); //draw strings
		draw_a<<"t1_a>>"<<oss_a.str(); //draw strings 
		//histos
		ht1_fft[i] = new TH1D(oss_f.str().c_str(),oss_f.str().c_str(),100, 0, 350);
		ht1_ak5[i] = new TH1D(oss_a.str().c_str(),oss_a.str().c_str(),100, 0, 350);
		ntup1->Draw(draw_f.str().c_str(),cut_f.str().c_str()); 
		ntup1->Draw(draw_a.str().c_str(),cut_a.str().c_str()); 
	}

	//maxdr plot
	TH1D* h_maxdrfj=new TH1D("h_maxdrfj","h_maxdrfj",50,0,1);
	TH1D* h_maxdrak=new TH1D("h_maxdrak","h_maxdrak",50,0,1);
  //wmass_pt vs wmass_reconstructed
	TH2D* h_ptvsmass=new TH2D("h_ptvsmass","h_ptvsmass",100,0,200,250,0,500);
  //Draw hists cuts
	ostringstream oss_cut;
	oss_cut<<"w1dau1_pt_fj>"<<ptCut
		<<"&&w1dau2_pt_fj>"<<ptCut;
	ntup1->Draw("maxdrfj>>h_maxdrfj",oss_cut.str().c_str());oss_cut.str("");
	oss_cut<<"w1dau1_pt_ak>"<<ptCut
		<<"&&w1dau2_pt_ak>"<<ptCut;
  oss_cut<<"&&fabs(w1dau1_eta_ak)>"<<0
      <<"&&fabs(w1dau2_eta_ak)<"<<1.3;
	ntup1->Draw("maxdrak>>h_maxdrak",oss_cut.str().c_str());
	ntup1->Draw("t1_pt_p:W1_a>>h_ptvsmass",oss_cut.str().c_str());
  // Draw hists
	TCanvas* cdr = new TCanvas("cdr","cdr",600,600);
	DrawHist1d(h_maxdrfj,h_maxdrak,"Max Dr Jet-Parton");
	cdr->SaveAs("Maxdr_JetPart.eps");
	TCanvas* cptvm = new TCanvas("cptvm","cptvm",600,600);
	DrawHist2d(h_ptvsmass,"W mass ak5","W pt Parton Level");
	cptvm->SaveAs("Wmass_vs_Wpt.eps");

	//Draw histograms
	for(unsigned i(0);i<numEta-1;++i)
	{
		const float etaLim[2] ={etaBins[i],etaBins[i+1]};
		//W
		ostringstream cn;
		cn<<"cw"<<i;
		cw[i] = new TCanvas(cn.str().c_str(),cn.str().c_str(),600,600);
		DrawHist1dRes(hW1_fft[i],hW1_ak5[i],"m^{W}(jj)",etaLim);
		//top
		cn.str(""); cn<<"ct"<<i;
		ct[i] = new TCanvas(cn.str().c_str(),cn.str().c_str(),600,600);
		DrawHist1dRes(ht1_fft[i],ht1_ak5[i],"m^{t}(jjb)",etaLim);
		//save histos
		cn.str(""); 
		cn<<"W_eta_"<<etaBins[i];
		cw[i]->SaveAs((cn.str()+".eps").c_str());
		cn.str("");
		cn<<"t_eta_"<<etaBins[i];
		ct[i]->SaveAs((cn.str()+".eps").c_str());
	}
}
