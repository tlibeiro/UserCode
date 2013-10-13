#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <iomanip>      // std::setw

struct legendVars {
legendVars(string l,float xmi,float ymi,
										float xma,float yma)
:label(l),xmin(xmi),ymin(ymi),xmax(xma),ymax(yma) {}
string label;
float xmin;
float xmax;
float ymin;
float ymax;
	};

struct Variable {
	Variable(std::string n, std::string xlab,unsigned bins,float xmi, float xma)
	:name(n), xlabel(xlab),nbins(bins),xmin(xmi), xmax(xma) {} 
  
  std::string name; 
  std::string xlabel;
  float nbins;
	float xmin;
	float xmax;
	};
  
  struct Sample {
    Sample(unsigned evProcessed,float xs,std::string n,std::string filetmp, int c, bool isdata,
					float datalumi=-1.)
      :nevents(evProcessed),xsec(xs),luminosity(datalumi),name(n),filename(filetmp),color(c),
			 data(isdata),f(0){}
		Sample()
      :nevents(0),xsec(0),luminosity(0),color(0),data(false),f(0) {}

    float lumi() {  if(data) return luminosity; else {return numberOfEvents()/xsec;} }
    float scale(float l) {return l/lumi();}
    TFile * file() { if(f) return f; else return f=new TFile(filename.c_str(),"read");}
    float numberOfEvents() 
    {
      if(nevents !=-1) return nevents;
 			else {cout<<"\nNo NumEvents found:"<<name<<'\n';}
    }  

    void dump(float l)
    {
       std::cout << name <<setw(25)<<  "\t& " << xsec << "\t& " <<  lumi()/1000 << "/fb \t& " << scale(l) << std::endl;
    }

    float nevents;
    float xsec;
    float luminosity;
    std::string name;
    std::string filename;
    int color;
    bool data;
    TFile * f;
  }; 

std::vector<Sample> samples()
{
  std::vector<Sample> s;
  //
//  const unsigned 
//								// GJets200Events (7907579),	
//								 GJets200Events (8068688),	//B2G  8158022
//								 GJets400Events(38601359),  //B2G 39084597
//                 WAEvents(94616), // B2G 95148
//                 DataEvents(4664618), // B2G ??
//								 VerifyLHEEvents(1000000)
//					;
 const unsigned 
								// GJets200Events (7907579),	
								 GJets200Events (8158022),	//B2G  8158022
								 GJets400Events(39084597),  //B2G 39084597
                 WAEvents         (95148), // B2G 95148
                 DataEvents     (0), // B2G RunA+B 4942917
								 VerifyLHEEvents(1000000)
					;
  const double 
 								 DataLumi((4.429+0.749439)*1000)
				  ;
  const string indir("/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput");
  s.push_back(Sample(GJets200Events,960.5,"GJets_HT_200_400",indir+"/VplusGTreeAnalyzer_GJets200.root", kGreen+1 , false));
  s.push_back(Sample(GJets400Events,107.5,"GJets_HT_400_Inf",indir+"/VplusGTreeAnalyzer_GJets400.root", kRed , false));
  s.push_back(Sample(WAEvents      ,0.0786,"WA_LO"          ,indir+"/VplusGTreeAnalyzer_WA.root"      , kBlue, false));
  s.push_back(Sample(DataEvents    ,1000  ,"data"           ,indir+"/VplusGTreeAnalyzer_Data.root"    ,kBlack, true, DataLumi));

////for test
//  s.push_back(Sample(1000000,0.0786,"WA_merged_w_lo","/uscms_data/d3/tlibeiro/GenerateWGamma/CMSSW_5_3_8_patch1/src/GeneratorInterface/LHEInterface/test/verifyLHE.root", kRed , false));
//  s.push_back(Sample(1000,1000,"data","/uscms_data/d3/tlibeiro/GenerateWGamma/CMSSW_5_3_8_patch1/src/GeneratorInterface/LHEInterface/test/verifyLHE.root",kRed , true, 19042.0));
//  //
  return s;
}
