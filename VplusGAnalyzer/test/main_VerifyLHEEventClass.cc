#include "VerifyLHEEventClass.h"
#include "TROOT.h"
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
  string indir = "/uscms_data/d3/tlibeiro/GenerateWGamma/CMSSW_5_3_8_patch1/src/GeneratorInterface/LHEInterface/test";
  string outdir= "/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput";
  string infile(indir+"/verifyLHE.root"), outfile(outdir+"/VerifyLHEEventClass.root");
  //signal files
  vector<string> infiles;
  infiles.push_back(infile);

  TROOT root("astring","bstring");
  root.SetBatch(kTRUE);

  cout<<"Verify LHE =========================================="<<endl;
	VerifyLHEEventClass vlhe(infiles,outfile);
  vlhe.Loop();

}

