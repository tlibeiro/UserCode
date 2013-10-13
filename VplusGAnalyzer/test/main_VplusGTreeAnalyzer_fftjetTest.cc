#include "VplusGTreeAnalyzer_fftjetTest.h"
#include "TROOT.h"
#include <string>
#include <stdlib.h>
#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
  string indir = "/uscms_data/d3/tlibeiro/VPlusGamma/condorfiles";
  string outdir= "/uscms_data/d3/tlibeiro/VPlusGamma/VplusGAnalyzerOutput";
  string infile_sg(indir+"/fftjetTest/VplusGTree_*.root"), outfile_sg(outdir+"/VplusGTreeAnalyzer_fftjetTest.root");
  //signal files
  vector<string> infiles_sg;
  infiles_sg.push_back(infile_sg);

  TROOT root("astring","bstring");
  root.SetBatch(kTRUE);
  
  VplusGTreeAnalyzer_fftjetTest test(infiles_sg,"VplusGTreeAnalyzer_fftjetTest_100Karis.root");
  test.Loop();

}

