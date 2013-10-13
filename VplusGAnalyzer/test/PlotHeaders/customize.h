#include<string>
#include<map>

struct Options
{
 Options():rebin(50), legendx1(0.70),legendy1(0.65),legendx2(0.90),legendy2(0.88),min(0.01),max(2700), yaxis("Events") {}
 int rebin;
 float legendx1;
 float legendy1;
 float legendx2;
 float legendy2;
 std::string xaxis;
 std::string yaxis;
 float min;
 float max;
};

std::map<std::string,Options> options;

void initOptions( )
{
Options o1;
o1.rebin=20;
options.insert(std::pair<std::string,Options>(std::string("VlightRegionHWen/SimpleJets_dPhiVlightRegionHWen"), o1)); 
}
