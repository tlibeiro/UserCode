#include <TH1F.h>
#include <iomanip.h>

double getRMS(const TH1F& h);
double getRMS(const TH1F& h)
{
	const unsigned numprob(2);
	double probsum[numprob];
        const double q[numprob] = {.1583,.8417};
	h.GetQuantiles(numprob,probsum,q);

        const double rms ((probsum[1]-probsum[0])/2);
	return rms;
}
inline double getRMS(const TH1D& h) {
	const unsigned numprob(2);
	double probsum[numprob];
        const double q[numprob] = {.1583,.8417};
	h.GetQuantiles(numprob,probsum,q);

        const double rms ((probsum[1]-probsum[0])/2);
	return rms;

};


inline double getMedian(const TH1F& h)
{
	const unsigned numprob(1);
	double probsum[numprob];
        const double q[numprob] = {.5};
	h.GetQuantiles(numprob,probsum,q);

	return probsum[0];
};

inline double getMedian(const TH1D& h)
{
	const unsigned numprob(1);
	double probsum[numprob];
        const double q[numprob] = {.5};
	h.GetQuantiles(numprob,probsum,q);

	return probsum[0];
};


inline double getQuantile(const TH1D& h,double quantile)
{
	const unsigned numprob(1);
	double probsum[numprob];
        const double q[numprob] = {quantile};
	h.GetQuantiles(numprob,probsum,q);
	//cout<<"getQuantile ouput "<<scientific<<setprecision(10)<<probsum[0]<<endl;
	return probsum[0];
};
inline double getdrFFTcut(const TNtuple& ntup,const string ak_cut ,const string other_fft_cuts(""))
{
	cout.precision(10);
	TH1D histo1("histo1","histo1",2000, 0, 20);
	TH1D histo2("histo2","histo2",2000, 0, 20);
	TH1D test("test","test",100, 0, 10);
	const double initialdrcut(0.25);
	char tmpstr[500]; std::string drcutstr(""), cut("");

	ntup.Draw("maxdrfj>>histo1",ak_cut.c_str());
	ntup.Draw("maxdrfj>>histo2",other_fft_cuts.c_str());
	double akevents = histo1.GetEntries(), ffevents = histo2.GetEntries();
	double newfraction = akevents/ffevents;	
	if (newfraction>1.0) newfraction=1.0; 
	double drfftcut = getQuantile(histo2,newfraction);
	sprintf(tmpstr,"%.10f",drfftcut);drcutstr.assign(tmpstr);
	cut = other_fft_cuts+"&& maxdrfj<"+drcutstr;
	ntup.Draw("maxdrfj>>histo2",cut.c_str());

	//double drcut(0), drstep(1e-4); 
	//unsigned count(0), maxcount(1e4);
	//while(akevents!=ffevents && count<maxcount) {
	//	if(akevents>ffevents)
	//		drcut = drfftcut + drstep; 
	//	else if(akevents<ffevents)
	//		drcut = drfftcut-drstep;
	//	sprintf(tmpstr,"%.10f",drcut);drcutstr.assign(tmpstr);
	//	cut = other_fft_cuts+"&& maxdrfj<"+drcutstr;

	//	ntup.Draw("maxdrfj>>histo2",cut.c_str());
	//	ffevents = histo2.GetEntries();
	//	newfraction = akevents/ffevents;
	//	cout<<"Selected fft  events "<<ffevents <<" Ak5 Entries "<<histo1.GetEntries()<<" new fraction "<<newfraction<<" drcut "<<drcut<<" count "<<count<<endl;
	//	++count;
	//}
	cout<<"drfftcut "<<drfftcut<<endl;


	return drfftcut;


};

inline TH1F& scaleData(const TH1F& h, const double setMedianValue)
{
	const median = getmedian(h);
	float* histocontent = h.GetArray();
	const unsigned nument = h->GetEntries(); 
	const unsigned sizearray = (sizeof (histocontent)/sizeof (histocontent[0]));
	for(unsigned i(0); i<sizearray; ++i)
		histocontent*=(setMedianValue/median);
}

