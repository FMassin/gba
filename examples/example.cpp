#include <gba.h>

#include <iostream>
#include <limits>


using namespace std;

/*
 * Compute expectation, MAP, and variance from marginal pdfs
 */
void get_moments(double *pdf, vector<double> samples,
		         double &mean, double &max, double &sigma){
	int i, imax=0, ndata;
	double pdfmax = std::numeric_limits<double>::min();
	double tmp1, tmp2, ds;

	ndata = samples.size();
	// Assuming equidistant samples
	ds = samples[1]-samples[0];

	// normalize pdf
	double evd = 0;
	for(i=0;i<ndata-1;i++)
		evd = evd + (pdf[i] + pdf[i+1])/2.*ds;
	for(i=0;i<ndata;i++)
		pdf[i] = pdf[i] / evd;

	mean = 0;
	for(i=0; i<ndata-1; i++){
		mean += (pdf[i]*samples[i] + pdf[i+1]*samples[i+1])/2.*ds;
		if (pdf[i] > pdfmax){
			pdfmax = pdf[i];
			imax = i;
		}
	}
	if(pdf[ndata-1] > pdfmax)
		imax = ndata-1;
	max = samples[imax];

	sigma=0;
	for(i=0; i<ndata; i++){
		tmp1 = (pdf[i]*(samples[i]-mean)*(samples[i]-mean))/2.;
		tmp2 = (pdf[i+1]*(samples[i+1]-mean)*(samples[i+1]-mean))/2.;
		sigma += (tmp1+tmp2)*ds;
	}
}

int main(){
	// Filterbank observation in 9 frequency bands
	double data[] = {0.00038221, 0.00067219, 0.00056877, 0.0014739, 0.0022047,
				   0.0031288, 0.0015537, 0.00038673, 0.00012311};

	// Number of frequency bands
	int nbands=9;

	// Number of most similar training data traces
	int nsim=30;

	// Arrays for the magnitudes and distances of the most similar traces
	// Note that they have to be of length 2*nsim as the hold the values for
	// both the vertical and the horizontal component. The values for the
	// vertical component have index 0:nsim-1 and for the horizontal
	// component index nsim::
	double mags[2*nsim], r[2*nsim];

	// Arrays for the mean and covariance matrix of magnitude and distance
	double mean[2], cov[2][2];

	GbA g(nbands,nsim);
	g.init();

	// The idea is that the following is run every time there is a new
	// filterbank observation available.
	g.process(data,nbands,3.0,TData::vertical);

	// The results can be queried anytime independently of the incoming data
	g.get_m_r(mags, 2*nsim, r, 2*nsim);
	g.get_mean_cov(mean, cov);
	cout << "Mean log distance: " << mean[0] << endl;
	cout << "Mean magnitude: " << mean[1] << endl;
	cout << "Variance log distance: " << cov[0][0] << endl;
	cout << "Variance magnitude: " << cov[1][1] << endl;
	cout << "Covariance: " << cov[0][1] <<endl;


	// Get the array of magnitude and distance samples
	vector<double> msmp = g.get_nms();
	vector<double> rsmp = g.get_nrs();
	double *m_marg, *r_marg, *pdf;
	double mbar,mhat,msig,rbar,rhat,rsig;
	m_marg = new double[msmp.size()];
	r_marg = new double[rsmp.size()];
	pdf = new double[msmp.size()*rsmp.size()];
	g.get_pdf(pdf,msmp.size(),rsmp.size(),&m_marg[0],msmp.size(),
			&r_marg[0],rsmp.size());
	get_moments(m_marg,msmp,mbar,mhat,msig);
	get_moments(r_marg,rsmp,rbar,rhat,rsig);
	cout << "E(log distance): " << rbar << endl;
	cout << "MAP(log distance): " << rhat << endl;
	cout << "Var(log distance): " << rsig << endl;
	cout << "E(magnitude): " << mbar << endl;
	cout << "MAP(magnitude): " << mhat << endl;
	cout << "Var(magnitude): " << msig << endl;



	return 0;
}
