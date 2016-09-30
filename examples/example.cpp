#include "gba.h"

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
	// The idea is that the following is run everytime there is a new
	// filterbank observation available.
	g.process(data,nbands,3.0,TData::vertical);
	g.get_m_r(mags, 2*nsim, r, 2*nsim);
	g.get_mean_cov(mean, cov);
	return 0;
}
