#ifndef GBA_H
#define GBA_H
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <netcdf>

#include <vector>
#include <map>


/*
 * Class to load and access training data stored as NetCDF4 file.
 */
class TData {
public:
	enum Component {
		vertical = 0,
		horizontal =1
	};
	// Load data stored in filename (relative or absolute path
	bool load(std::string filename);

	// Return the number of traces in the training data
	int get_noftraces();

	// Return the number of frequency bands in the training data.
	int get_nofbands();

	// Return the number of time steps in the training data.
	int get_noftimes();

	// Return amplitudes at timestep index timeidx and for the vertical
	// or horizontal component.
	gsl_matrix * get_amps(int timeidx, Component compnt);

	// Return the epicentral distances of the training data.
	gsl_vector * get_dist();

	// Return the magnitudes of the training data.
	gsl_vector * get_mag();

	// Return the time steps of the training data.
	gsl_vector * get_times();

private:
	typedef std::map<Component, gsl_matrix *> Amps;
	typedef std::map<Component, netCDF::NcVar> AmpsNC;

	int _ntraces, _nbands, _ntimes;
	int _currenttime_idx;
	std::vector<size_t> _start, _count;
	netCDF::NcFile *_nid;
	netCDF::NcVar _var_zamp, _var_hamp, _var_mag, _var_dist, _var_time;
	netCDF::NcDim _dim_f, _dim_nt, _dim_t;
	AmpsNC _amps_var;
	Amps _amps;
	gsl_vector * _m_gsl, *_r_gsl, *_t_gsl;
};


/*
 * Class to compute the magnitude-distance probability density function (PDF)
 * based on observed filterbank values and the training data.
 */
class GbA {

private:
	typedef std::map<TData::Component, gsl_vector_view> MRresult;
	typedef std::map<TData::Component, bool> UseCmp;

	TData *_td;
	int _nbands,_nsim;
	double _cov[2][2];
	double _mn[2];
	MRresult _m, _r;
	gsl_vector *_mags, *_dists;
	gsl_vector_view _mv, _rv;
	UseCmp _status;
	double *_msamples, *_rsamples;
	int _nms, _nrs;
	pthread_mutex_t _process_lock;

public:
	GbA(int nb, int ns){
		_nbands = nb;
		_nsim = ns;
		_status[TData::vertical] = false;
		_status[TData::horizontal] = false;
		_mags = gsl_vector_calloc(2*_nsim);
		_dists = gsl_vector_calloc(2*_nsim);
		_m[TData::vertical] = gsl_vector_subvector(_mags,0,_nsim);
		_m[TData::horizontal] = gsl_vector_subvector(_mags,_nsim,_nsim);
		_r[TData::vertical] = gsl_vector_subvector(_dists,0,_nsim);
		_r[TData::horizontal] = gsl_vector_subvector(_dists,_nsim,_nsim);
		pthread_mutex_init(&_process_lock,NULL);

		/* Note that the magnitude and logarithmic distance samples are
		 * initially set to M=[2.0, 2.2, ..., 8.0] and R to [0.0, 0.1, ...,2.0]
		 */
		_nms = 31;
		_nrs = 21;
		_msamples = new double[_nms];
		_rsamples = new double[_nrs];
		for(int i=0;i<_nms;i++)
			_msamples[i] = 2.0 + 0.2*i;
		for(int i=0;i<_nrs;i++)
			_rsamples[i] = 0.0 + 0.1*i;
	};

	// Change the default values for magnitude and distance samples
	void set_samples(double *msamples, int nms, double *rsamples, int nrs);

	// Load the training data
	void init(const std::string &filename = "./data/GbA_training.nc");

	// Compute mean and covariance matrix for magnitude and epicentral distance
	// based on the given filterbank data
	void process(double *data, int nbands, float time, int cmpnt);

	// Return the magnitude and epicentral distance values from the training
	// data that were the closest match to the given filterbank observation
	void get_m_r(double *m, int nm, double *r, int nr);

	// Return the mean and covariance matrix for magnitude and epicentral
	// distance
	void get_mean_cov(double mean[2], double cov[2][2]);

	// Return the bivariate Gaussian distribution sampled at the magnitude and
	// distance samples as well as the marginal PDFs for magnitude and distance.
	void get_pdf(double *pdf, int nm, int nr, double *mmarg, int nms,
			          double *rmarg, int nrs);
};
#endif
