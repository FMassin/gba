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
	TData();

	~TData();

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
	GbA(int nb, int ns);

	~GbA();

	// Return the magnitude samples. This is done by converting the values to
	// a vector to return both values and size
	std::vector<double> get_nms();

	// Return the distance samples. This is done by converting the values to
	// a vector to return both values and size
	std::vector<double> get_nrs();

	// Change the default values for magnitude and distance samples
	void set_samples(double *msamples, int nms, double *rsamples, int nrs);

	// Get the path to the default training data
	std::string get_default_tdata();

	// Load the training data
	void init(const std::string &filename = "");

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
