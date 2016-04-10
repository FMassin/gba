#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstring>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <netcdfcpp.h>



class TData {
private:
	typedef std::map<char, gsl_matrix *> Amps;
	typedef std::map<char, NcVar *> AmpsNC;

	int _ntraces, _nbands, _ntimes;
	int _currenttime_idx;
	NcFile *_nid;
	NcVar *_var_zamp, *_var_hamp, *_var_mag, *_var_dist, *_var_time;
	NcDim *_dim_f, *_dim_nt, *_dim_t;
	AmpsNC _amps_var;
	Amps _amps;
	gsl_vector * _m_gsl, *_r_gsl, *_t_gsl;


public:
	bool load(std::string filename);
	int get_noftraces();
	int get_nofbands();
	int get_noftimes();
	gsl_matrix * get_amps(int timeidx, char compnt);
	gsl_vector * get_dist();
	gsl_vector * get_mag();
	gsl_vector * get_times();

};

class GbA {

private:
	typedef std::map<char, gsl_vector_view> MRresult;
	typedef std::map<char, bool> UseCmp;

	TData *_td;
	int _nbands,_nsim;
	double _cov[2][2];
	double _mn[2];
	MRresult _m, _r;
	gsl_vector *_mags, *_dists;
	gsl_vector_view _mv, _rv;
	UseCmp _status;

public:
	GbA(int nb, int ns):
		_nbands(nb),
		_nsim(ns){
		_status['z'] = false;
		_status['h'] = false;
		_mags = gsl_vector_calloc(2*_nsim);
		_dists = gsl_vector_calloc(2*_nsim);
		_m['z'] = gsl_vector_subvector(_mags,0,_nsim);
		_m['h'] = gsl_vector_subvector(_mags,_nsim,_nsim);
		_r['z'] = gsl_vector_subvector(_dists,0,_nsim);
		_r['h'] = gsl_vector_subvector(_dists,_nsim,_nsim);

	};
	void init(const std::string &filename = "./data/GbA_training.nc");
	void process(double *data, int nbands, float time, char cmpnt);
	void get_m_r(double *m, int nm, double *r, int nr);
	void get_mean_cov(double mean[2], double cov[2][2]);
	void get_pdf(double *pdf, int nm, int nr, double *msamples, int nms,
			          double *rsamples, int nrs);
};
