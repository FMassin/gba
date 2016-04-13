#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cassert>
#include <cstring>
#include <pthread.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>
#include <netcdfcpp.h>



class TData {
public:
	enum Component {
		vertical = 0,
		horizontal =1
	};
	bool load(std::string filename);
	int get_noftraces();
	int get_nofbands();
	int get_noftimes();
	gsl_matrix * get_amps(int timeidx, Component compnt);
	gsl_vector * get_dist();
	gsl_vector * get_mag();
	gsl_vector * get_times();

private:
	typedef std::map<Component, gsl_matrix *> Amps;
	typedef std::map<Component, NcVar *> AmpsNC;

	int _ntraces, _nbands, _ntimes;
	int _currenttime_idx;
	NcFile *_nid;
	NcVar *_var_zamp, *_var_hamp, *_var_mag, *_var_dist, *_var_time;
	NcDim *_dim_f, *_dim_nt, *_dim_t;
	AmpsNC _amps_var;
	Amps _amps;
	gsl_vector * _m_gsl, *_r_gsl, *_t_gsl;



};

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
	GbA(int nb, int ns):
		_nbands(nb),
		_nsim(ns){
		_status[TData::vertical] = false;
		_status[TData::horizontal] = false;
		_mags = gsl_vector_calloc(2*_nsim);
		_dists = gsl_vector_calloc(2*_nsim);
		_m[TData::vertical] = gsl_vector_subvector(_mags,0,_nsim);
		_m[TData::horizontal] = gsl_vector_subvector(_mags,_nsim,_nsim);
		_r[TData::vertical] = gsl_vector_subvector(_dists,0,_nsim);
		_r[TData::horizontal] = gsl_vector_subvector(_dists,_nsim,_nsim);
		pthread_mutex_init(&_process_lock,NULL);
		_nms = 31;
		_nrs = 21;
		_msamples = new double[_nms];
		_rsamples = new double[_nrs];
		for(int i=0;i<_nms;i++)
			_msamples[i] = 2.0 + 0.2*i;
		for(int i=0;i<_nrs;i++)
			_rsamples[i] = 0.0 + 0.1*i;
	};
	void set_samples(double *msamples, int nms, double *rsamples, int nrs);
	void init(const std::string &filename = "./data/GbA_training.nc");
	void process(double *data, int nbands, float time, int cmpnt);
	void get_m_r(double *m, int nm, double *r, int nr);
	void get_mean_cov(double mean[2], double cov[2][2]);
	void get_pdf(double *pdf, int nm, int nr, double *mmarg, int nms,
			          double *rmarg, int nrs);
};
