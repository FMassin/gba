#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <netcdfcpp.h>


typedef std::map<char, gsl_matrix *> Amps;
typedef std::map<char, NcVar *> AmpsNC;

class TData {
private:
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
	TData *_td;
	double _cov[2][2];
	double _mn[2];

public:
	void init(const std::string &filename = "./data/GbA_training.nc");
	void compute_likelihood(double *data, int nbands, float time, char cmpnt,
								int nsim, double mean[2], double cov[2][2]);
};
