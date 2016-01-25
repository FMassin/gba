#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <cmath>

typedef std::vector<double> row;
typedef std::vector<row> mat ;
typedef std::map<std::string, std::string> datafiles;

class TData {
private:
	mat _tdata;
	row _m, _r;
	int _ntraces, _nbands, _ntimes;
	std::istream * _isP;
	std::filebuf * _fbP;

public:
	bool read(std::string dtype);
	bool read_traces();
	bool read_magnitudes();
	bool read_distances();
	bool load(std::string filename, std::string dtype);
	bool open(std::string filename);
	void close();
	int get_noftraces();
	int get_nofbands();

	gsl_vector * _m_gsl, *_r_gsl;
	gsl_matrix * _tdata_gsl;
};

class GbA {

private:
	TData *_td;
	double _cov[2][2];
	double _mn[2];

public:
	void init();
	void compute_likelihood(double *data, int nstats, int nbands, int nsim,
			double mean[2], double cov[2][2]);

};
