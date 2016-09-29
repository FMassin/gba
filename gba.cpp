#include "gba.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_randist.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <pthread.h>

TData::TData(){
	_nid = 0;
	_m_gsl = 0;
	_r_gsl = 0;
	_t_gsl = 0;
	_amps[TData::vertical] = 0;
	_amps[TData::horizontal] = 0;
}

TData::~TData(){
	delete _nid;
	gsl_vector_free(_m_gsl);
	gsl_vector_free(_r_gsl);
	gsl_vector_free(_t_gsl);
	gsl_matrix_free(_amps[TData::vertical]);
	gsl_matrix_free(_amps[TData::horizontal]);
}

bool TData::load(std::string filename) {
#ifdef DEBUG
		std::cout << "Reading training data set from: " << filename << std::endl;
#endif
	_start.push_back(0);
	_count.push_back(0);
	try{
	_nid = new netCDF::NcFile(filename, netCDF::NcFile::read);
	_dim_nt = _nid->getDim("traces");
	_ntraces = _dim_nt.getSize();
	_dim_t = _nid->getDim("time");
	_ntimes = _dim_t.getSize();
	_dim_f = _nid->getDim("filter");
	_nbands = _dim_f.getSize();
	_m_gsl = gsl_vector_alloc (_dim_nt.getSize());
	_r_gsl = gsl_vector_alloc (_dim_nt.getSize());
	_t_gsl = gsl_vector_alloc (_dim_t.getSize());
	_var_mag = _nid->getVar("magnitude");
	_var_dist = _nid->getVar("epicdist");
	_var_time = _nid->getVar("time");
	_count[0] = _dim_nt.getSize();
	_var_mag.getVar(_start, _count, _m_gsl->data);
	_var_dist.getVar(_start, _count, _r_gsl->data);
	_count[0] = _dim_t.getSize();
	_var_time.getVar(_start, _count, _t_gsl->data);
	_amps_var[vertical] = _nid->getVar("z");
	_amps_var[horizontal] = _nid->getVar("h");

#ifdef DEBUG
	std::cout << "Number of traces: " << _dim_nt.getSize() << std::endl;
	std::cout << "Number of filters: " << _dim_f.getSize() << std::endl;
	std::cout << "Number of timesteps: " << _dim_t.getSize() << std::endl;
#endif

	_amps[vertical] = gsl_matrix_alloc(_dim_nt.getSize(),_dim_f.getSize());
	_amps[horizontal] = gsl_matrix_alloc(_dim_nt.getSize(),_dim_f.getSize());

	// Preload data for first time step
	_start.push_back(0);
	_start.push_back(0);
	_count[0] = 1;
	_count.push_back(_dim_nt.getSize());
	_count.push_back(_dim_f.getSize());
	_amps_var[vertical].getVar(_start,_count,_amps[vertical]->data);
	_amps_var[horizontal].getVar(_start,_count,_amps[horizontal]->data);
	_currenttime_idx = 1;
	}catch(netCDF::exceptions::NcException &e){
		e.what();
		return false;
	}
	return true;
}

gsl_matrix * TData::get_amps(int timeidx, TData::Component cmpnt){
	if (timeidx != _currenttime_idx){
		_start[0] = timeidx;
		try{
			_amps_var[cmpnt].getVar(_start,_count,_amps[cmpnt]->data);
		}catch(netCDF::exceptions::NcException &e){
			e.what();
		}
		return _amps[cmpnt];
	}else{
		return _amps[cmpnt];
	}
}

int TData::get_noftraces(){
	return _ntraces;
}

int TData::get_nofbands(){
	return _nbands;
}

int TData::get_noftimes(){
	return _ntimes;
}

gsl_vector * TData::get_dist(){
	return _r_gsl;
}

gsl_vector * TData::get_mag(){
	return _m_gsl;
}

gsl_vector * TData::get_times(){
	return _t_gsl;
}

GbA::GbA(int nb, int ns){
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


GbA::~GbA(){
	delete[] _msamples;
	delete[] _rsamples;
	gsl_vector_free(_mags);
	gsl_vector_free(_dists);
	delete _td;
}

void GbA::init(const std::string &filename){
	_td = new TData;
	_td->load(filename);
}

void GbA::set_samples(double *msamples, int nms, double *rsamples, int nrs){
	delete[] _msamples;
	delete[] _rsamples;
	_msamples = new double[nms];
	_rsamples = new double[nrs];
	memcpy(_msamples,msamples,nms*sizeof(double));
	memcpy(_rsamples,rsamples,nrs*sizeof(double));
	_nms = nms;
	_nrs = nrs;
}

std::vector<double> GbA::get_nms(){
	std::vector<double> m(_nms);
	for(int i=0;i<_nms;i++){
		m[i] = _msamples[i];
	}
	return m;
}

std::vector<double> GbA::get_nrs(){
	std::vector<double> r(_nrs);
	for(int i=0;i<_nrs;i++){
		r[i] = _rsamples[i];
	}
	return r;
}

void GbA::process(double *data, int nbands, float time, int cmpnt){

	assert(nbands == _nbands);
	pthread_mutex_lock(&_process_lock);
	int i, j, k, timeidx=0;
	double timemin=std::numeric_limits<double>::max();
	gsl_vector *lsq = gsl_vector_alloc (_td->get_noftraces());
	gsl_vector *input = gsl_vector_alloc (nbands);
	gsl_vector *trace = gsl_vector_alloc (nbands);
	gsl_matrix *tdata;
	double misfit_l2, terror;
	size_t *indices = new size_t[_nsim];
	TData::Component _c = static_cast<TData::Component>(cmpnt);
	_status[_c] = true;
	// Find time index
	for(i=0;i<_td->get_noftimes();i++){
		terror = gsl_vector_get(_td->get_times(),i) - time;
		if(abs(terror)< timemin){
			timeidx = i;
			timemin = abs(terror);
		}
	}

	tdata = _td->get_amps(timeidx,_c);


	for (i=0; i<nbands; i++){
			gsl_vector_set(input, i, std::log10(data[i]));
	}
	for (j=0; j<_td->get_noftraces(); j++){
		gsl_matrix_get_row(trace, tdata, j);
		gsl_vector_sub(trace,input);
		misfit_l2 = gsl_blas_dnrm2(trace);
		gsl_vector_set(lsq, j, misfit_l2*misfit_l2);
	}

	gsl_sort_vector_smallest_index(indices, _nsim,lsq);

	for(k=0; k<_nsim; k++){
		gsl_vector_set(&_r[_c].vector,k, std::log10(gsl_vector_get(_td->get_dist(), indices[k])));
		gsl_vector_set(&_m[_c].vector,k, gsl_vector_get(_td->get_mag(), indices[k]));
	}

	if (_status[TData::vertical] && _status[TData::horizontal]){
		_mv = gsl_vector_subvector(_mags,0,_mags->size);
		_rv = gsl_vector_subvector(_dists,0,_dists->size);
	} else if(_status[TData::vertical]){
		_mv = gsl_vector_subvector(_mags,0,_nsim);
		_rv = gsl_vector_subvector(_dists,0,_nsim);
	} else if(_status[TData::horizontal]){
		_mv = gsl_vector_subvector(_mags,_nsim,_nsim);
		_rv = gsl_vector_subvector(_dists,_nsim,_nsim);
	}

	_mn[0] = gsl_stats_mean(_rv.vector.data,1,_rv.vector.size);
	_mn[1] = gsl_stats_mean(_mv.vector.data,1,_mv.vector.size);
	_cov[0][0] = gsl_stats_variance(_rv.vector.data, 1, _rv.vector.size);
	_cov[1][1] = gsl_stats_variance(_mv.vector.data, 1, _mv.vector.size);
	_cov[0][1] = gsl_stats_covariance_m(_rv.vector.data,1,_mv.vector.data,1,
			                            _rv.vector.size,_mn[0],_mn[1]);
	_cov[1][0] = _cov[0][1];
	delete[] indices;
	pthread_mutex_unlock(&_process_lock);
}

void GbA::get_m_r(double *m, int nm, double *r, int nr){
	assert(nm == 2*_nsim);
	for(size_t i=0; i<_mv.vector.size; i++){
		m[i] = gsl_vector_get(&_mv.vector,i);
		r[i] = gsl_vector_get(&_rv.vector,i);
	}
}

void GbA::get_mean_cov(double mean[2], double cov[2][2]){
	mean[0] = _mn[0];
	mean[1] = _mn[1];
	for(int i=0;i<2;i++){
		for(int j=0;j<2;j++){
			cov[i][j] = _cov[i][j];
		}
	}
}

void GbA::get_pdf(double *pdf, int nm, int nr, double *mmarg, int nms,
		          double *rmarg, int nrs){
	assert(nms==nm && nrs==nr);
	assert(nr > 1 && nm > 1);
	assert(nms == _nms && nrs == _nrs);

	pthread_mutex_lock(&_process_lock);
	double dm, dr, evd;
	int i,j;
	double *_m_tmp, *_r_tmp, *_r1, *_m1;
	_m_tmp = new double[nm];
	_r_tmp = new double[nr];
	_r1 = new double[nr];
	_m1 = new double[nm];
	memcpy(_r1,_rsamples,nr*sizeof(double));
	memcpy(_m1,_msamples,nm*sizeof(double));

	// assume regular samples
	dr = _rsamples[1] - _rsamples[0];
	dm = _msamples[1] - _msamples[0];

	// Compute pdf
	for(i=0; i<nm; i++){
		for(j=0; j<nr; j++){
			pdf[i*nr+j] = gsl_ran_bivariate_gaussian_pdf(_msamples[i]-_mn[1], _rsamples[j]-_mn[0],
							std::sqrt(_cov[1][1]), std::sqrt(_cov[0][0]),
							_cov[0][1]/(std::sqrt(_cov[0][0])*std::sqrt(_cov[1][1])));
		}
	}

	for(i=0;i<nm;i++){
		_m_tmp[i] = 0;
		for(j=0; j<nr-1; j++){
			_m_tmp[i] = _m_tmp[i] + (pdf[i*nr+j] + pdf[i*nr+j+1])/2.*dr;
		}
	}
	evd = 0;
	for(i=0;i<nm-1;i++)
		evd = evd + (_m_tmp[i] + _m_tmp[i+1])/2.*dm;
	for(i=0;i<nm;i++)
		_m_tmp[i] = _m_tmp[i] / evd;

	for(j=0;j<nr;j++){
		_r_tmp[j] = 0;
		for(i=0; i<nm-1; i++){
			_r_tmp[j] = _r_tmp[j] + (pdf[i*nr+j] + pdf[(i+1)*nr+j])/2.*dm;
		}
	}
	evd=0;
	for(i=0;i<nr-1;i++)
		evd = evd + (_r_tmp[i] + _r_tmp[i+1])/2.*dr;
	for(i=0;i<nr;i++)
		_r_tmp[i] = _r_tmp[i] / evd;

	memcpy(mmarg,_m_tmp,nm*sizeof(double));
	memcpy(rmarg,_r_tmp,nr*sizeof(double));
	delete[] _r1;
	delete[] _m1;
	delete[] _m_tmp;
	delete[] _r_tmp;
	pthread_mutex_unlock(&_process_lock);
}
