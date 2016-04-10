#include "gba.h"
#include <iostream>


bool TData::load(std::string filename) {
#ifdef DEBUG
		std::cout << "Reading training data set from: " << filename << std::endl;
#endif
	_nid = new NcFile(filename.c_str());
	if(!_nid->is_valid())
		return false;
	_dim_nt = _nid->get_dim("traces");
	_ntraces = _dim_nt->size();
	_dim_t = _nid->get_dim("time");
	_ntimes = _dim_t->size();
	_dim_f = _nid->get_dim("filter");
	_nbands = _dim_f->size();
	_m_gsl = gsl_vector_alloc (_dim_nt->size());
	_r_gsl = gsl_vector_alloc (_dim_nt->size());
	_t_gsl = gsl_vector_alloc (_dim_t->size());
	_var_mag = _nid->get_var("magnitude");
	_var_dist = _nid->get_var("epicdist");
	_var_time = _nid->get_var("time");
	_var_mag->get(_m_gsl->data,_dim_nt->size());
	_var_dist->get(_r_gsl->data,_dim_nt->size());
	_var_time->get(_t_gsl->data,_dim_t->size());
	_amps_var['z'] = _nid->get_var("z");
	_amps_var['h'] = _nid->get_var("h");

#ifdef DEBUG
	std::cout << "Number of traces: " << _dim_nt->size() << std::endl;
	std::cout << "Number of filters: " << _dim_f->size() << std::endl;
	std::cout << "Number of timesteps: " << _dim_t->size() << std::endl;
#endif

	_amps['z'] = gsl_matrix_alloc(_dim_nt->size(),_dim_f->size());
	_amps['h'] = gsl_matrix_alloc(_dim_nt->size(),_dim_f->size());

	// Preload data for first time step
	_amps_var['z']->get(_amps['z']->data,1,_dim_nt->size(),_dim_f->size());
	_amps_var['h']->get(_amps['h']->data,1,_dim_nt->size(),_dim_f->size());
	_currenttime_idx = 1;
	return true; // if successful
}

gsl_matrix * TData::get_amps(int timeidx, char cmpnt){
	if (timeidx != _currenttime_idx){
		_amps_var[cmpnt]->set_cur(timeidx,0,0);
		_amps_var[cmpnt]->get(_amps[cmpnt]->data,1,_dim_nt->size(),_dim_f->size());
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

void GbA::init(const std::string &filename){
	_td = new TData;
	_td->load(filename);
}

void GbA::process(double *data, int nbands, float time, char cmpnt){

	assert(nbands == _nbands);

	int i, j, k, timeidx=0;
	double timemin=std::numeric_limits<double>::max();
	gsl_vector *lsq = gsl_vector_alloc (_td->get_noftraces());
	gsl_vector *input = gsl_vector_alloc (nbands);
	gsl_vector *trace = gsl_vector_alloc (nbands);
	gsl_matrix *tdata;
	double misfit_l2;
	size_t *indices = new size_t[_nsim];

	_status[cmpnt] = true;
	// Find time index
	for(i=0;i<_td->get_noftimes();i++){
		if(abs(gsl_vector_get(_td->get_times(),i) - time) < timemin){
			timeidx = i;
			timemin = abs(gsl_vector_get(_td->get_times(),i) - time);
		}
	}

	tdata = _td->get_amps(timeidx,cmpnt);


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
		gsl_vector_set(&_r[cmpnt].vector,k, std::log10(gsl_vector_get(_td->get_dist(), indices[k])));
		gsl_vector_set(&_m[cmpnt].vector,k, gsl_vector_get(_td->get_mag(), indices[k]));
	}

	if (_status['z'] && _status['h']){
		_mv = gsl_vector_subvector(_mags,0,_mags->size);
		_rv = gsl_vector_subvector(_dists,0,_dists->size);
	} else if(_status['z']){
		_mv = gsl_vector_subvector(_mags,0,_nsim);
		_rv = gsl_vector_subvector(_dists,0,_nsim);
	} else if(_status['h']){
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

void GbA::get_pdf(double *pdf, int nm, int nr, double *msamples, int nms,
		          double *rsamples, int nrs){
	assert(nms==nm && nrs==nr);
	assert(nr > 1 && nm > 1);

	double dm, dr, evd;
	int i,j;
	double *_m_tmp, *_r_tmp, *_r1, *_m1;
	_m_tmp = new double[nm];
	_r_tmp = new double[nr];
	_r1 = new double[nr];
	_m1 = new double[nm];
	memcpy(_r1,rsamples,nr*sizeof(double));
	memcpy(_m1,msamples,nm*sizeof(double));

	// assume regular samples
	dr = rsamples[1] - rsamples[0];
	dm = msamples[1] - msamples[0];

	// Compute pdf
	for(i=0; i<nm; i++){
		for(j=0; j<nr; j++){
			pdf[i*nr+j] = gsl_ran_bivariate_gaussian_pdf(msamples[i]-_mn[1], rsamples[j]-_mn[0],
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

	memcpy(msamples,_m_tmp,nm*sizeof(double));
	memcpy(rsamples,_r_tmp,nr*sizeof(double));
	delete[] _r1;
	delete[] _m1;
	delete[] _m_tmp;
	delete[] _r_tmp;
}
