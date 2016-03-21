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

void GbA::compute_likelihood(double *data, int nbands, float time, char cmpnt,
							int nsim, double mean[2], double cov[2][2]){
	int i, j, k, timeidx=0;
	double timemin=std::numeric_limits<double>::max();
	gsl_vector *lsq = gsl_vector_alloc (_td->get_noftraces());
	gsl_vector *input = gsl_vector_alloc (nbands);
	gsl_vector *trace = gsl_vector_alloc (nbands);
	gsl_matrix *tdata;
	double *mags, *dists;
	double misfit_l2;
	size_t *indices = new size_t[nsim];
	mags = new double[nsim];
	dists = new double[nsim];
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
	gsl_sort_vector_smallest_index(indices, nsim,lsq);
	mean[0] = 0;
	mean[1] = 0;
	for(k=0; k<nsim; k++){
		dists[k] = std::log10(gsl_vector_get(_td->get_dist(), indices[k]));
		mags[k] = gsl_vector_get(_td->get_mag(), indices[k]);
	}
	mean[0] = gsl_stats_mean(dists,1,nsim);
	mean[1] = gsl_stats_mean(mags,1,nsim);
	cov[0][0] = gsl_stats_variance(dists, 1, nsim);
	cov[1][1] = gsl_stats_variance(mags, 1, nsim);
	cov[0][1] = gsl_stats_covariance_m(dists,1,mags,1,nsim,mean[0],mean[1]);
	cov[1][0] = cov[0][1];
#ifdef DEBUG
		std::cout << "Distance: " << mean[0] << "; Magnitude: " << mean[1] << std::endl;
		std::cout << "Covariance matrix:" << std::endl;
		for(i=0; i<2; i++){
			for(j=0;j<2;j++){
				printf("%d, %d, %g\n",i,j,cov[i][j]);
			}
		}
#endif

}
