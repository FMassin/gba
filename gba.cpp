#include "gba.h"
#include <iostream>


bool TData::open(std::string filename) {
	_isP = NULL; //&std::cin;
	_fbP = NULL;

	_fbP = new std::filebuf(); // NEW_MEM
	if ( !_fbP->open(filename.c_str(), std::ios::in) ) {
		std::cerr << "couldn't open " << filename << std::endl;
		return false;
	}
	_isP = new std::istream(_fbP); // NEW_MEM
	return true;
}

void TData::close() {
	delete _isP;
	_fbP->close();
	delete _fbP;
}

bool TData::load(std::string filename, std::string dtype) {
	if ( !open(filename) ) {
		std::cerr << "could not open file " << filename << std::endl;
		return false;
	} else {
#ifdef DEBUG
		std::cout << "Reading: " << filename << std::endl;
#endif
		if ( !read(dtype) ) {
			std::cerr << "errors while reading or interpreting file " << filename << std::endl;
			return false;
		}
		close();
	}
	return true; // if successful
}

bool TData::read_traces(){
	row r (9);
	int linecnt, charcnt, ret;
	char buf[400];
	size_t i,j;

	for (linecnt = 0; !(*_isP).eof(); linecnt++ ) {
		(*_isP).getline(buf, 400);
		charcnt = (*_isP).gcount();
		if ( 3 > charcnt ) {
			continue; // empty line found, skip to next
		}
		if ( !(400 > charcnt) )
			return false; // line too long

		ret = sscanf(buf, "%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
				&r[0], &r[1], &r[2], &r[3], &r[4], &r[5], &r[6], &r[7], &r[8]);
		if (9 > ret)
			std::cerr << "Reading failed for line: " << buf << std::endl;
		_tdata.push_back(r);
	}
#ifdef DEBUG
	std::cout << "Number of traces in training data set: " << _tdata.size() << std::endl;
#endif
	_tdata_gsl = gsl_matrix_alloc(_tdata.size(),9);
	for(i=0; i < _tdata.size(); i++){
		for (j=0; j < 9; j++){
			gsl_matrix_set(_tdata_gsl,i,j,std::log10(_tdata[i][j]));
		}
	}
	_ntraces = _tdata.size();
	_nbands = 9;
	return true;
}

bool TData::read_magnitudes(){
	int linecnt, charcnt, ret;
	double magnitude;
	char buf[10];

	for (linecnt = 0; !(*_isP).eof(); linecnt++ ) {
		(*_isP).getline(buf, 10);
		charcnt = (*_isP).gcount();
		if ( 1 > charcnt ) {
			continue; // empty line found, skip to next
		}
		if ( !(10 > charcnt) )
			return false; // line too long

		ret = sscanf(buf, "%lg",&magnitude);
		if (1 > ret){
			std::cerr << "Reading failed for line: " << buf << std::endl;
			return false;
		}
		_m.push_back(magnitude);
	}
#ifdef DEBUG
	std::cout << "Number of magnitudes in training data set: " << _m.size() << std::endl;
#endif
	_m_gsl = gsl_vector_alloc (_m.size());
	for (size_t i=0; i<_m.size();i++){
		gsl_vector_set(_m_gsl,i,_m[i]);
	}
	return true;
}

bool TData::read_distances(){
	int linecnt, charcnt, ret;
	double distance;
	char buf[10];

	for (linecnt = 0; !(*_isP).eof(); linecnt++ ) {
		(*_isP).getline(buf, 10);
		charcnt = (*_isP).gcount();
		if ( 1 > charcnt ) {
			continue; // empty line found, skip to next
		}
		if ( !(10 > charcnt) )
			return false; // line too long

		ret = sscanf(buf, "%lg",&distance);
		if (1 > ret){
			std::cerr << "Reading failed for line: " << buf << std::endl;
			return false;
		}
		_r.push_back(distance);
	}
#ifdef DEBUG
	std::cout << "Number of distances in training data set: " << _r.size() << std::endl;
#endif
	_r_gsl = gsl_vector_alloc (_r.size());
	for (size_t i=0; i<_r.size();i++){
		gsl_vector_set(_r_gsl,i,_r[i]);
	}
	return true;
}

bool TData::read(std::string dtype){
	bool success = false;
	switch(dtype[0]){
	case 't': success = read_traces(); break;
	case 'm': success = read_magnitudes(); break;
	case 'd': success = read_distances(); break;
	default : std::cerr << "Datatype has to be either 'traces', 'magnitudes', or 'distances'" << std::endl;
	}
	return success;
}

int TData::get_noftraces(){
	return _ntraces;
}

int TData::get_nofbands(){
	return _nbands;
}

void GbA::init(){
	datafiles df;
	datafiles::iterator dfit;
	_td = new TData;
	df["./data/az_training.txt"] = "traces";
	df["./data/m.txt"] = "magnitudes";
	df["./data/r.txt"] = "distances";
	for (dfit=df.begin(); dfit != df.end(); dfit++ )
		_td->load(dfit->first,dfit->second);
}

void GbA::compute_likelihood(double *data, int nstats, int nbands, int nsim,
						double mean[2], double cov[2][2]){
	int i, j, k, l;
	gsl_vector *lsq = gsl_vector_alloc (_td->get_noftraces());
	gsl_vector *input = gsl_vector_alloc (nbands);
	gsl_vector *trace = gsl_vector_alloc (nbands);
	double *mags, *dists;
	double misfit_l2;
	size_t *indices = new size_t[nsim];
	mags = new double[nsim];
	dists = new double[nsim];

	for (i=0; i<nstats; i++){
		for (j=0; j<nbands; j++){
			gsl_vector_set(input, j, std::log10(data[i*nbands+j]));
		}
		for (k=0; k<_td->get_noftraces(); k++){
			gsl_matrix_get_row(trace, _td->_tdata_gsl, k);
			gsl_vector_sub(trace,input);
			misfit_l2 = gsl_blas_dnrm2(trace);
			gsl_vector_set(lsq, k, misfit_l2*misfit_l2);
		}
		gsl_sort_vector_smallest_index(indices, nsim,lsq);
		mean[0] = 0;
		mean[1] = 0;
		for(l=0; l<nsim; l++){
			dists[l] = std::log10(gsl_vector_get(_td->_r_gsl, indices[l]));
			mags[l] = gsl_vector_get(_td->_m_gsl, indices[l]);
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
}
