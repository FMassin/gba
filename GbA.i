/* File : gba.i */
%module gba
%include "typemaps.i"
%include "std_string.i"


%{
#define SWIG_FILE_WITH_INIT
#include "gba.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* data, int nbands)};
%apply (double INPLACE_ARRAY1[ANY]) {double mean[2]};
%apply (double INPLACE_ARRAY2[ANY][ANY]) {double cov[2][2]};

/* Let's just grab the original header file here */
%include "gba.h"

