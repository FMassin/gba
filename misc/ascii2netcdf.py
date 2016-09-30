#!/usr/bin/env python

import numpy as np
from Scientific.IO.NetCDF import NetCDFFile

tdata = np.loadtxt('az_training.txt',delimiter=',')
mags = np.loadtxt('m.txt')
dists = np.loadtxt('r.txt')

fout = 'az_training.nc'

nid = NetCDFFile(fout,'w')
nid.createDimension('time',1);
nid.createDimension('traces',tdata.shape[0]);
nid.createDimension('filter',tdata.shape[1]);

t = nid.createVariable('time',np.dtype(float).char,('time',))
t.units = 's'
f = nid.createVariable('filter',np.dtype(float).char,('filter',))
f.units = 'Hz'
m = nid.createVariable('magnitude',np.dtype(float).char,('traces',))
m[:] = mags
ed = nid.createVariable('epicdist',np.dtype(float).char,('traces',))
ed.units = 'km'
ed[:] = dists
z = nid.createVariable('z',np.dtype(float).char,('time','traces','filter'))
z[0,:,:] = np.log10(tdata)
h = nid.createVariable('h',np.dtype(float).char,('time','traces','filter'))
h[0,:,:] = np.log10(tdata)
nid.close()


