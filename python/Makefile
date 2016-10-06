all: interface library

CC=g++
CCFLAGS=-g -Wall -DDEBUG
LIBS=-lgsl -lcblas -latlas -lm -lnetcdf_c++4 
NUMPYDEV=/usr/local/lib/python2.7/dist-packages/numpy/core/include/
PYDEV=/usr/include/python2.7

interface:
	swig -python -c++ GbA.i

library: ../gba.cpp ../gba.h
	g++ -fPIC $(CCFLAGS) -c $(abspath ../gba.cpp) -I../;\
	g++ -fPIC  $(CCFLAGS) -c GbA_wrap.cxx -I${PYDEV} -I${NUMPYDEV} -I../ ;\
 	g++ -shared $(FLAGS) gba.o GbA_wrap.o -o _gba.so $(LIBS)

clean:
	-rm *.o *.so *.pyc *.cxx *~
	