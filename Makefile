all: interface library cppexample

CC=g++
CCFLAGS=-g -Wall -DDEBUG
LIBS=-lgsl -lcblas -latlas -lm -lnetcdf_c++4 
NUMPYDEV=/usr/local/lib/python2.7/dist-packages/numpy/core/include/
PYDEV=/usr/include/python2.7

ofiles= gba.o

%.o:%.cpp
	$(CC) -c $(CCFLAGS) $<
	 
interface:
	swig -python -c++ GbA.i

library:
	g++ -fPIC $(CCFLAGS) -c gba.cpp ;\
	g++ -fPIC  $(CCFLAGS) -c GbA_wrap.cxx -I${PYDEV} -I${NUMPYDEV} ;\
 	g++ -shared $(FLAGS) gba.o GbA_wrap.o -o _gba.so $(LIBS)

cppexample: $(ofiles)
	$(CC) $(CCFLAGS) -o gba_example example.cpp $(ofiles) $(LIBS)
	
clean:
	-rm *.o *.so *.pyc *.cxx *~
	