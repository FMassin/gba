all: library

CC=g++
CCFLAGS=-g -Wall -DDEBUG
LIBS=-lgsl -lcblas -latlas -lm -lnetcdf_c++4 
	 
library:
	g++ -fPIC $(CCFLAGS) -c $(abspath gba.cpp) ;\
 	g++ -shared $(FLAGS) gba.o -o gba.so $(LIBS)

clean:
	-rm *.o *.so
	