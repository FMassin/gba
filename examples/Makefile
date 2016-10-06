CC=g++
CCFLAGS=-g -Wall -DDEBUG
LIBS=-lgsl -lcblas -latlas -lm -lnetcdf_c++4 

cppexample: ../gba.cpp ../gba.h
	$(CC) $(CCFLAGS) -o gba_example example.cpp $(abspath ../gba.cpp) -I../ $(LIBS)