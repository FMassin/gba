all: interface library gut

CC=g++
#CCFLAGS= -g -Wall -DDEBUG
CCFLAGS=-g -Wall
LIBS=-lgsl -lcblas -latlas -lm -lnetcdf_c++ 

ofiles= gba.o

%.o:%.cpp
	$(CC) -c $(CCFLAGS) $<
	 
interface:
	swig -python -c++ GbA.i

library:
	g++ -fPIC $(CCFLAGS) -c  gba.cpp ;\
	g++ -fPIC  $(CCFLAGS) -c GbA_wrap.cxx -I/usr/include/python2.7;\
 	g++ -shared $(FLAGS) gba.o GbA_wrap.o -o _gba.so $(LIBS)

gut: $(ofiles)
	$(CC) $(CCFLAGS) -o gba main.cpp $(ofiles) $(LIBS)
	
clean:
	-rm *.o *.so *.pyc *.cxx *~
	