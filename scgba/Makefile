LIBS=-lboost_filesystem -lboost_iostreams -lboost_program_options -lboost_regex -lboost_signals -lboost_system -lboost_thread -lbson-1.0

scgba: main.cpp ../gba.cpp
	g++ -g -o scgba main.cpp ../gba.o -I../ -I/home/behry/eewamps/include \
	-L/usr/local/src/eewamps/build/lib -lseiscomp3_client -lseiscomp3_core \
	-lseiscomp3_eewamps $(LIBS) -lgsl -lcblas -latlas -lm -lnetcdf_c++ 