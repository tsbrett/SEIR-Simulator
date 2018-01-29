CC=g++
LIBS=-lgsl -lgslcblas -std=c++11



seir_simulator: seir_simulator.cpp seir_parameters.h
	@echo Compiling $<...
	$(CC) -o $@ $< $(LIBS)


		
