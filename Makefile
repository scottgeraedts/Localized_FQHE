#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=g++
HOME=/home/geraedts
LIBDIR=$(HOME)/Localized_FQHE/
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
MYCDIR=$(HOME)/myClibrary/
ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
EIGEN_DIR= $(HOME)/sources/eigen-eigen-bdd17ee3b1b3
CFLAGS=-c -O3 -Wall -I$(LIBDIR) -I$(MYDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR) -I$(EIGEN_DIR)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
LIBS=-larpack $(HOME)/sources/CBLAS/cblas_feynman.a -lblas -lgfortran $(HOME)/myClibrary/utils.o
SOURCES= main.cpp SingleSolver.cpp Potential.cpp SphereSolver.cpp sphere.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=torus.out sphere.out

all: $(SOURCES) $(EXECUTABLE)
	
singletest.out: singletest.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) singletest.o SingleSolver.o Potential.o $(LIBS) -o singletest.out

torus.out: main.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) main.o SingleSolver.o Potential.o $(LIBS) -o torus.out

sphere.out: SphereSolver.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) SphereSolver.o SingleSolver.o Potential.o $(LIBS) -o sphere.out

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
