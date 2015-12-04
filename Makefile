#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=g++
HOME=/home/sgeraedt
LIBDIR=$(HOME)/Localized_FQHE/
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
MYCDIR=$(HOME)/myClibrary/
ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
CFLAGS=-c -O3 -Wall -I$(LIBDIR) -I$(MYDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
LIBS=-larpack -lcblas -lrefblas -lgfortran $(HOME)/myClibrary/utils.o
SOURCES= main.cpp SingleSolver.cpp Potential.cpp SphereSolver.cpp sphere.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=torus.out sphere.out

all: $(SOURCES) $(EXECUTABLE)
	
torus.out: main.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) main.o SingleSolver.o Potential.o $(LIBS) -o torus.out

sphere.out: SphereSolver.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) SphereSolver.o SingleSolver.o Potential.o $(LIBS) -o sphere.out

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
