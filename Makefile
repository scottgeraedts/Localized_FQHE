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
LIBS=-larpack -llapack -lcblas -lrefblas -lgfortran $(HOME)/myClibrary/utils.o
SOURCES= main.cpp SingleSolver.cpp Potential.cpp SphereSolver.cpp sphere.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=torus.out sphere.out

all: $(SOURCES) $(EXECUTABLE)
	
singletest.out: singletest.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) singletest.o SingleSolver.o Potential.o $(LIBS) -o singletest.out

torus.out: main.o SingleSolver.o Potential.o new_coulomb_m.o wf_tools.o z_function_m.o
	$(CC) -I$(MYDIR) main.o SingleSolver.o Potential.o new_coulomb_m.o wf_tools.o z_function_m.o $(LIBS) -o torus.out

sphere.out: SphereSolver.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) SphereSolver.o SingleSolver.o Potential.o $(LIBS) -o sphere.out

twobody.out: twobody_test.o wf_tools.o z_function_m.o new_coulomb_m.o 
	$(CC) twobody_test.o wf_tools.o z_function_m.o new_coulomb_m.o $(LIBS) -o twobody.out
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
%.o: %.f90
	gfortran $(CFLAGS) $< -o $@
clean:
	rm *.o
