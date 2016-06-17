
CC=icpc
HOME=/home/geraedts
LIBDIR=$(HOME)/Localized_FQHE/
MYDIR=$(HOME)/sources/lapack-3.5.0/lapacke/include/
MYCDIR=$(HOME)/myClibrary/
ARPACKPP_DIR = $(HOME)/sources/arpack++/include/
EIGEN_DIR= $(HOME)/sources/eigen-eigen-bdd17ee3b1b3
MKL_LINKLINE= -larpack -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
#MKL_LINKLINE= -larpack -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lpthread -lm
#MKL_COMPILE= -openmp 
MKL_COMPILE= -openmp -I$(MKLROOT)/include
CFLAGS=-O3 -Wall -I$(LIBDIR) -I$(ARPACKPP_DIR) -I$(MYCDIR) -I$(EIGEN_DIR) $(MKL_COMPILE)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
LIBS=-larpack $(HOME)/myClibrary/utils.o $(MKL_LINKLINE)
SOURCES= main.cpp SingleSolver.cpp Potential.cpp SphereSolver.cpp sphere.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=torus.out sphere.out

all: $(SOURCES) $(EXECUTABLE)
	
singletest.out: singletest.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) singletest.o SingleSolver.o Potential.o $(LIBS) -o singletest.out

torus.out: main.o SingleSolver.o Potential.o 
	$(CC) $(CFLAGS) main.o SingleSolver.o Potential.o  $(LIBS) -o torus.out

sphere.out: SphereSolver.o SingleSolver.o Potential.o
	$(CC) -I$(MYDIR) SphereSolver.o SingleSolver.o Potential.o $(LIBS) -o sphere.out

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
	
%.o: %.f90
	gfortran $(CFLAGS) -c $< -o $@
clean:
	rm *.o
