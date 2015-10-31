#load makefile from arpack
#include /home/sgeraedt/sources/arpack++/Makefile.inc

CC=g++
LIBDIR=/home/sgeraedt/Localized_FQHE/
MYDIR=/home/sgeraedt/sources/lapack-3.5.0/lapacke/include/
ARPACKPP_DIR = /home/sgeraedt/sources/arpack++/include/
CFLAGS=-c -O3 -Wall -I$(LIBDIR) -I$(MYDIR) -I$(ARPACKPP_DIR)
#LDFLAGS=-I$(LIBDIR) -I$(MYDIR) 
SOURCES= main.cpp SingleSolver.cpp Potential.cpp
OBJECTS=$(SOURCES:.cpp=.o) 
EXECUTABLE=a.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) -I$(MYDIR) $(OBJECTS) -larpack -llapacke -llapack -lcblas -lrefblas -lgfortran -o $@ 
	
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
clean:
	rm *.o
