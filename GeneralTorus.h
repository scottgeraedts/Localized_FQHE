#ifndef GENERAL_TORUS_H
#define GENERAL_TORUS_H

#include "ManySolver.h"

class GeneralTorus:public ManySolver<complex<double> >{

public:
	GeneralTorus(int tcharge);
	GeneralTorus(int tcharge, double alpha, double theta);
	void ground_state(int charge);
	Eigen::SparseMatrix< complex<double> > shrinkMatrix;
	Eigen::SparseMatrix<complex<double> > density_operator(int mx, int my);
	void makeShrinker(int nx);
	
private:
	void make_Vmk();
	void interaction_cache();
	int get_charge(int);

	
	complex<double> two_body(int a, int b);
	complex<double> four_body(int a, int b, int c, int d);
	
	double Lx,Ly,LDelta;
	complex<double> tau;
	vector< vector< complex<double> > > Vmk;	

};


#endif
