#ifndef SINGLE_SOLVER
#define SINGLE_SOLVER

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <math.h>
#include <fstream>
#include "MersenneTwister.h"
#include <map>

using namespace std;
using Eigen::MatrixXcd;

const int DT=1; //fineness of grid used when threading boudary conditions
const int CUTOFF=2; //max number of times you can hop around the torus in the y direction

class Potential{
public:
	Potential(string _type, int _qbounds, int seed, double _Lx, double _Ly, map <string, double> &params);
	complex<double> get_potential(int mx,int my);
	void plot_potential();
	int get_qbounds();
	vector<double> xloc,yloc,sign;//locations of delta function impurities
	
private:
	int qbounds;
	string type;
	double Lx,Ly;
	Eigen::MatrixXcd V;
	MTRand ran;
	void make_potential_delta(map <string, double> &params);
	void make_potential_lattice(map <string, double> &params);
	void make_potential_gaussian(map <string, double> &params);
	void make_potential_test(map <string, double> &params);
};

class SingleSolver{
public:
	SingleSolver(int NPhi,int N_deltas, double Lx=0, double Ly=0, int NROD=1, int shift=0, double Vstrength=1. );
	void testDistance();
	void visualizer();
	void run();
	complex<double> getH(int,int);
	void switchVstrength();
	complex<double> evec(int,int);
	void init(int, double, int, int);
	void printEnergy(int,int);

	SingleSolver () {} 
	Eigen::SelfAdjointEigenSolver<MatrixXcd> proj;	
//	SingleSolver& operator=(const SingleSolver& other);
		

private:
	int NPhi,qbounds,N_trunc,energy_mesh,NROD,shift;
	double Vstrength,Lx,Ly;
	Eigen::MatrixXcd dVdX,dVdY,Hnn;
	Eigen::VectorXd energies,spacings,energy_hist;

	void make_Hamiltonian(Potential &pot, double theta_x=0., double theta_y=0.);
	void plot_density(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, double theta_x=0., double theta_y=0.);
	void density_near_x(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, double xloc, double yloc, double theta_x=0., double theta_y=0.);
	Eigen::MatrixXcd Hnn_to_Hmm(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, int nLow);
	Eigen::VectorXd conductivity(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es_m);
	Eigen::VectorXd ThoulessNumber(Potential &deltas, Potential &gaussian, Eigen::VectorXd &oldEnergies);
	void print(Eigen::VectorXd &full_sigma, Eigen::VectorXd &full_thouless);
	
};

#endif
