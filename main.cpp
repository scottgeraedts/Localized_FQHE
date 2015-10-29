#include "TorusSolver.h"

using namespace std;

//The main used for running the program
int main(int argc, char** argv) {

	int Ne,charge;
	double V;
	ifstream params;
	params.open("params");
	params>>Ne;
	params>>charge;
	params>>V;
	TorusSolver<complex<double> > H(3,charge,0.1,3*3);
	for(int i=4;i<7;i++) H=TorusSolver<complex<double> >(i,charge,i*0.1,3*i);
//	TorusSolver<double> H(Ne,charge);

}

