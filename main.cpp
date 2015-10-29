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
	TorusSolver<complex<double> > H(Ne,charge,0.1,3*Ne);
	//for(int i=1;i<Ne;i++) H=TorusSolver<complex<double> >(Ne,i,0.1,3*Ne);
//	TorusSolver<double> H(Ne,charge);

}

