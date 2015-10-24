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
	TorusSolver<complex<double> > H(Ne,charge,0.1);
//	for(int i=1;i<20;i++) H=TorusSolver<complex<double> >(Ne,charge,i*0.1);
//	TorusSolver<double> H(Ne,charge);

}

