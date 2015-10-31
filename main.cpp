#include "TorusSolver.h"

using namespace std;

//The main used for running the program
int main(int argc, char** argv) {

	int Ne,charge,NPhi;
	double V;
	ifstream params;
	params.open("params");
	params>>Ne;
	params>>NPhi;
	params>>charge;
	params>>V;
	TorusSolver<complex<double> > H(Ne,2,V,NPhi*Ne,"test");
//	for(int e=2;e<=3;e++){
		for(int n=4;n<=8;n++){ 
			for(int i=0;i<n;i++)
				H=TorusSolver<complex<double> >(n,i,0,2*n,"2");
		}
//	}

}
