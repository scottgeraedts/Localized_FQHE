#ifndef SPHERESOLVER_H
#define SPHERESOLVER_H

#include "ManySolver.h"
class SphereSolver:ManySolver<double>{
public:
	SphereSolver(int);
	void energy_spectrum();
	void ground_state();
	
private:
	int dQ;
	vector<double> VL;

	void make_VL_coulomb();
	void make_VL_haldane();	
	void make_VL_Tpfaff();	
//	double ClebschGordan(int m1,int m2,int L);
	double two_body(int a,int b);
	double four_body(int a,int b,int c,int d);
	double six_body(int a, int b, int c, int d, int e, int f);
	int get_charge(int b);
//	int adjust_sign(int a,int b,bitset<NBITS> state);
//	int hasbit(int i,int a);
//	int lookup_flipped(bitset<NBITS> i,int a,int b,int c,int d);
//	int state_to_index(int state);
//	int getLz(bitset<NBITS>);
//	long int comb(int,int);
//	double factorial(int,int);
//	double factorial(int);
//	int intpow(int,int);

};



#endif
