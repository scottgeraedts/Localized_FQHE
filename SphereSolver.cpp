#include "SphereSolver.h"
SphereSolver::SphereSolver(int x):ManySolver<double>(){
	//stuff unique to the sphere
//	NPhi=3*Ne-tshift+1;

	periodic=0;
	make_states();
	if (NPhi%2==1) Qeven=1;
	else Qeven=0;

	make_VL_haldane();
	cout<<nStates<<" "<<Ne<<" "<<NPhi<<" "<<charge<<" "<<has_charge<<endl;
	make_Hnn();
	EigenDenseEigs();
	for(int i=0;i<nStates;i++) 	
		cout<<setprecision(7)<<"energies: "<<eigvals[i]<<endl;
//		cout<<setprecision(7)<<"energies: "<<(eigvals[i]-pow(Ne,2)/sqrt((NPhi-1)*2))/(1.*Ne)<<endl;
//	for(int j=0;j<nStates;j++){
//		for(int i=0;i<nStates;i++) cout<<eigvecs[j][i]<<" "<<(bitset<10>)states[i]<<endl;
//		cout<<endl;
//	}
}

//compute formula 3.224 from jain. 2Q=NPhi-1, R=sqrt(Q). NPhite only need elements with (L-2Q)%2==1
//the factorials have all been reorganized to try to prevent integer overflow
void SphereSolver::make_VL_coulomb(){
	VL=vector<double>(NPhi,0);
	double part1,part2,part3,part4;
	for (int L=NPhi%2;L<NPhi;L+=2){
//		VL[L]=(2.*comb(2*(NPhi-1)-2*L,NPhi-1-L)*comb(2*(NPhi-1)+2*L+2,NPhi-1+L+1))/(1.*intpow(comb(2*(NPhi-1)+2,(NPhi-1)+1),2)*sqrt( (1.0*(NPhi-1) )/2. )) ;
		part1=factorial(2*(NPhi-1)+2*L+2,2*(NPhi-1)+2);
		part2=factorial((NPhi-1)+1,(NPhi-1)-L);
		part3=factorial(2*(NPhi-1)+2,2*(NPhi-1)-2*L);
		part4=factorial((NPhi-1)+L+1,(NPhi-1)+1);
		VL[L]=(2.*part1*pow(part2,2))/(sqrt( (1.0*(NPhi-1))/2.)*part3*pow(part4,2));
	}
}
void SphereSolver::make_VL_haldane(){
	VL=vector<double>(NPhi,0);
	VL[NPhi-1-1]=1.;
}	
double SphereSolver::two_body(int a,int b){
	cout<<"two body was called!"<<endl;
	double out=0;
	for(int L=abs(a+b-(NPhi-1));L<=(NPhi-1);L++){
		if( (L-(NPhi-1))%2==0) continue;
//		cout<<"***"<<a<<" "<<b<<" "<<L<<" "<<ClebschGordan(a,b,L)<<endl;
		out+=2.*VL[L]*pow(ClebschGordan(a,b,L,NPhi),2);
	}	
	return out;
}
double SphereSolver::four_body(int a,int b,int c,int d){
	if (d != a+b-c) return 0.; //hack so that I can use makeHnn for both periodic and finite bc
	double out=0;
	for(int L=abs(a+b-(NPhi-1));L<=(NPhi-1);L++){
		if( (L-(NPhi-1))%2==1) continue;
		out+=2*VL[L]*ClebschGordan(a,b,L,NPhi)*ClebschGordan(d,c,L,NPhi);
	}		
	return out;
}

int SphereSolver::get_charge(int b){
	int out=0;
	for (int i=0;i<NPhi;i++)
		if(b & 1<<i)  out+=2*i-NPhi+1;	
	return out/2;
}
int main(int argc, char** argv) {

	SphereSolver H(1);
}
