#include "SphereSolver.h"
SphereSolver::SphereSolver(int x):ManySolver<double>(){
	//stuff unique to the sphere
//	NPhi=3*Ne-tshift+1;

	periodic=0;
	dQ=NPhi-1; //2*Q, useful in a lot of formulae
//	make_VL_Tpfaff();
	make_VL_haldane();
//	make_VL_coulomb();

	ground_state();		
}
void SphereSolver::energy_spectrum(){
	for(int c=Ne%2;c<=Ne*(Ne+1-Ne%2);c+=2){
		charge=c;
		make_states();

		make_Hnn();
//		make_Hnn_six();
//		ph_symmetrize();

		EigenDenseEigs();
		for(int i=0;i<nStates;i++) 	
			cout<<setprecision(7)<<charge<<" "<<eigvals[i]<<endl;

	}
}

void SphereSolver::ground_state(){
	make_states();

	make_Hnn();
		make_Hnn_six();
//		ph_symmetrize();
	EigenDenseEigs();
//	eigenvalues(10,-100);
	for(int i=0;i<eigvals.size();i++) 	
		cout<<setprecision(7)<<"energy: "<<eigvals[i]<<endl;

		Eigen::VectorXd eigen_eigvec(eigvecs[0].size());
		for(int i=0;i<eigvecs[0].size();i++) eigen_eigvec(i)=eigvecs[0][i];
		cout<<"eigenstate precision: "<<calcVarEigen(eigen_eigvec)<<endl;
		cout<<"eigenstate: "<<endl;
		double factor,smallest=1000;
		vector<double> factored_vec(nStates);
		for(int i=0;i<nStates;i++){
			factor=1;
			for(int s=0;s<NPhi;s++) 
				if(states[i] & 1<<s) factor*=sqrt( comb(NPhi-1,s) );	
			factored_vec[i]=eigvecs[0][i]*factor;
			cout<<eigvecs[0][i]<<" "<<(bitset<24>)states[i]<<endl;
			if(abs(factored_vec[i])<smallest && abs(factored_vec[i])>1e-6) smallest=abs(factored_vec[i]);
		}
		cout<<"factored eigenstate"<<endl;
		for(int i=0;i<nStates;i++){
			cout<<factored_vec[i]/smallest<<"\t"<<(bitset<24>)states[i]<<endl;
		}

///	plot_spectrum("spect");

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
	VL[NPhi-1-1]=0.;
}
void SphereSolver::make_VL_Tpfaff(){
// \sum_M C^{JQ}_{m_1+m_2,M}**2, assuming m1+m2=0 (though it should be the same for any m1+m2
	VL=vector<double>(NPhi,0);
	double V1=0, V3=0,temp1=0,temp3=0,temp1_2=0,temp3_2=0;	
	for(int dJ=2*dQ-6;dJ<2*dQ;dJ+=4){
		temp1+=sqrt( (dJ+1)*(2*dQ-2+1))*Wigner6j(dQ,dQ,dJ,dQ,3*dQ-6,2*dQ-2);			
		temp3+=sqrt( (dJ+1)*(2*dQ-6+1))*Wigner6j(dQ,dQ,dJ,dQ,3*dQ-6,2*dQ-6);			
		for(int dJ_2=2*dQ-6;dJ_2<2*dQ;dJ_2+=4){
			temp1_2+=4*sqrt( (dJ+1)*(dJ_2+1))*(2*dQ-2+1)*Wigner6j(dQ,dQ,dJ,dQ,3*dQ-6,2*dQ-2)*Wigner6j(dQ,dQ,dJ_2,dQ,3*dQ-6,2*dQ-2);
			temp3_2+=4*sqrt( (dJ+1)*(dJ_2+1))*(2*dQ-6+1)*Wigner6j(dQ,dQ,dJ,dQ,3*dQ-6,2*dQ-6)*Wigner6j(dQ,dQ,dJ_2,dQ,3*dQ-6,2*dQ-6);
		}
	}
	for(int M=-dQ;M<=dQ;M+=2){
		if(abs(M) > 3*dQ-6) continue;
		V1+=pow(ClebschGordan(2*dQ-2,dQ,0,M,3*dQ-6),2)*(1+4*temp1*lil_sign(dQ+1) + temp1_2);
		V3+=pow(ClebschGordan(2*dQ-6,dQ,0,M,3*dQ-6),2)*(1+4*temp3*lil_sign(dQ+1) + temp3_2);
	}
//	cout<<setprecision(12)<<"calculated pseudopotentials: V1="<<V1<<" V3="<<V3<<endl;
	VL[dQ-1]=V1;
	VL[dQ-3]=V3;
}	
double SphereSolver::two_body(int a,int b){
	cout<<"two body was called!"<<endl;
	double out=0;
	for(int L=abs(a+b-(NPhi-1));L<=(NPhi-1);L++){
		if( (L-(NPhi-1))%2==0) continue;
//		cout<<"***"<<a<<" "<<b<<" "<<L<<" "<<ClebschGordan(a,b,L)<<endl;
	//	out+=2.*VL[L]*pow(ClebschGordan(a,b,L,NPhi),2);
	}	
	return out;
}
double SphereSolver::four_body(int a,int b,int c,int d){
	if (d != a+b-c) return 0.; //hack so that I can use makeHnn for both periodic and finite bc
	double out=0;
	for(int L=abs(a+b-(NPhi-1));L<=(NPhi-1);L++){
		if( (L-(NPhi-1))%2==1) continue;
		out+=2*VL[L]*ClebschGordan(NPhi-1, NPhi-1, 2*a-(NPhi-1),2*b-(NPhi-1), 2*L)*ClebschGordan(NPhi-1, NPhi-1, 2*d-(NPhi-1), 2*c-(NPhi-1) , 2*L);
	}		
	return out;
}
double SphereSolver::six_body(int a, int b, int c, int d, int e, int f){
	double outfirst=0, outsecond=0;

	vector<int> invals(6,0);
	invals[0]=a; invals[1]=b; invals[2]=c; invals[3]=d; invals[4]=e; invals[5]=f;
	for(int i=0;i<6;i++) invals[i]=2*invals[i]-(NPhi-1); //convert an m from 0 to 2Q+1 to one that goes from -Q to Q

	double temp;
//	for(int dJ1=2*NPhi-8;dJ1<2*(NPhi-1);dJ1+=4){
//		for(int x1=0;x1<3;x1++){
//			if(dJ1>=abs(invals[x1]+invals[(x1+1)%3]) ){
//				temp=ClebschGordan(NPhi-1, NPhi-1, invals[x1],invals[(x1+1)%3], dJ1)*
//					ClebschGordan(dJ1, dQ, invals[x1]+invals[(x1+1)%3], invals[(x1+2)%3], 3*dQ-6 );
////				cout<<"first: "<<dJ1<<" "<<x1<<" "<<temp<<endl;
//				outfirst+=temp;
//			}
//			if(dJ1>=abs(invals[x1+3]+invals[(x1+1)%3+3]) ){
//				temp=ClebschGordan(NPhi-1, NPhi-1, invals[x1+3],invals[(x1+1)%3+3], dJ1)*
//					ClebschGordan(dJ1, dQ, invals[x1+3]+invals[(x1+1)%3+3], invals[(x1+2)%3+3], 3*dQ-6 );
////				cout<<"second: "<<dJ1<<" "<<x1<<" "<<temp<<endl;
//				outsecond+=temp;	
//			}
//		}
//	}
//	cout<<"testing six body basis conversion"<<endl;
//	int x1=0; temp=0;
//	for(int dJ1=2*dQ-6;dJ1<=2*dQ-2;dJ1+=4){
//		//first term
//		if(dJ1>=abs(invals[x1]+invals[(x1+1)%3]) ){
//			temp=ClebschGordan(NPhi-1, NPhi-1, invals[0],invals[1], dJ1)*
//				ClebschGordan(dJ1, dQ, invals[0]+invals[1], invals[2], 3*dQ-6 );	
////			cout<<"first: "<<dJ1<<" "<<0<<" "<<temp<<endl;
//			outfirst+=temp;
//		}
//		temp=0;
//		for(int dJ12=2*dQ-6;dJ12<=2*dQ;dJ12+=2){
//			if(dJ12>=abs(invals[0]+invals[1])){
//				temp+=sign*lil_sign((dJ1+dJ12)/2)*ClebschGordan(NPhi-1, NPhi-1, invals[0],invals[1], dJ12)*
//					ClebschGordan(dJ12, dQ, invals[0]+invals[1], invals[2], 3*dQ-6 )*
//					sqrt( (dJ12+1)*(dJ1+1) )*Wigner6j(dQ,dQ,dJ1,dQ,3*dQ-6,dJ12);
//			}	
//		}
////		cout<<"first: "<<dJ1<<" "<<2<<" "<<temp<<endl;
//		outfirst+=temp;
//		temp=0;
//		for(int dJ12=2*dQ-6;dJ12<=2*dQ;dJ12+=2){
//			if(dJ12>=abs(invals[0]+invals[1]) ){
//				temp+=sign*ClebschGordan(NPhi-1, NPhi-1, invals[0],invals[1], dJ12)*
//					ClebschGordan(dJ12, dQ, invals[0]+invals[1], invals[2], 3*dQ-6 )*
//					sqrt( (dJ12+1)*(dJ1+1) )*Wigner6j(dQ,dQ,dJ1,dQ,3*dQ-6,dJ12);
//			}
//		}	
////		cout<<"first: "<<dJ1<<" "<<1<<" "<<temp<<endl;
//		outfirst+=temp;	
//		temp=0;
//		//second term
//		if(  dJ1>=abs(invals[3]+invals[4]) ){
//			temp=ClebschGordan(NPhi-1, NPhi-1, invals[3],invals[4], dJ1)*
//				ClebschGordan(dJ1, dQ, invals[3]+invals[4], invals[5], 3*dQ-6 );		
////			cout<<"second: "<<dJ1<<" "<<0<<" "<<temp<<endl;
//			outsecond+=temp;
//		}
//		temp=0;
//		for(int dJ12=2*dQ-6;dJ12<=2*dQ;dJ12+=2){
//			if(dJ12>=abs(invals[4]+invals[3]) ){
//				temp+=sign*lil_sign((dJ1+dJ12)/2)*ClebschGordan(NPhi-1, NPhi-1, invals[3],invals[4], dJ12)*
//					ClebschGordan(dJ12, dQ, invals[4]+invals[3], invals[5], 3*dQ-6 )*
//					sqrt( (dJ12+1)*(dJ1+1) )*Wigner6j(dQ,dQ,dJ1,dQ,3*dQ-6,dJ12);
//			}
//		}
////		cout<<"second: "<<dJ1<<" "<<2<<" "<<temp<<endl;
//		outsecond+=temp;
//		temp=0;
//		for(int dJ12=2*dQ-6;dJ12<=2*dQ;dJ12+=2){			
//			if(dJ12>=abs(invals[4]+invals[3])){
//			temp+=sign*ClebschGordan(NPhi-1, NPhi-1, invals[3],invals[4], dJ12)*
//				ClebschGordan(dJ12, dQ, invals[4]+invals[3], invals[5], 3*dQ-6 )*
//				sqrt( (dJ12+1)*(dJ1+1) )*Wigner6j(dQ,dQ,dJ1,dQ,3*dQ-6,dJ12);
//			}
//		}	
////		cout<<"second: "<<dJ1<<" "<<1<<" "<<temp<<endl;
//		outsecond+=temp;
//	}
//	
//	cout<<"another test of six body"<<endl;
	outfirst=0; outsecond=0;
	double temp2;
	for(int dJ12=2*dQ-6; dJ12<2*dQ; dJ12+=4){
		temp=0;
		for(int dJ=2*dQ-6;dJ<2*dQ;dJ+=4)
			temp+=sqrt( (dJ+1)*(dJ12+1) ) * Wigner6j(dQ,dQ,dJ,dQ,3*dQ-6,dJ12);
		if(dJ12>=abs(invals[0]+invals[1])){
			temp2=ClebschGordan(dQ,dQ,invals[0],invals[1],dJ12)*ClebschGordan(dJ12,dQ,invals[0]+invals[1],invals[2],3*dQ-6)*(1+2*lil_sign(dQ+1)*temp);
			outfirst+=temp2;
	//		cout<<dJ12<<" "<<temp2<<endl;
		}
		if(dJ12>=abs(invals[3]+invals[4]) ){
			temp2=ClebschGordan(dQ,dQ,invals[3],invals[4],dJ12)*ClebschGordan(dJ12,dQ,invals[3]+invals[4],invals[5],3*dQ-6)*(1+2*lil_sign(dQ+1)*temp);
			outsecond+=temp2;
	//		cout<<dJ12<<" "<<temp2<<endl;
		}
	}	
	return outfirst*outsecond;	
}
//actually returns twice the charge to prevent me having to deal with fractions
int SphereSolver::get_charge(int b){
	int out=0;
	for (int i=0;i<NPhi;i++)
		if(b & 1<<i)  out+=2*i-NPhi+1;	
	return out;
}
int main(int argc, char** argv) {

	SphereSolver H(1);
}
