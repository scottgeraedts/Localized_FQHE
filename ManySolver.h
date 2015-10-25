#ifndef MANYSOLVER_H
#define MANYSOLVER_H

#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>
#include "math.h"
#include <ctime>
#include "lapacke.h"
#include "matrixcontainer.h"
#include "SingleSolver.h"

//will use std::bitset library to store bits, which needs to know the number of bits at compile time
//making this a bit too big is probably fine, some required values
#define NBITS 30
using namespace std;

template <class ART>
class ManySolver{
public:
	ManySolver(int Ne,int charge,int periodic);
	void print_H();

	//could make the below virtual if I ever implement another geometry
	void make_Hnn();
	void ZeroHnn();
	void disorderHnn();
//	virtual void make_Hmm();
	
protected:
	int Ne,NPhi,nStates;
	int charge,has_charge,disorder,periodic,cache,project;
	vector<bitset<NBITS> > states;
	Eigen::Matrix<ART,Eigen::Dynamic, Eigen::Dynamic> Hnn;//ManySolver in Landau basis
	Eigen::Matrix<ART,Eigen::Dynamic, Eigen::Dynamic> Hmm;//ManySolver in localized basis
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> > es;
	
	virtual ART two_body(int a,int b) =0;//could make virtual
	virtual ART four_body(int a,int b,int c,int d) =0;//could make virtual
	virtual int get_charge(bitset<NBITS> state)=0;//could make virtual
	
	void make_cache();
	vector<ART> four_body_cache, two_body_cache;
	ART four_body_project(int,int,int,int);
	ART two_body_project(int,int);
	Eigen::Matrix<ART,Eigen::Dynamic, Eigen::Dynamic> basis;//used to test projection
	ART get_interaction(int,int,int,int);
	ART get_interaction(int,int);
	ART get_disorder(int,int);
	int four_array_map(int,int,int,int);
	
	//stuff to deal with disorder
	SingleSolver single; //for other geometries, this would need to be promoted to a single-electron parent class
//	virtual void add_disorder()=0;
	
	void init();
	void make_states();
	int adjust_sign(int a,int b,bitset<NBITS> state);
	int hasbit(int i,int a);
	int lookup_flipped(bitset<NBITS> i,int a,int b,int c,int d);
	int lookup_flipped(bitset<NBITS> i,int a,int b);
	int state_to_index(int state);
	double V_Coulomb(double qx,double qy);

	//math functions
	double ClebschGordan(int,int,int);
	long int comb(int,int);
	double factorial(int,int);
	double factorial(int);
	int intpow(int,int);

};

///**************DEFINITIONS HERE************///

template<class ART>
ManySolver<ART>::ManySolver(int tNe,int tcharge,int tperiodic):Ne(tNe),charge(tcharge),periodic(tperiodic){
	disorder=0;
	cache=0;
	project=0;
	if (charge==-1) has_charge=0;
	else has_charge=1;
}
template<class ART>
void ManySolver<ART>::init(){
	make_states();
//	for(int i=0;i<nStates;i++) cout<<states[i]<<" ";
//	cout<<Hnn<<endl;

	if(cache){
		four_body_cache.resize(pow(NPhi,4));
		two_body_cache.resize(pow(NPhi,2));
		make_cache();
	}
	
	ZeroHnn();	
	make_Hnn();
//	if(disorder) disorderHnn();	

//	clock_t t;
//	t=clock();
	es.compute(Hnn);
//	t=clock()-t;
//	cout<<"eigen: "<<((float)t)/CLOCKS_PER_SEC<<" seconds"<<endl;	

//***a call to lapack, this can't handle complex numbers, changes Hnn, and is slower than Eigen. Probably should never be used
//	double *matrix=new double[nStates*nStates];
//	matrix=Hnn.data();
//	double *w=new double[nStates];
//	t=clock();
//	LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','U',nStates,matrix,nStates,w);
//	t=clock()-t;
//	cout<<"lapack: "<<((float)t)/CLOCKS_PER_SEC<<" seconds"<<endl;	
//	cout<<w[0]/Ne<<endl;
//	delete [] w;
//	delete [] matrix;

//****A call to ARPACK++. The fastest of all methods, but I don't know how to get it to work for complex numbers right now
//	MatrixContainer<ART> mat(nStates,Hnn);
//	t=clock();
//    ARSymStdEig<ART, MatrixContainer<ART> >  dprob(mat.ncols(), 9, &mat, &MatrixContainer<ART>::MultMv,"LM");//someday put this part into matprod?
//    dprob.FindEigenvalues();
//	for(int i=0;i<dprob.ConvergedEigenvalues();i++) cout<<dprob.Eigenvalue(i)<<endl;

//	t=clock()-t;
//	cout<<"arpack: "<<((float)t)/CLOCKS_PER_SEC<<" seconds"<<endl;	
}


/// A state can be represented by a bit string, 
//ex. 7=11100000 has electrons in positions 0,1,2, 13=1011 has electron in 0,2,3
//this function starts at 0 and sees if its bitstring representation has the right number of electrons
//it it does, it adds it to a vector called states
template<class ART>
void ManySolver<ART>::make_states(){
	states=vector<bitset<NBITS> >(comb(NPhi,Ne),bitset<NBITS>());
	int j=0;
	bitset<NBITS> s;
	for(int i=0;i<intpow(2,NPhi);i++){
		s=bitset<NBITS>(i);
		if (s.count()==(unsigned) Ne && (has_charge==0 || get_charge(s)==charge) ){
			states[j]=s;
			j++;
		}
	}
	states.resize(j);
	nStates=j;		 
//	cout<<"nStates: "<<nStates<<endl;
}
template<class ART>
void ManySolver<ART>::make_Hnn(){
	int d,j;
	ART temp;
	for(signed int a=0;a<NPhi;a++){

		//terms from the disorder potential
		if(disorder){
			for(int b=0;b<NPhi;b++){
				temp=get_disorder(a,b);
				for(int i=0;i<(signed)nStates;i++){
					if (states[i].test(a) && (b==a || !states[i].test(b)) ){
						j=lookup_flipped(states[i],a,b);
						Hnn(i,j)+=(double)adjust_sign(a,b,states[i]) * temp;
					}
				}
			}
		}

		for(signed int b=a+1;b<NPhi;b++){
			if(b==a) continue;
			
			//this is for `two body' terms which are four-body interaction terms where the indices come in pairs
			temp=get_interaction(a,b);
			for(int i=0;i<(signed)nStates;i++){
				if(states[i].test(a) && states[i].test(b)){
					Hnn(i,i)+=temp;
				}
			}
			
			//a+b=c+d, this sets some restrictions on which c and d need to be summed over
			//on a sphere, there are already restrictions on the possible c, though on the torus more c's are possible (its difficult to a priori determine which ones)
			int c_upper,c_lower;
			if(periodic){
				c_upper=NPhi; c_lower=0;
			}else{
				if (a+b>=NPhi){
					c_upper=NPhi; c_lower=a+b-NPhi+1;
				}else{
					c_upper=a+b; c_lower=0;
				}
			}
			for(int c=c_lower;c<c_upper;c++){
				if (c==a || c==b) continue;				
				if ((a+b-c)<0) d=a+b-c+NPhi;
				else d=(a+b-c)%NPhi;
				if (d==a || d==b || d<=c) continue;
				temp=get_interaction(a,b,c,d);
//				cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<temp<<endl;
				for(int i=0;i<nStates;i++){
					if(states[i].test(a) && states[i].test(b) && !states[i].test(c) && !states[i].test(d)){
						j=lookup_flipped(states[i],a,b,c,d);
						Hnn(i,j)+=(double)(adjust_sign(a,b,states[i])*adjust_sign(c,d,states[i])) * temp;
					}
				}	
			}	
		}
	}
}

template<class ART>
void ManySolver<ART>::make_cache(){
	int d;
	ART temp;
	
	for(signed int a=0;a<NPhi;a++){
		//terms from the disorder potential
		if(disorder){
			for(signed int b=0;b<NPhi;b++){		
				two_body_cache[four_array_map(a,b,0,0)]=single.getH(a,b);
			}
		}
		for(signed int b=a+1;b<NPhi;b++){

			//this is for `two body' terms which are four-body interaction terms where the indices come in pairs
			four_body_cache[four_array_map(a,b,a,b)]=two_body(a,b);
			

			//a+b=c+d, this sets some restrictions on which c and d need to be summed over
			//on a sphere, there are already restrictions on the possible c, though on the torus more c's are possible (its difficult to a priori determine which ones)
			int c_upper,c_lower;
			if(periodic){
				c_upper=NPhi; c_lower=0;
			}else{
				if (a+b>=NPhi){
					c_upper=NPhi; c_lower=a+b-NPhi+1;
				}else{
					c_upper=a+b; c_lower=0;
				}
			}
			for(int c=c_lower;c<c_upper;c++){
				if (c==a || c==b) continue;				
				if ((a+b-c)<0) d=a+b-c+NPhi;
				else d=(a+b-c)%NPhi;
				if (d==a || d==b || d<=c) continue;
				four_body_cache[four_array_map(a,b,c,d)]=four_body(a,b,c,d);
			}	
		}
	}
}

//adds only the disorder terms to the Hamiltonian, no longer used
template<class ART>
void ManySolver<ART>::disorderHnn(){
	int j;
	if (!disorder){
		cout<<"Disordered Hamiltonian requested, but no disorder has been set"<<endl;
		exit(0);
	}
	for(signed int a=0;a<NPhi;a++){
		//terms from the disorder potential that act on the same site
		for(int i=0;i<(signed)nStates;i++){
			if (states[i].test(a)) Hnn(i,i)+=get_disorder(a,a);
		}
		
		//terms from the disorder potential that act on different sites
		for(signed int b=a+1;b<NPhi;b++){			
			for(int i=0;i<(signed)nStates;i++){
				if (states[i].test(a) && !states[i].test(b)){
					j=lookup_flipped(states[i],a,b);
					Hnn(i,j)+=get_disorder(a,b);
				}
			}
		}
	}
}
template<class ART>
void ManySolver<ART>::ZeroHnn(){ Hnn=Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>::Zero(nStates,nStates);}

//returns the interaction terms from a Hamiltonian
template<class ART>
ART ManySolver<ART>::get_interaction(int a,int b){
	if(!cache) return two_body(a,b);
	else if(!project) return four_body_cache[four_array_map(a,b,a,b)];
	else return four_body_project(a,b,a,b);
}
template<class ART>
ART ManySolver<ART>::get_interaction(int a,int b,int c, int d){
	if(!cache) return four_body(a,b,c,d);
	else if(!project) return four_body_cache[four_array_map(a,b,c,d)];
	else return four_body_project(a,b,c,d);
}

//returns the disorder terms from a Hamiltonian
template<class ART>
ART ManySolver<ART>::get_disorder(int a, int b){
	if(!cache) return single.getH(a,b);
	else if(!project) return two_body_cache[four_array_map(a,b,0,0)];
	else return two_body_project(a,b);
}

//maps indices from a four array to indices from a 1 array
template<class ART>
int ManySolver<ART>::four_array_map(int a,int b,int c, int d){ return a+b*NPhi+c*NPhi*NPhi+d*pow(NPhi,3); }

template<class ART>
ART ManySolver<ART>::four_body_project(int a,int b,int c,int d){
	ART out=0;
	int id,sa,sb,sc,sd;
	for(int ia=0;ia<NPhi;ia++){
		for(int ib=0;ib<NPhi;ib++){
			if (ib==ia) continue;

			if (ib<ia){ sa=ib;	sb=ia; //we only saved terms with a<b in the cache, so need to swap a and b if this isn't true
			}else{ sa=ia; sb=ib; }

			//pick bounds on c
			int c_upper,c_lower;
			if(periodic){
				c_upper=NPhi; c_lower=0;
			}else{
				if (a+b>=NPhi){
					c_upper=NPhi; c_lower=a+b-NPhi+1;
				}else{
					c_upper=a+b; c_lower=0;
				}
			}
			for(int ic=c_lower;ic<c_upper;ic++){
				if ((ia+ib-ic)<0) id=ia+ib-ic+NPhi;
				else id=(ia+ib-ic)%NPhi;				
				if (id==ic) continue;

				if (id<ic){ sc=id;	sd=ic;
				}else{ sc=ic; sd=id; }
					
				out+=basis(a,ia)*basis(b,ib)*basis(c,ic)*basis(d,id)*four_body_cache[four_array_map(sa,sb,sc,sd)];
			}
		}
	}
	return out;
}
template<class ART>
ART ManySolver<ART>::two_body_project(int a,int b){
	ART out=0;
	for(int ia=0;ia<NPhi;ia++){
		for(int ib=0;ib<NPhi;ib++){
			out+=basis(a,ia)*basis(b,ib)*two_body_cache[four_array_map(ia,ib,0,0)];
		}
	}
	return out;
}

//a Coulomb interaction
template<class ART>
double ManySolver<ART>::V_Coulomb(double qx,double qy){
	return 1./(sqrt(qx*qx+qy*qy)*NPhi);
}

//puts in minus signs for correct Fermi statistics
//for four-body terms, we have something like c^dagger_a c_d. This gives a (-1) for every electron between a and d
//this also works fine when you consider that we are performing two hops
template<class ART>
int ManySolver<ART>::adjust_sign(int a,int b,bitset<NBITS> state){
	int sign=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state.test(i)) sign*=-1;
	return sign;	
}
template<class ART>
int ManySolver<ART>::lookup_flipped(bitset<NBITS> state,int a,int b,int c,int d){
//given a bitstring, finds its index in the bitstring array
//this is currently done by a linear search, which is a terrible way to do it, and likely needs to be improved
	bitset<NBITS> compare=state.flip(a);
	compare=compare.flip(b);
	compare=compare.flip(c);
	compare=compare.flip(d);
	for (int i=0;i<nStates;i++){
		if (compare==states[i]) return i;
	}
	cout<<"error in lookup_flipped: "<<state<<" "<<compare<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
	exit(0);	
}

template<class ART>
int ManySolver<ART>::lookup_flipped(bitset<NBITS> state,int a,int b){
//given a bitstring, finds its index in the bitstring array
//this is currently done by a linear search, which is a terrible way to do it, and likely needs to be improved (the really sexy thing to do would be to combine all these functions)
	bitset<NBITS> compare=state.flip(a);
	compare=compare.flip(b);
	for (int i=0;i<nStates;i++){
		if (compare==states[i]) return i;
	}
	cout<<"error in lookup_flipped: "<<state<<" "<<compare<<" "<<a<<" "<<b<<endl;
	exit(0);	
}

template<class ART>
void ManySolver<ART>::print_H(){
	ofstream Hout;
	Hout.open("Hout");
	Hout<<Ne<<" "<<NPhi<<" "<<nStates<<endl;
	for(int i=0;i<nStates;i++){
		for(int j=0;j<nStates;j++){
			if ( abs(Hnn(i,j)) > 0.00000000001) Hout<<setprecision(12)<<i<<" "<<j<<"  "<<Hnn(i,j)<<endl;
		}
	}
	Hout.close();
}		

//////////////////////utility functions, generally	
//someday should probably make a math library that contains all these functions

//x choose p
template<class ART>
long int ManySolver<ART>::comb(int x,int p){
	return factorial(x,x-p)/factorial(p);	
}
//computes the products of numbers from start to stop (a partial factorial)
template<class ART>
double ManySolver<ART>::factorial(int start,int stop){
	if (stop==start) return 1;
	else return start*factorial(start-1,stop);
}	
template<class ART>
double ManySolver<ART>::factorial(int n){
	if (n==1 || n==0) return 1;    
	return n*factorial(n-1);
}	
//a function for computing powers of integers
template<class ART>
int ManySolver<ART>::intpow(int x, int p)
{
  if (p == 0) return 1;
  if (p == 1) return x;

  int tmp = intpow(x, p/2);
  if (p%2 == 0) return tmp * tmp;
  else return x * tmp * tmp;
}
template<class ART>
double ManySolver<ART>::ClebschGordan(int a,int b,int L){
//calculate CG coefficients from the formula on wikipedia
//m1=a-Q,m2=b-Q, 2Q=NPhi-1,these are stored this way because sometimes they are half-integer
//may eventually need to tabulate these since they are pretty slow to calculate
	double prefactor1=sqrt((2*L+1)*factorial(L)*factorial(L)*factorial((NPhi-1)-L)/(1.*factorial((NPhi-1)+L+1)) );
	double prefactor2=sqrt(factorial(L+a+b-(NPhi-1))*factorial(L-a-b+(NPhi-1))*factorial((NPhi-1)-a)*factorial((NPhi-1)-b)*factorial(a)*factorial(b));
	double sum=0.;
	int sign=1;
	for(int k=0;k<=b;k++){
		if (k%2==1) sign=-1;
		else sign=1;
		if(L-b+k<0) continue;
		if((NPhi-1)-L-k<0 || L-(NPhi-1)+a+k<0 || (NPhi-1)-a-k<0) continue;
		sum+=(sign*1.)/(1.*factorial((NPhi-1)-L-k)*factorial((NPhi-1)-a-k)*factorial(b-k)*factorial(L-(NPhi-1)+a+k)*factorial(L-b+k)*factorial(k));
	}
	return prefactor1*prefactor2*sum;
}

#endif
