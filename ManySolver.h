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
#include "SingleSolver.h"
#include "wavefunction.h"
#include "utils.h"
#include <string>
#include "matprod.h"
extern "C"{
#include "cblas.h"
}
//#include "arscomp.h"

//will use std::bitset library to store bits, which needs to know the number of bits at compile time
//making this a bit too big is probably fine, some required values
using namespace std;

template <class ART>
class ManySolver:public MatrixWithProduct<ART>, public Wavefunction<ART>{
public:
	ManySolver();
	void print_H();

	void ZeroHnn();
	void make_Hnn();
	void make_Hnn_six();
	void disorderHnn();
	void make_states();
	
protected:
	int Ne,NPhi,nStates,NPhi2,NPhi3;
	int nHigh,nLow,NROD,random_offset;
	int charge,has_charge,disorder,periodic,cache,project,lookups, disordered_projection;
	double disorder_strength,V3overV1;
	vector<double> HaldaneV;
	string outfilename;

	vector<int > states;
	
	virtual ART two_body(int a,int b) =0;//could make virtual
	virtual ART four_body(int a,int b,int c,int d) =0;//could make virtual
	virtual ART six_body(int a, int b, int c);
	virtual int get_charge(int state)=0;//could make virtual

	void init();
	void interaction_cache();
	void disorder_cache();
	void projected_interaction_cache();

	vector<ART> four_body_cache, two_body_cache;
	vector<ART> projected_four_body_cache;
	ART four_body_project(int,int,int,int);
	ART two_body_project(int,int);
	ART get_interaction(int,int,int,int);
	ART get_disorder(int,int);
	int four_array_map(int,int,int,int);
	
	//stuff to deal with disorder
	SingleSolver single; //for other geometries, this would need to be promoted to a single-electron parent class
	
	//matvecs
	void explicitMultMv(ART *v, ART *w);
//	void MultMv(ART *v, ART *w);
	Eigen::Matrix<ART,-1,1> MultEigen(Eigen::Matrix<ART,-1,1>);
	void MatvecToDense();

	//lookup stuff
	vector< vector<int> > lookup_table_four,lookup_table_two;
	void make_lookups();
	
	//handy hamiltonian makign stuff
	int adjust_sign(int a,int b,int state);
	int adjust_sign(int a,int b, int c, int d, int state);
	int adjust_sign(int a,int b, int c, int d, int e, int f, int state);
	int hasbit(int i,int a);
	int state_to_index(int state);
	double V_Coulomb(double qx,double qy);
	
	//converting between truncated and original basis
	void plot_spectrum(string name);
	vector<int> orig_states;
	void basis_convert(vector<ART> &evec);
	void expand(int,ART,int, vector<ART> &new_evec);

};

///**************DEFINITIONS HERE************///

template<class ART>
ManySolver<ART>::ManySolver():MatrixWithProduct<ART>(),Wavefunction<ART>(){
	//read stuff for this run from the parameters file
	ifstream infile;
	infile.open("params");
	string temp("out");
	if (infile.is_open()){
		Ne=value_from_file(infile,-2);
		NPhi=value_from_file(infile,-2);
		charge=value_from_file(infile,-1);
		nLow=value_from_file(infile,0);
		nHigh=value_from_file(infile,0);
		outfilename=value_from_file(infile,temp);
		disorder_strength=value_from_file(infile,0.);
		V3overV1=value_from_file(infile,0.);
		NROD=value_from_file(infile,1);
		random_offset=value_from_file(infile,0);
		
		infile.close();
	}

	else cout << "Unable to open file"; 	
	
	this->init_wavefunction(NPhi);
	
	//set flags and do some simple sanity checks
	if(disorder_strength>0) disorder=1;
	else disorder=0;
	if(nHigh>0 || nLow>0) project=1;
	else project=0;
	if(disorder || project) cache=1;
	else cache=0;
	lookups=0;
	if (charge==-1) has_charge=0;
	else has_charge=1;
	if(project && has_charge){
		cout<<"you can't specify a charge if you want to project, so I removed the charge conservation"<<endl;
		has_charge=0;
		exit(0);
	}
	if(disorder && has_charge){
		cout<<"you can't specify a charge if you want disorder, so I removed the charge conservation"<<endl;
		has_charge=0;
		exit(0);
	}
	disordered_projection=0;
		
	HaldaneV=vector<double>(4,0);
	HaldaneV[1]=1/(1.+V3overV1);
	HaldaneV[3]=V3overV1/(1.+V3overV1);

	//you wouldn't believe how much faster this makes the code
	NPhi2=NPhi*NPhi;
	NPhi3=NPhi*NPhi2;
}	

template<class ART>
void ManySolver<ART>::init(){
	make_states(); 
	interaction_cache();
	if(project && !disordered_projection){
		projected_interaction_cache();
	}
}
////*********************SETUP FUNCTIONS
/// A state can be represented by a bit string, 
//ex. 7=11100000 has electrons in positions 0,1,2, 13=1011 has electron in 0,2,3
//this function starts at 0 and sees if its bitstring representation has the right number of electrons
//it it does, it adds it to a vector called states
template<class ART>
void ManySolver<ART>::make_states(){
	states.clear();
	int j=0,skip,count;
	for(int i=0;i<intpow(2,NPhi);i++){
		skip=0;
		//count number of bits in the integer
		count=0;
		for(int n=0;n<NPhi;n++){
			if(i& 1<<n) count++;
		}
		if (count==Ne && (has_charge==0 || get_charge(i)==charge) ){
			for(int k=0;k<nLow;k++) //these 4 lines are all that is needed for projection!
				if( !(i & 1<<k) ) skip=1;
			for(int k=NPhi-1;k>NPhi-1-nHigh;k--)
				if(i & 1<<k) skip=1;
			if(!skip){
				states.push_back(i);
				j++;
			}
		}
	}
	nStates=j;		 
	this->setrows(nStates);
	
	//this code is just a repeat of make_states, but it makes the untruncated basis
	if(project){
		for(int i=0;i<intpow(2,NPhi);i++){
			if (count_bits(i)==Ne){
				orig_states.push_back(i);
			}
		}
	}	
	cout<<"nStates: "<<nStates<<endl;	
//	for(int i=0;i<nStates;i++)
//		cout<<(bitset<9>)states[i]<<endl;
}
template<class ART>
void ManySolver<ART>::make_Hnn(){
	if(cache && disorder) disorder_cache();
	this->EigenDense=Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>::Zero(nStates,nStates);
	int j;
	ART temp;
	for(signed int a=0;a<NPhi;a++){

		//terms from the disorder potential
		if(disorder){
			for(int b=0;b<NPhi;b++){
//				cout<<"get disorder "<<a<<" "<<b<<"....";
				temp=get_disorder(a,b);
				for(int i=0;i<(signed)nStates;i++){
					if ((states[i] & 1<<a) && (b==a || (!( states[i] & 1<<b) && a>=nLow && a<NPhi-nHigh && b>=nLow && b<NPhi-nHigh) ) ){
						j=lookup_flipped(i,states,2,a,b);
						this->EigenDense(i,j)+=(double)adjust_sign(a,b,states[i]) * temp;
					}
				}
			}
		}

		for(signed int b=a+1;b<NPhi;b++){
						
			for(int c=0;c<NPhi;c++){
				for(int d=0;d<c;d++){
					temp=get_interaction(a,b,c,d);
					if(!project) //if we are not projecting into a different subspace, then this term conserves momentum and we can skip some elements in the sum
						if( (periodic && (a+b)%NPhi != (c+d)%NPhi) || (!periodic && a+b!=c+d) ) continue;
	//				cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<temp<<endl;
					for(int i=0;i<nStates;i++){
						if( (states[i] & 1<<a) && (states[i] & 1<<b) &&
						 ( (!(states[i] & 1<<c) && c<NPhi-nHigh && c>=nLow) || c==a || c==b) && 
						 ( (!(states[i] & 1<<d) && d>=nLow && d<NPhi-nHigh) || a==d || b==d) && 
						 ( (a>=nLow && a<NPhi-nHigh) || a==c ||a==d) &&
						 ( (b<NPhi-nHigh && b>=nLow) || b==c ||b==d)  )  {
							j=lookup_flipped(i,states,4,a,b,c,d);
							this->EigenDense(i,j)+=(double)(adjust_sign(a,b,c,d,states[i]) ) * temp;
						}
					}
				}//d	
			}//c	
		}//b
	}//a
}

//same as above, but does a six-body term for getting pfaffians
//note: this function is not compatible with projection or periodic BC
template<class ART>
void ManySolver<ART>::make_Hnn_six(){
	int f,j;
	ART temp;
	for(int a=0;a<NPhi;a++){
		for(int b=a+1;b<NPhi;b++){
			for(int c=b+1;c<NPhi;c++){
				for(int d=0;d<NPhi;d++){
					for(int e=d+1;e<NPhi;e++){
						f=a+b+c-d-e;
						if(f>=NPhi || f<=e) continue;
						temp=(six_body(a,b,c)-six_body(b,a,c)+six_body(b,c,a)-six_body(c,b,a)+six_body(c,a,b)-six_body(a,c,b));
						temp*=(six_body(d,e,f)-six_body(e,d,f)+six_body(e,f,d)-six_body(f,e,d)+six_body(f,d,e)-six_body(d,f,e));
						for(int i=0;i<nStates;i++){
							//cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<(bitset<7>)states[i]<<" "<<endl;
							if( (states[i] & 1<<a) && (states[i] & 1<<b) && (states[i] & 1<<c) && 
							(!(states[i] & 1<<d) || d==a || d==b || d==c) && 
							(!(states[i] & 1<<e) || e==a || e==b || e==c) && 
							(!(states[i] & 1<<f) || f==a || f==b || f==c) ){
								j=lookup_flipped(i,states,6,a,b,c,d,e,f);
								this->EigenDense(i,j)+=(double)(adjust_sign(a,b,c,d,e,f,states[i]) ) * temp;
							}
						}//for loop
					}//e
				}//d
			}//c
		}//b
	}//a
}
template<class ART> 
void ManySolver<ART>::interaction_cache(){
	int d;
	ART temp;
	four_body_cache.resize(pow(NPhi,4));
	
	for(int a=0;a<NPhi;a++){
		for(int b=a+1;b<NPhi;b++){

//			if(b>a) four_body_cache[four_array_map(a,b,b,a)]=two_body(a,b);		
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
				if ((a+b-c)<0) d=a+b-c+NPhi;
				else d=(a+b-c)%NPhi;
				if (d>=c) continue;
				four_body_cache[four_array_map(a,b,c,d)]=four_body(a,b,c,d);
			}	
		}
	}
}
	
template<>
inline void ManySolver< complex<double> >::disorder_cache(){
	complex<double> temp;
	two_body_cache.resize(pow(NPhi,2));
	
	for(signed int a=0;a<NPhi;a++){
		//terms from the disorder potential
		for(signed int b=0;b<NPhi;b++){		
			two_body_cache[four_array_map(a,b,0,0)]=single.getH(b,a);
		}
	}
}
template<>
inline void ManySolver< double >::disorder_cache(){
	cout<<"can't do disorder on a sphere!"<<endl;
	exit(0);
}
template<class ART>
void ManySolver<ART>::projected_interaction_cache(){
	ART temp;
	projected_four_body_cache.resize(pow(NPhi,4));
	
	for(int a=0;a<NPhi;a++){
		for(int b=a+1;b<NPhi;b++){

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
				for(int d=0;d<c;d++){
					projected_four_body_cache[four_array_map(a,b,c,d)]=four_body_project(a,b,c,d);
				}
			}	
		}
	}
}

////*****FUNCTIONS WHICH HELP WITH INTERACTION TERMS
template<>
inline double ManySolver<double>::get_disorder(int a,int b){
	cout<<"can't add disorder to a non-complex valued object"<<endl;
	exit(0);
}

//maps indices from a four array to indices from a 1 array
template<class ART>
int ManySolver<ART>::four_array_map(int a,int b,int c, int d){ return a+b*NPhi+c*NPhi2+d*NPhi3; }

template<>
inline complex<double> ManySolver< complex<double> >::four_body_project(int a,int b,int c,int d){
	complex<double> out=0;
	int id, sa,sb,sc,sd,signAB, signCD;
	for(int ia=0;ia<NPhi;ia++){
		for(int ib=0;ib<NPhi;ib++){
			if (ib==ia) continue;

			if(ib<ia){
				sa=ib; sb=ia; signAB=-1;
			}else{
				sa=ia; sb=ib; signAB=1;
			}
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
				if(id>ic){
					sc=id; sd=ic; signCD=-1;
				}else{
					sc=ic; sd=id; signCD=1;
				}
				out+=(double)(signAB*signCD)*single.evec(a,ia)*single.evec(b,ib)*conj(single.evec(c,ic))*conj(single.evec(d,id))*four_body_cache[four_array_map(sa,sb,sc,sd)];
			}
		}
	}
	return out;
}

template<>
inline double ManySolver<double>::four_body_project(int a, int b, int c, int d){
	cout<<"can't project in a non-complex object"<<endl;
	exit(0);
}
template<>
inline complex<double> ManySolver< complex<double> >::two_body_project(int a,int b){
	complex<double> out=0;
	for(int ia=0;ia<NPhi;ia++){
		for(int ib=0;ib<NPhi;ib++){
	//		cout<<ia<<" "<<ib<<" "<<single.evec(a+lowDeltas, ia)<<" "<<
			out+=single.evec(a,ia)*conj(single.evec(b,ib))*two_body_cache[four_array_map(ia,ib,0,0)];
		}
	}
//	cout<<a<<" "<<b<<" **************"<<out<<endl;
	return out;
}
template<>
inline double ManySolver<double>::two_body_project(int a, int b){
	cout<<"can't project in a non-complex object"<<endl;
	exit(0);
}

template<class ART>
ART ManySolver<ART>::get_interaction(int a,int b,int c, int d){
	if(!cache) return four_body(a,b,c,d);
	else if(!project) return four_body_cache[four_array_map(a,b,c,d)];
	else if(!disordered_projection) return projected_four_body_cache[four_array_map(a,b,c,d)];
	else return four_body_project(a,b,c,d);
}

//returns the disorder terms from a Hamiltonian
template<>
inline complex<double> ManySolver< complex<double> >::get_disorder(int a, int b){
	if(!cache) return single.getH(b,a);
	else if(!project) return two_body_cache[four_array_map(a,b,0,0)];
	else return two_body_project(a,b);
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
int ManySolver<ART>::adjust_sign(int a,int b,int state){
	int sign=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i) sign*=-1;
	return sign;	
}
//puts in minus signs for correct Fermi statistics
//for four-body terms, we have something like c^dagger_a c_d. This gives a (-1) for every electron between a and d
//this also works fine when you consider that we are performing two hops
template<class ART>
int ManySolver<ART>::adjust_sign(int a,int b,int c,int d,int state){
	int sign=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i && i!=c && i!=d) sign*=-1;
	start=c, end=d;
	if (c>d){ start=d; end=c;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i && i!=a && i!=b) sign*=-1;
	return sign;	
}
template<class ART>
int ManySolver<ART>::adjust_sign(int a, int b, int c, int d, int e, int f, int state){
	int sign=1;
	int start=a, end=b;
	if (a>b){ start=b; end=a;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i && i!=c && i!=d && i!=e && i!=f) sign*=-1;
	start=c, end=d;
	if (c>d){ start=d; end=c;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i && i!=a && i!=b && i!=e && i!=f) sign*=-1;
	start=e, end=f;
	if (e>f){ start=f; end=e;}
	for(int i=start+1;i<end;i++)
		if(state & 1<<i && i!=a && i!=b && i!=c && i!=d) sign*=-1;

	return sign;	
}

template<class ART>
ART ManySolver<ART>::six_body(int a, int b, int c){
	cout<<"six body terms not implemented in this geometry"<<endl;
	exit(0);
	return 0;
}
template<class ART>
void ManySolver<ART>::make_lookups(){
	vector<int> temp(NPhi3*NPhi,0);
	lookup_table_four=vector<vector<int> >(nStates,temp);
	vector<int> temp2(NPhi2,0);
	lookup_table_two=vector<vector<int> >(nStates,temp);
	for(int i=0;i<nStates;i++){
		for(int a=nLow;a<NPhi-nHigh;a++){
			for(int b=nLow;b<NPhi-nHigh;b++){
				lookup_table_two[i][four_array_map(a,b,0,0)]=lookup_flipped(i,states,2,a,b);
			}
			for(int b=a+1;b<NPhi-nHigh;b++){
				for(int c=nLow;c<NPhi-nHigh;c++){
					for(int d=nLow;d<c;d++){
						if(!project && (a+b)%NPhi != (c+d) % NPhi) continue;
						lookup_table_four[i][four_array_map(a,b,c,d)]=lookup_flipped(i,states,4,a,b,c,d);
					}
				}
			}
		}
	}
	lookups=1;
}
	
//////******MATVECS*****////
//template<class ART>
//void ManySolver<ART>::explicitMultMv(ART* v, ART* w){
//	cout<<"this function is broken!"<<endl;
//	exit(0);
////an old matvec using the matrix Hnn
////	double alpha=1., beta=0.;
////	cblas_zgemv(CblasColMajor, CblasNoTrans, nStates, nStates, &alpha, Hnn.data(), nStates, v, 1, &beta, w, 1);
//	vector<int> filled;
//	vector<int> empty;
//	int j;
//	for(int i=0;i<nStates;i++) w[i]=0;
//	for(int i=0;i<nStates;i++){ // i is the input state
//		//find all filled and empty terms in input state
//		filled.clear(); empty.clear();
//		for(int m=0;m<NPhi;m++){
//			if (states[i] & 1<<m) filled.push_back(m);
//			else empty.push_back(m);
//		}		
//		//loop through filled and empty states to generate disorder-tunnelling term
//		if(disorder){
//			for(int f=0;f<filled.size();f++){
//				for(int e=0;e<empty.size();e++){
//					if(filled[f]<nLow || empty[e] > NPhi-nHigh-1) continue;
//					j=lookup_flipped(i,filled[f],empty[e]);
//					w[j]+=(double)adjust_sign(filled[f],empty[e],states[i]) * final_two_body_cache[four_array_map(filled[f],empty[e],0,0)]*v[i];
//				}
//			}
//		}
//		//loop through filled states to generate disorder diagonal term
//		if(disorder){
//			for(int f=0;f<filled.size();f++){
//				w[i]+=final_two_body_cache[four_array_map(filled[f],filled[f],0,0)]*v[i];
//			}
//		}
//		//loop through all pairs of filled states to generate diagonal term
//		for(int f1=0;f1<filled.size();f1++){
//			for(int f2=f1+1;f2<filled.size();f2++){
//				w[i]+=final_four_body_cache[four_array_map(filled[f1],filled[f2],filled[f2],filled[f1])]*v[i];
//			}
//		}
//		//loop through filled states^2 + empty states to generate four-body one-particle tunnelling terms
//		int a,b,c,d;
//		if(project){
//			for(int f1=0;f1<filled.size();f1++){
//				if(filled[f1]<nLow) continue;
//				for(int f2=0;f2<filled.size();f2++){
//					if(f1==f2) continue;
//					for(int e=0;e<empty.size();e++){
//						if(empty[e] > NPhi-nHigh-1) continue;
//						j=lookup_flipped(i,filled[f1],empty[e]);
//						if(filled[f2]<filled[f1]){ a=filled[f2]; b=filled[f1]; }
//						else{ a=filled[f1]; b=filled[f2]; }
//						if(filled[f2]<empty[e]){ d=filled[f2]; c=empty[e]; }
//						else{ d=empty[e]; c=filled[f2]; }
//				
//						w[j]+=(double)adjust_sign(a,b,c,d,states[i])*final_four_body_cache[four_array_map(a,b,c,d)]*v[i];
//					}
//				}
//			}
//		}
//		//loop through filled states^2 empty states^2 to generate four body tunnelling
//		for(int f1=0;f1<filled.size();f1++){
//			if(filled[f1]<nLow) continue;
//			for(int f2=f1+1;f2<filled.size();f2++){
//				if(filled[f2]<nLow) continue;
//				for(int e1=0;e1<empty.size();e1++){
//					if(empty[e1]>NPhi-nHigh-1) continue;
//					for(int e2=e1+1;e2<empty.size();e2++){
//						if(empty[e2] > NPhi-nHigh-1) continue;
//						if(!project && (filled[f1]+filled[f2])%NPhi!=(empty[e1]+empty[e2])%NPhi) continue;
//						j=lookup_flipped(i,filled[f1],filled[f2],empty[e2],empty[e1]);
//						//cout<<i<<" "<<j<<" "<<filled[f1]<<" "<<filled[f2]<<" "<<empty[e1]<<" "<<empty[e2]<<" "<<final_four_body_cache[four_array_map(filled[f1],filled[f2],empty[e2],empty[e1])]<<endl;
//						w[j]+=(double)adjust_sign(filled[f1],filled[f2],empty[e2],empty[e1],states[i])*final_four_body_cache[four_array_map(filled[f1],filled[f2],empty[e2],empty[e1])]*v[i];
//					}
//				}
//			}
//		}
//	}
//		
//} //  MultMv.

//template<class ART>
//void ManySolver<ART>::MultMv(ART * v, ART * w){
////	for(int i=0;i<nStates;i++) w[i]=0;
////	for(int i=0;i<nStates;i++){
////		for(int j=0;j<nStates;j++)
////			w[i]+=v[j]*this->EigenDense(i,j);
////	}	
//	Eigen::Map <Eigen::Matrix<ART, Eigen::Dynamic, 1> > mapped_v(v,nStates);
////	Eigen::Matrix<ART,-1,1> mapped_v(n);
////	for(int i=0;i<n;i++) mapped_v(i)=v[i];	
//	mult_out=this->EigenDense*mapped_v;
//	//for(int i=0;i<n;i++) w[i]=out(i);	
//	Eigen::Map <Eigen::Matrix<ART, -1, 1> > (w,nStates,1)=mult_out; //using just out.data() fails for an unknown reason
//		
////	double alpha=1., beta=0.;
////	cblas_zgemv(CblasColMajor, CblasNoTrans, nStates, nStates, &alpha, this->EigenDense.data(), nStates, v, 1, &beta, w, 1);
//}

template<class ART>
Eigen::Matrix<ART,-1,1> ManySolver<ART>::MultEigen(Eigen::Matrix<ART,-1,1> invec){
	return this->EigenDense*invec;
}

//convert a wavefunction in the truncated-orbitals basis to one in the Landau basis
template<class ART>
void ManySolver<ART>::basis_convert(vector<ART> &evec){
	vector<ART> new_evec(orig_states.size(),0);

	for(int i=0;i<nStates;i++){
		expand(states[i],evec[i],0,new_evec);
	}
	evec=new_evec;
	
}
//a recursive function used to help with the above basis change
//it finds the first set bit in a state, and expands that bit in the old basis
//then it chops off that part of the bit and recurses
template< >
inline void ManySolver<complex<double> >::expand(int state, complex<double>  coeff, int orig_state, vector<complex<double> > &new_evec){
	vector<int>::const_iterator low;
	int largest_bit=-1;
	complex<double>  temp;
	//find largest set bit
	for(int i=NPhi;i>-1;i--){//speedup idea: start at previously found largest bit
		if(state & 1<<i) {
			largest_bit=i;
			break;
		}
	}
//	cout<<largest_bit<<endl;
	//if can't find largest set bit, return
	if(largest_bit==-1){
		//this code is basically a copy of lookup_flipped, but it acts on the original basis
//		cout<<coeff<<endl;
		low=lower_bound(orig_states.begin(),orig_states.end(),orig_state);
		if(low!=orig_states.end()){
			new_evec[low-orig_states.begin()]+=coeff;
		}else{
			cout<<"error in expand: "<<(bitset<30>)orig_state<<endl; exit(0);
		}
	}
	//else recurse
	else{
		for(int i=0;i<NPhi;i++){
			if(orig_state & 1<<i) continue; //if that bit is already set, then this state doesn't contribute and we are done
			
			//to get the sign right, remember that we are expanding bits from largest to smallest. So if the expanded bit is not the smallest, we need to do 
			//(-1)^(number of smaller bits)
			if(count_bits(orig_state%(1<<i))%2==1) temp=(complex<double> )(-1);
			else temp=(complex<double> )1;
			temp*=conj(single.evec(largest_bit,i));
	//		cout<<i<<" "<<largest_bit<<" "<<temp<<endl;
			expand(state%(1<<largest_bit),coeff*temp,orig_state| 1<<i,new_evec);
		}
	}
}
template<>
inline void ManySolver<double>::expand(int state, double coeff, int o, vector<double> &evec){
	cout<<"can't project in a non-complex object"<<endl;
	exit(0);
}

#endif
