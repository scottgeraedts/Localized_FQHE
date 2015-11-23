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
#include <string>
extern "C"{
#include "cblas.h"
}
#include "arscomp.h"

//will use std::bitset library to store bits, which needs to know the number of bits at compile time
//making this a bit too big is probably fine, some required values
using namespace std;

template <class ART>
class ManySolver{
public:
	ManySolver();
	void print_H();

	void ZeroHnn();
	void make_Hnn();
	void disorderHnn();
	
protected:
	int Ne,NPhi,nStates,NPhi2,NPhi3;
	int nHigh,nLow,NROD;
	int charge,has_charge,disorder,periodic,cache,project,lookups;
	double disorder_strength,V3overV1;
	vector<double> HaldaneV;
	string outfilename;

	vector<int > states;
	Eigen::Matrix<ART,Eigen::Dynamic, Eigen::Dynamic> Hnn;//ManySolver in Landau basis
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> > es;
	
	virtual ART two_body(int a,int b) =0;//could make virtual
	virtual ART four_body(int a,int b,int c,int d) =0;//could make virtual
	virtual int get_charge(int state)=0;//could make virtual

	void init(int, double, int, int, double, double);
	template<class T> T value_from_file(ifstream &infile, T def);
	void make_states();
	void make_cache();
	void second_cache();

	vector<ART> four_body_cache, two_body_cache;
	vector<ART> final_four_body_cache, final_two_body_cache;
	ART four_body_project(int,int,int,int);
	ART two_body_project(int,int);
	ART get_interaction(int,int,int,int);
	ART get_disorder(int,int);
	int four_array_map(int,int,int,int);
	
	//stuff to deal with disorder
	SingleSolver single; //for other geometries, this would need to be promoted to a single-electron parent class
	
	void MultMv(ART *v, ART *w);
	void Hnn_matvec(ART *v, ART *w);
	void MatvecToDense();

	vector< vector<int> > lookup_table_four,lookup_table_two;
	int lookup_flipped(int i,int a,int b,int c,int d);
	int lookup_flipped(int i,int a,int b);
	void make_lookups();
	
	void print_eigenstates();
	int adjust_sign(int a,int b,int state);
	int adjust_sign(int a,int b, int c, int d, int state);
	int hasbit(int i,int a);
	int state_to_index(int state);
	double V_Coulomb(double qx,double qy);
	
	double entanglement_entropy(const vector< complex<double> > &evec);
	void basis_convert(vector<ART> &evec);
	void expand(int,ART,int, vector<ART> &new_evec, const vector<int> &orig_states);

	//math functions
	double ClebschGordan(int,int,int);
	int count_bits(int);
	long int comb(int,int);
	double factorial(int,int);
	double factorial(int);
	int intpow(int,int);

};

///**************DEFINITIONS HERE************///

template<class ART>
ManySolver<ART>::ManySolver(){
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
		
		infile.close();
	}

	else cout << "Unable to open file"; 	

	cache=1;
	//set flags and do some simple sanity checks
	if(disorder_strength>0) disorder=1;
	else disorder=0;
	if(nHigh>0 || nLow>0) project=1;
	else project=0;
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
	HaldaneV=vector<double>(4,0);
	HaldaneV[1]=1/(1.+V3overV1);
	HaldaneV[3]=V3overV1/(1.+V3overV1);

	//you wouldn't believe how much faster this makes the code
	NPhi2=NPhi*NPhi;
	NPhi3=NPhi*NPhi2;
}

template<class ART> template<class T>
T ManySolver<ART>::value_from_file(ifstream &infile, T def){
	string line;
	stringstream convert;
	T out;
	if(getline(infile,line)){
		convert.str("");
		convert<<line;
		convert>>out;
		return out;
	}
	else return def;
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
				if(i & i<<k) skip=1;
			if(!skip){
				states.push_back(i);
				j++;
			}
		}
	}
	nStates=j;		 
//	cout<<"nStates: "<<nStates<<endl;	
//	for(int i=0;i<nStates;i++)
//		cout<<(bitset<9>)states[i]<<endl;
}
template<class ART>
void ManySolver<ART>::make_Hnn(){
	if(cache) make_cache();
	make_states();
	Hnn=Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>::Zero(nStates,nStates);
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
						j=lookup_flipped(i,a,b);
						Hnn(i,j)+=(double)adjust_sign(a,b,states[i]) * temp;
					}
				}
			}
		}

		for(signed int b=a+1;b<NPhi;b++){
//			if (b>a){
//				temp=get_interaction(a,b);
//				for(int i=0;i<nStates;i++)
//					if (states[i].test(a) && states[i].test(b)) Hnn(i,i)+=temp;
//			}
						
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
				//	cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<i<<" "<<endl;
							j=lookup_flipped(i,a,b,c,d);
							Hnn(i,j)+=(double)(adjust_sign(a,b,c,d,states[i]) ) * temp;
						}
					}
				}//d	
			}//c	
		}//b
	}//a
}

template<class ART>
void ManySolver<ART>::make_cache(){
	int d;
	ART temp;
	four_body_cache.resize(pow(NPhi,4));
	two_body_cache.resize(pow(NPhi,2));
	
	for(signed int a=0;a<NPhi;a++){
		//terms from the disorder potential
		if(disorder){
			for(signed int b=0;b<NPhi;b++){		
				two_body_cache[four_array_map(a,b,0,0)]=single.getH(b,a);
			}
		}
		for(signed int b=a+1;b<NPhi;b++){

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

template<class ART>
void ManySolver<ART>::second_cache(){
	ART temp;
	final_four_body_cache.resize(pow(NPhi,4));
	final_two_body_cache.resize(pow(NPhi,2));
	
	for(signed int a=0;a<NPhi;a++){
		//terms from the disorder potential
		if(disorder){
			for(signed int b=0;b<NPhi;b++){		
				final_two_body_cache[four_array_map(a,b,0,0)]=get_disorder(a,b);
			}
		}
		for(signed int b=a+1;b<NPhi;b++){

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
				for(int d=0;d<c;d++){
					if(!project && (a+b)%NPhi != (c+d)%NPhi ) continue;
					final_four_body_cache[four_array_map(a,b,c,d)]=get_interaction(a,b,c,d);
				}
			}	
		}
	}
}

////*****FUNCTIONS WHICH HELP WITH INTERACTION TERMS
template<class ART>
ART ManySolver<ART>::get_interaction(int a,int b,int c, int d){
	if(!cache) return four_body(a,b,c,d);
	else if(!project) return four_body_cache[four_array_map(a,b,c,d)];
	else return four_body_project(a,b,c,d);
}

//returns the disorder terms from a Hamiltonian
template<class ART>
ART ManySolver<ART>::get_disorder(int a, int b){
	if(!cache) return single.getH(b,a);
	else if(!project) return two_body_cache[four_array_map(a,b,0,0)];
	else return two_body_project(a,b);
}

//maps indices from a four array to indices from a 1 array
template<class ART>
int ManySolver<ART>::four_array_map(int a,int b,int c, int d){ return a+b*NPhi+c*NPhi2+d*NPhi3; }

template<class ART>
ART ManySolver<ART>::four_body_project(int a,int b,int c,int d){
	ART out=0;
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
template<class ART>
ART ManySolver<ART>::two_body_project(int a,int b){
	ART out=0;
	for(int ia=0;ia<NPhi;ia++){
		for(int ib=0;ib<NPhi;ib++){
	//		cout<<ia<<" "<<ib<<" "<<single.evec(a+lowDeltas, ia)<<" "<<
			out+=single.evec(a,ia)*conj(single.evec(b,ib))*two_body_cache[four_array_map(ia,ib,0,0)];
		}
	}
//	cout<<a<<" "<<b<<" **************"<<out<<endl;
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
void ManySolver<ART>::make_lookups(){
	vector<int> temp(NPhi3*NPhi,0);
	lookup_table_four=vector<vector<int> >(nStates,temp);
	vector<int> temp2(NPhi2,0);
	lookup_table_two=vector<vector<int> >(nStates,temp);
	for(int i=0;i<nStates;i++){
		for(int a=nLow;a<NPhi-nHigh;a++){
			for(int b=nLow;b<NPhi-nHigh;b++){
				lookup_table_two[i][four_array_map(a,b,0,0)]=lookup_flipped(i,a,b);
			}
			for(int b=a+1;b<NPhi-nHigh;b++){
				for(int c=nLow;c<NPhi-nHigh;c++){
					for(int d=nLow;d<c;d++){
						if(!project && (a+b)%NPhi != (c+d) % NPhi) continue;
						lookup_table_four[i][four_array_map(a,b,c,d)]=lookup_flipped(i,a,b,c,d);
					}
				}
			}
		}
	}
	lookups=1;
}
	
template<class ART>
int ManySolver<ART>::lookup_flipped(int i,int a,int b,int c,int d){
//given a bitstring, finds its index in the bitstring array
	if(lookups) return lookup_table_four[i][four_array_map(a,b,c,d)];
	else{
		int compare=states[i];
		compare=compare ^ 1<<a;
		compare=compare ^ 1<<b;
		compare=compare ^ 1<<c;
		compare=compare ^ 1<<d;
		vector<int>::iterator low;
		low=lower_bound(states.begin(),states.end(),compare);
		if(low!=states.end()) return (low-states.begin());
		else{
			cout<<"error in lookup_flipped: "<<(bitset<30>)states[i]<<" "<<(bitset<30>)compare<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
			exit(0);	
			return 0;
		}
	}
}

template<class ART>
int ManySolver<ART>::lookup_flipped(int i,int a,int b){
//given a bitstring, finds its index in the bitstring array
	if(lookups) return lookup_table_two[i][four_array_map(a,b,0,0)];
	else{
		int compare=states[i];
		compare=compare ^ 1<<a;
		compare=compare ^ 1<<b;
		vector<int>::iterator low;
		low=lower_bound(states.begin(),states.end(),compare);
		if(low!=states.end()) return (low-states.begin());
		else{
			cout<<"error in lookup_flipped: "<<(bitset<30>)states[i]<<" "<<(bitset<30>)compare<<" "<<a<<" "<<b<<endl;
			exit(0);
			return 0;	
		}
	}
}

////******MATVECS*****////
template<class ART>
void ManySolver<ART>::MultMv(ART* v, ART* w){
//an old matvec using the matrix Hnn
//	double alpha=1., beta=0.;
//	cblas_zgemv(CblasColMajor, CblasNoTrans, nStates, nStates, &alpha, Hnn.data(), nStates, v, 1, &beta, w, 1);
	vector<int> filled;
	vector<int> empty;
	int j;
	for(int i=0;i<nStates;i++) w[i]=0;
	for(int i=0;i<nStates;i++){ // i is the input state
		//find all filled and empty terms in input state
		filled.clear(); empty.clear();
		for(int m=0;m<NPhi;m++){
			if (states[i] & 1<<m) filled.push_back(m);
			else empty.push_back(m);
		}		
		//loop through filled and empty states to generate disorder-tunnelling term
		if(disorder){
			for(int f=0;f<filled.size();f++){
				for(int e=0;e<empty.size();e++){
					if(filled[f]<nLow || empty[e] > NPhi-nHigh-1) continue;
					j=lookup_flipped(i,filled[f],empty[e]);
					w[j]+=(double)adjust_sign(filled[f],empty[e],states[i]) * final_two_body_cache[four_array_map(filled[f],empty[e],0,0)]*v[i];
				}
			}
		}
		//loop through filled states to generate disorder diagonal term
		if(disorder){
			for(int f=0;f<filled.size();f++){
				w[i]+=final_two_body_cache[four_array_map(filled[f],filled[f],0,0)]*v[i];
			}
		}
		//loop through all pairs of filled states to generate diagonal term
		for(int f1=0;f1<filled.size();f1++){
			for(int f2=f1+1;f2<filled.size();f2++){
				w[i]+=final_four_body_cache[four_array_map(filled[f1],filled[f2],filled[f2],filled[f1])]*v[i];
			}
		}
		//loop through filled states^2 + empty states to generate four-body one-particle tunnelling terms
		int a,b,c,d;
		if(project){
			for(int f1=0;f1<filled.size();f1++){
				if(filled[f1]<nLow) continue;
				for(int f2=0;f2<filled.size();f2++){
					if(f1==f2) continue;
					for(int e=0;e<empty.size();e++){
						if(empty[e] > NPhi-nHigh-1) continue;
						j=lookup_flipped(i,filled[f1],empty[e]);
						if(filled[f2]<filled[f1]){ a=filled[f2]; b=filled[f1]; }
						else{ a=filled[f1]; b=filled[f2]; }
						if(filled[f2]<empty[e]){ d=filled[f2]; c=empty[e]; }
						else{ d=empty[e]; c=filled[f2]; }
				
						w[j]+=(double)adjust_sign(a,b,c,d,states[i])*final_four_body_cache[four_array_map(a,b,c,d)]*v[i];
					}
				}
			}
		}
		//loop through filled states^2 empty states^2 to generate four body tunnelling
		for(int f1=0;f1<filled.size();f1++){
			if(filled[f1]<nLow) continue;
			for(int f2=f1+1;f2<filled.size();f2++){
				if(filled[f2]<nLow) continue;
				for(int e1=0;e1<empty.size();e1++){
					if(empty[e1]>NPhi-nHigh-1) continue;
					for(int e2=e1+1;e2<empty.size();e2++){
						if(empty[e2] > NPhi-nHigh-1) continue;
						if(!project && (filled[f1]+filled[f2])%NPhi!=(empty[e1]+empty[e2])%NPhi) continue;
						j=lookup_flipped(i,filled[f1],filled[f2],empty[e2],empty[e1]);
						//cout<<i<<" "<<j<<" "<<filled[f1]<<" "<<filled[f2]<<" "<<empty[e1]<<" "<<empty[e2]<<" "<<final_four_body_cache[four_array_map(filled[f1],filled[f2],empty[e2],empty[e1])]<<endl;
						w[j]+=(double)adjust_sign(filled[f1],filled[f2],empty[e2],empty[e1],states[i])*final_four_body_cache[four_array_map(filled[f1],filled[f2],empty[e2],empty[e1])]*v[i];
					}
				}
			}
		}
	}
		
} //  MultMv.

template<class ART>
void ManySolver<ART>::Hnn_matvec(ART * v, ART * w){
	for(int i=0;i<nStates;i++) w[i]=0;
	for(int i=0;i<nStates;i++){
		for(int j=0;j<nStates;j++)
			w[i]+=v[j]*Hnn(i,j);
	}	
		
//	double alpha=1., beta=0.;
//	cblas_zgemv(CblasColMajor, CblasNoTrans, nStates, nStates, &alpha, Hnn.data(), nStates, v, 1, &beta, w, 1);
}
//turn the matvec into a dense matrix
template<class ART>
void ManySolver<ART>::MatvecToDense(){
	ART *dense=new ART[nStates*nStates];
	ART *v=new ART[nStates];
	ART *w=new ART[nStates];
	for(int i=0;i<nStates;i++){
		for(int j=0; j<nStates; j++){
			if(i==j || j==0) v[j]=1;
			else v[j]=0;
			w[j]=1;
		}
		Hnn_matvec(v,w);
		for(int j=0; j<nStates; j++){
			dense[i+j*nStates]=w[j];
			cout<<w[j]<<" ";
		}
		cout<<endl;
	}
	delete [] v;
	delete [] w;
	delete [] dense;
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

/////********MEASUREMENT FUNCTIONS*************////
template<class ART>
double ManySolver<ART>::entanglement_entropy(const vector< complex<double> > &evec){
	double out=0;
	int trace_start,trace_end;
	ofstream spectout;
	stringstream filename;
	vector<int> rows; 
	int nTrunc, trunc_part;
	vector<int> trunc_states;
	bool found;
//	for(int j=0;j<nStates;j++) cout<<(bitset<9>)states[j]<<" "<<(*evec)[j]<<endl;
	Eigen::Matrix<double,-1,-1> rhoblock;
	for(int orb=0;orb<NPhi;orb++){
		trace_start=orb; trace_end=(orb+NPhi/2)%NPhi; //which states to trace out

		//compute trunc_part
		trunc_part=0; //bitwise-AND with this gives just the component in the traced over states
		for(int i=0;i<NPhi;i++){
			if ((i>=trace_start && i <trace_end && trace_end>trace_start) || (trace_end<trace_start && (i < trace_end || i >= trace_start) ) )
				trunc_part=trunc_part | 1<<i;
		}

		//figure out how many reduced states there are		
		trunc_states.clear();
		nTrunc=0;
		for(int i=0;i<evec.size();i++){
			found=false;
			for(int j=0;j<trunc_states.size();j++)
				if( (states[i] & ~trunc_part) ==trunc_states[j]) found=true;
			if(!found) trunc_states.push_back(states[i] & ~trunc_part);
		}
		nTrunc=trunc_states.size();
		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rho=Eigen::Matrix<double,-1,-1>::Zero(nTrunc,nTrunc);
		//make matrix by looping over all states that aren't traced over
		int ti,tj;
		for(int i=0;i<evec.size();i++){
			for(int j=0;j<evec.size();j++){
				if( (states[i] & trunc_part) == (states[j] & trunc_part) ){
					for(ti=0;ti<nTrunc;ti++)
						if( (states[i] & ~trunc_part) == trunc_states[ti]) break;
					for(tj=0;tj<nTrunc;tj++)
						if( (states[j] & ~trunc_part) == trunc_states[tj]) break;
				
					rho(ti,tj)+=(evec[i]*conj(evec[j])).real();
				}
			}
		}
		//diagonalize matrix
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,-1,-1> > rs(rho);
		//output sum
		for(int i=0;i<nTrunc;i++) 
			if(rs.eigenvalues()(i)>0) out-=rs.eigenvalues()(i)*log(rs.eigenvalues()(i));
			
		///if the model has a charge, plot the entanglement spectrum			
		if(has_charge && periodic){
			filename.str("");
			filename<<"spectrum"<<orb;
			spectout.open(filename.str().c_str());			
			for(int c=0;c<NPhi;c++){
				for(int n=0;n<Ne;n++){
					//sort trunc_states based on their charge
					rows.clear();
					for(int j=0;j<nTrunc;j++)
						if(get_charge(trunc_states[j])==c && count_bits(trunc_states[j])==n) rows.push_back(j);
					if (rows.size()	== 0) continue;
	
					//use the same keys to get block-diagonal rho
					rhoblock=Eigen::Matrix<double,-1,-1>::Zero(rows.size(),rows.size());
					for(int j=0;j<rows.size();j++){
						for(int k=0;k<rows.size();k++){
							rhoblock(j,k)=rho(rows[j],rows[k]);
						}
					}
	
					//diagonalize the rho and print the results
					rs.compute(rhoblock);			
					for(int j=0;j<rows.size();j++) spectout<<c<<" "<<n<<" "<<rs.eigenvalues()(j)<<endl;
				}
			}
		}
		spectout.close();	
	}
	return out/(1.*NPhi);
}
//convert a wavefunction in the truncated-orbitals basis to one in the Landau basis
template<class ART>
void ManySolver<ART>::basis_convert(vector<ART> &evec){
	//this code is just a repeat of make_states, but it makes the untruncated basis
	vector<int> orig_states;
	for(int i=0;i<intpow(2,NPhi);i++){
		if (count_bits(i)==Ne){
			orig_states.push_back(i);
		}
	}
	vector<ART> new_evec(orig_states.size(),0);

	for(int i=0;i<nStates;i++){
		expand(states[i],evec[i],0,new_evec,orig_states);
	}
	evec=new_evec;
	
}
//a recursive function used to help with the above basis change
//it finds the first set bit in a state, and expands that bit in the old basis
//then it chops off that part of the bit and recurses
template<class ART>
void ManySolver<ART>::expand(int state, ART coeff, int orig_state, vector<ART> &new_evec, const vector<int> &orig_states){
	vector<int>::const_iterator low;
	int largest_bit=-1;
	int orig_state2;
	ART temp;
	//cout<<(bitset<6>)orig_state<<" "<<(bitset<6>)state<<" "<<coeff<<endl;
	//find largest set bit
	for(int i=NPhi;i>-1;i--){
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
			else orig_state2=orig_state | 1<<i;
			
			//to get the sign right, remember that we are expanding bits from largest to smallest. So if the expanded bit is not the smallest, we need to do 
			//(-1)^(number of smaller bits)
			if(count_bits(orig_state%(1<<i))%2==1) temp=(ART)(-1);
			else temp=(ART)1;
			temp*=conj(single.evec(largest_bit,i));
	//		cout<<i<<" "<<largest_bit<<" "<<temp<<endl;
			expand(state%(1<<largest_bit),coeff*temp,orig_state2,new_evec,orig_states);
		}
	}
}

//////////////////////utility functions, generally	
//someday should probably make a math library that contains all these functions

//counts the number of set bits in an integer
template<class ART>
int ManySolver<ART>::count_bits(int x){
	int out=0, i=0,found_bits=0;
	while(x!=found_bits){
		if(1<<i & x){
			out++;
			found_bits+=1<<i;
		}
		i++;
	}
	return out;
}
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
