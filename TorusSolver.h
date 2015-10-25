#ifndef TORUS_SOLVER_H
#define TORUS_SOLVER_H
#include "ManySolver.h"

template <class ART>
class TorusSolver:public ManySolver<ART>{

public:

	TorusSolver(int Ne,int charge, double V);
	void make_Hnn();
	double numerical_semidefinite_integral(double dx,double start,double tol,double z);
	void add_disorder(int seed, double V);

private:

	double Lx,Ly,V;
	int count;
	double self_energy();
	double Misra_onehalf(double t,double z);

	ART two_body(int a,int b);//could make virtual
	ART four_body(int a,int b,int c,int d);//could make virtual
	int get_charge(bitset<NBITS> state);//could make virtual

};

template <class ART>
TorusSolver<ART>::TorusSolver(int tNe,int tcharge, double V):ManySolver<ART>(tNe,tcharge,1){
	//stuff unique to the torus

	this->NPhi=3*this->Ne;
	Ly=sqrt(M_PI*this->NPhi*this->Ne/2.);//aspect ratio is Lx/Ly=Ne/4
	Lx=(4/(1.*this->Ne))*Ly;

	//add_disorder(0,0.1);
	clock_t t;
	t=clock();
	this->init();
	cout<<this->es.eigenvalues().head(10).array()/this->Ne+self_energy()<<endl;
	t=clock()-t;
	cout<<"time: "<<((float)t)/CLOCKS_PER_SEC<<endl;

	this->cache=1;
	t=clock();
	this->init();
	cout<<this->es.eigenvalues().head(10).array()/this->Ne+self_energy()<<endl;
	t=clock()-t;
	cout<<"time: "<<((float)t)/CLOCKS_PER_SEC<<endl;
//	cout<<this->Hnn<<endl;

	this->basis=Eigen::Matrix<ART, -1, -1>::Identity(this->NPhi, this->NPhi);
	this->basis(0,0)=1/sqrt(2);
	this->basis(1,0)=-1/sqrt(2);
	this->basis(0,1)=1/sqrt(2);
	this->basis(1,1)=1/sqrt(2);
	
	this->project=1;
	t=clock();
	this->init();
	cout<<this->es.eigenvalues().head(10).array()/this->Ne+self_energy()<<endl;
	t=clock()-t;
	cout<<"time: "<<((float)t)/CLOCKS_PER_SEC<<endl;
//	cout<<this->Hnn<<endl;
	
//	count=0;
//	Eigen::VectorXd evals;
//	int N=10;
//	Eigen::VectorXd sum=Eigen::VectorXd::Zero(N);
//	int NROD=1;
//	for(int i=0;i<NROD;i++){
//		add_disorder(i,V);		
//		this->init();
//		sum.array()+=this->es.eigenvalues().head(N).array();
//		this->single.switchVstrength();
//		this->init();
//		sum.array()+=this->es.eigenvalues().head(N).array();
//	}
//	sum/=(NROD*2.*this->Ne);
//	sum.array()+=self_energy();
//	cout<<V<<" "<<sum.transpose()<<endl;

//	Eigen::VectorXd evec=es.eigenvectors().col(0);
//	for(int i=0;i<nStates;i++) cout<<states[i]<<" "<<evec(i)<<endl;

//	cout<<self_energy()<<endl;
	//cout<<count<<endl;
}
///functions which are definitely unique to the torus
template <class ART>
int TorusSolver<ART>::get_charge(bitset<NBITS> s){
	int out=0;
	for(unsigned int i=0;i<NBITS;i++){
		if (s.test(i)) out+=i;
	}
	return out%this->NPhi;
}
//compute two-body terms for torus case
template <class ART>
ART TorusSolver<ART>::two_body(int a,int b){
	//return 0;
	ART out=0.; double tol=1e-17;
	double qy,V,expqy,expp;
	double kdiff=2*M_PI/Lx*(a-b);
	int qfactor=1;
	for(int my=0;my>-1;my++){
		qy=2*M_PI*my/Ly;
		if (my>0) qfactor=2;
		expqy=exp(-0.5*qy*qy);
		if (expqy<tol) break;
//		cout<<expqy<<endl;
		//piece for a=n1=n4, b=n2=n3
		for(int p=0;p>-1;p++){
			if(p==0 && my==0) continue;
			expp=exp(-0.5*p*p*Ly*Ly);
			if(expqy*expp<tol) break;

			if(p==0) V=this->V_Coulomb(0,qy);
			else V=2*this->V_Coulomb(p*Ly,qy);

			out+=qfactor*cos(qy*kdiff)*expqy*expp*V;
			//cout<<p<<" "<<expp<<" "<<qfactor*cos(qy*kdiff)*expqy*expp*V<<endl;
		}

		//piece for a=n1=n3, b=n2=n4, non-negative p
		for(int p=1;p>-1;p++){
			expp=exp(-0.5*pow(p*Ly+kdiff,2));
			if (expp*expqy<tol) break;
			out-=qfactor * expqy*expp* this->V_Coulomb(p*Ly+kdiff,qy);
			//cout<<p<<" "<<expp<<" "<<qfactor*cos(qy*kdiff)*expqy*expp*V_Coulomb(p*Ly+kdiff,qy)<<endl;			
		}
		//piece for n1=n3, negative p
		for(int p=0;p<1;p--){
			expp=exp(-0.5*pow(p*Ly+kdiff,2));
			if (expp*expqy<tol) break;
			out-=qfactor * expqy*expp* this->V_Coulomb(p*Ly+kdiff,qy);
			//cout<<p<<" "<<expp<<" "<<qfactor*cos(qy*kdiff)*expqy*expp*V_Coulomb(p*Ly+kdiff,qy)<<endl;			
		}
	}
	return out;	
}
//for the torus case electrons interact with their own images, leading to an energy offset given in Yoshioka 2.7
//the answer involves computing an infinite sum of numerical integrals, so there are three numerical cutoffs that need to be implemented
//1: how close to infinity to set the upper bound of the integral
//	I don't know how to determine what the error will be from a given cutoff, but it takes O(-log(cutoff)) time to do this, so I'll just set it to machine precision
//	In practice this seems to give an error ~1e-10, so no point in other tolerances being smaller than this
//2: mesh of numerical integration, the error from a given mesh is O(z*e^{-z}dx^2), so can calculate dx given z and tol
//3: where to cut off the infinite sum. Each term is less than e^-z/z, so keep only terms with e^-z/z>tol
template <class ART>
double TorusSolver<ART>::self_energy(){
	double tol=1e-10;
	double out=0,z,temp;
	for(int l1=0;l1>-1;l1++){
		if (exp(-l1*l1*M_PI)/(M_PI*l1*l1)<tol) break;
		for(int l2=0;l2>-1;l2++){
			if (l1==0 && l2==0) continue;
			z=M_PI*(l1*l1*Lx/Ly+l2*l2*Ly/Lx);
			if (exp(-z)/z<tol) break;
			temp=numerical_semidefinite_integral(1.,sqrt(tol*exp(z)/z),1e-17,z);
			if (l1==0 || l2==0) out+=2*temp;
			else out+=4*temp;
		}
	}	
	return -(2-out)/sqrt(Lx*Ly);
}
//compute four-body terms for torus case
template <class ART>
ART TorusSolver<ART>::four_body(int a,int b,int c,int d){
	int print=0,start=0;
	ART out=0; double tol=1e-17;
	double qx,qy,expqy,expqx,temp,qfactor;
	
	for(int my=0;my>-1;my++){
		if(my==0) qfactor=1.;
		else qfactor=2.;
		
		qy=2*M_PI/Ly*my;
		expqy=exp(-0.5*qy*qy);
		if (expqy<tol) break;
		//a=n1, d=n4, b=n2, c=n3 term, note that p's are done seperately
		if (a<d) start=1; 
		else start=0; //depending on the sign of (a-d), the largest contribution comes from p=1, 0 or -1. start makes sure 1 or -1 is always included
		for(int p1=start;p1>-1;p1++){
			qx=2*M_PI/Lx*(a-d)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			//p2=((c-b)-(a-d)-p1*this->NPhi)/this->NPhi; //don't need this since i'm not including theta_y terms
			temp=qfactor*expqx*expqy* cos(qy*2*M_PI/Lx*(a-c))*this->V_Coulomb(qx,qy);//used V(qx,qy)=V(qx,-qy), this would need to be more complicated if that weren't true
			out-=temp;
			if(print) cout<<"p+ ad "<<p1<<" "<<2*M_PI/Lx*(a-d)<<" "<<p1*Ly<<" "<<temp<<endl;
		}
		for(int p1=start-1;p1<1;p1--){
			qx=2*M_PI/Lx*(a-d)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor*expqx*expqy* cos(qy*2*M_PI/Lx*(a-c))*this->V_Coulomb(qx,qy);
			out-=temp;
			if(print) cout<<"p- ad "<<p1<<" "<<temp<<endl;
		}			

		//a=n1,b=n2,c=n4,d=n3 term
		if(a<c) start=1;
		else start=0;
		for(int p1=start;p1>-1;p1++){
			qx=2*M_PI/Lx*(a-c)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor * expqx*expqy* cos(qy*2*M_PI/Lx*(a-d))* this->V_Coulomb(qx,qy);
			out+=temp;
			if(print) cout<<"p+ ac "<<p1<<" "<<temp<<endl;
		}
		for(int p1=start-1;p1<1;p1--){
			qx=2*M_PI/Lx*(a-c)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor * expqx*expqy* cos(qy*2*M_PI/Lx*(a-d)) * this->V_Coulomb(qx,qy);
			out+=temp;
			if(print) cout<<"p- ac "<<p1<<" "<<temp<<endl;
		}			
		
	}
	return out;
}
//loads single particle class that provides disorder potentials
template<class ART>
void TorusSolver<ART>::add_disorder(int seed, double V){
	if(this->has_charge){
		cout<<"can't add disorder to a system with conserved charge!"<<endl;
		exit(0);
	}
	this->disorder=1;
	this->single=SingleSolver(this->NPhi,0);
	this->single.init(seed,V);
}

//right now this calculates the integral of a Misra function (though you can change the code to make it do any integral from start to infinity)
//if a general solver is desired would need function pointers, and then z would be passed to the called function, not this function
//tol tells you how small a term can get before you stop including terms
template <class ART>
double TorusSolver<ART>::numerical_semidefinite_integral(double start,double dx,double inf_tol,double z){
	double out=0,last;
	int i=0;
	while(1){
		last=dx*0.5*(Misra_onehalf( start+dx*i ,z) +Misra_onehalf( start+dx*(i+1) ,z) );
		out+=last;
		if (abs(last)<inf_tol) break;
		i++;
	}
	return out;
}	
template <class ART>
double TorusSolver<ART>::Misra_onehalf(double t,double z){ return exp(-z*t)/sqrt(t);}		


#endif
