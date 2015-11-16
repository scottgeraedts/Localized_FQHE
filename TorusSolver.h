#ifndef TORUS_SOLVER_H
#define TORUS_SOLVER_H
#include "ManySolver.h"
#include "arscomp.h"

template <class ART>
class TorusSolver:public ManySolver<ART>{

public:

	TorusSolver(int Ne,int charge, double V, int invfilling, string name,int);
	void make_Hnn();
	double numerical_semidefinite_integral(double dx,double start,double tol,double z);
	void add_disorder(int seed, double V,int,int);

private:

	double Lx,Ly,V;
	vector<double> HaldaneV;
	int count;
	double self_energy();
	double Misra_onehalf(double t,double z);

	ART two_body(int a,int b);//could make virtual
	ART four_body(int a,int b,int c,int d);//could make virtual
	ART four_body_coulomb(int a,int b, int c, int d);
	ART four_body_haldane(int a,int b, int c, int d);
	double Hermite(double x,int n);
	int get_charge(int state);//could make virtual

};

template <class ART>
TorusSolver<ART>::TorusSolver(int tNe,int tcharge, double V, int _NPhi, string name, int N):ManySolver<ART>(tNe,tcharge,1){
	//stuff unique to the torus
	this->NPhi=_NPhi;
	Ly=sqrt(M_PI*this->NPhi*this->Ne/2.);//aspect ratio is Lx/Ly=Ne/4
	Lx=(4/(1.*this->Ne))*Ly;

	HaldaneV=vector<double>(4,0);
	HaldaneV[1]=cos(V);
	HaldaneV[3]=sin(V);
	double se=self_energy();
	double ee=0;
	int stop=15,NROD=10;	
	Eigen::VectorXd sum=Eigen::VectorXd::Zero(stop);
//	cout<<"run with Ne="<<this->Ne<<endl;
	stringstream filename;
	filename.str("");
	filename<<"V3_1good/"<<name<<"_"<<this->Ne;
	ofstream cfout;
	if(tcharge==-1) cfout.open(filename.str().c_str());
	else cfout.open(filename.str().c_str(),ofstream::app);

	for(int i=0;i<NROD;i++){
		this->init(i,V,0,N,Lx,Ly);
		if(this->nStates<=stop) stop=this->nStates-1;
	
//		this->es.compute(this->Hnn);
//		sum=this->es.eigenvalues();
//		sum=sum/(1.*this->Ne);
//		sum=sum.array()+se;
//		for(int j=0;j<stop;j++) cfout<<tcharge<<" "<<sum(j)<<endl;

	//****A call to ARPACK++. The fastest of all methods		
//		clock_t t;
//		t=clock();
		ARCompStdEig<double, TorusSolver<ART> >  dprob(this->nStates, stop, this, &TorusSolver<ART>::Hnn_matvec,"SR",(int)0, 1e-10,1e6);//someday put this part into matprod?
		dprob.FindEigenvalues();
		dprob.FindEigenvectors();
//		cout<<"old time: "<<(float)(clock()-t)/CLOCKS_PER_SEC<<endl;
		
		for(int i=0;i<dprob.ConvergedEigenvalues();i++) sum(i)+=dprob.Eigenvalue(i).real()/(1.*this->Ne)+se;
		ee+=this->entanglement_entropy(dprob.StlEigenvector(0));
	}
	sum/=(1.*NROD);
	cfout<<sum<<endl;

//	for(int i=0;i<stop;i++) cfout<<tcharge<<" "<<sum(i)<<endl;

	cfout.close();

}
///functions which are definitely unique to the torus
template <class ART>
int TorusSolver<ART>::get_charge(int s){
	int out=0;
	for(unsigned int i=0;i<this->NPhi;i++){
		if (s & 1<<i) out+=i;
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
ART TorusSolver<ART>::four_body_coulomb(int a,int b,int c,int d){
	int print=0,start=0;
	ART out=0; double tol=1e-17;
	double qx,qy,expqy,expqx,temp,qfactor;
	
	for(int my=0;my>-1;my++){
		if(my==0) qfactor=1.;
		else qfactor=2.;
		
		qy=2*M_PI/Ly*my;
		expqy=exp(-0.5*qy*qy);
		if (expqy<tol) break;

		//a=n1,b=n2,c=n4,d=n3 term
		if(a<c) start=1;
		else start=0;
		for(int p1=start;p1>-1;p1++){
			qx=2*M_PI/Lx*(a-c)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor * expqx*expqy* cos(qy*2*M_PI/Lx*(a-d))* this->V_Coulomb(qx,qy);
			out-=temp;
			if(print) cout<<"p+ ac "<<p1<<" "<<temp<<endl;
		}
		for(int p1=start-1;p1<1;p1--){
			qx=2*M_PI/Lx*(a-c)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor * expqx*expqy* cos(qy*2*M_PI/Lx*(a-d)) * this->V_Coulomb(qx,qy);
			out-=temp;
			if(print) cout<<"p- ac "<<p1<<" "<<temp<<endl;
		}		
		
		if(a==d && my==0) continue;

		//a=n1, d=n4, b=n2, c=n3 term, note that p's are done seperately
		if (a<d) start=1; 
		else start=0; //depending on the sign of (a-d), the largest contribution comes from p=1, 0 or -1. start makes sure 1 or -1 is always included
		for(int p1=start;p1>-1;p1++){
			qx=2*M_PI/Lx*(a-d)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			//p2=((c-b)-(a-d)-p1*this->NPhi)/this->NPhi; //don't need this since i'm not including theta_y terms
			temp=qfactor*expqx*expqy* cos(qy*2*M_PI/Lx*(a-c))*this->V_Coulomb(qx,qy);//used V(qx,qy)=V(qx,-qy), this would need to be more complicated if that weren't true
			out+=temp;
			if(print) cout<<"p+ ad "<<p1<<" "<<2*M_PI/Lx*(a-d)<<" "<<p1*Ly<<" "<<temp<<endl;
		}
		for(int p1=start-1;p1<1;p1--){
			qx=2*M_PI/Lx*(a-d)+p1*Ly;
			expqx=exp(-0.5*qx*qx);
			if (expqx*expqy<tol) break;			
			temp=qfactor*expqx*expqy* cos(qy*2*M_PI/Lx*(a-c))*this->V_Coulomb(qx,qy);
			out+=temp;
			if(print) cout<<"p- ad "<<p1<<" "<<temp<<endl;
		}			

			
		
	}
	return out;
}

template<class ART>
ART TorusSolver<ART>::four_body_haldane(int a, int b, int c, int d){
	int m=a-d, k=a-c; //these are the indices defined in Roger and Mike's paper they are the two variables left after COM and momentum conservation
	//V=H((m+k)/2)H((m-k)/2), but need to take the PBC into account
	//I'm using the definition in Zlatko's paper, which is off by Mikes definition by a mysterious factor of 1/(2*I)
	double bigout=0,out=0,expm,expk,qm,qk,tol=1e-12;
	int startm=0, startk=0;
	if(m<0) startm=1;
	if(k<0) startk=1;
	
	for(int i=0;i<HaldaneV.size();i++){
		if(HaldaneV[i]!=0){
			for(int pm=startm; pm>-1;pm++){
				qm=2.*M_PI*(m+pm*this->NPhi)/Lx;
				expm=exp(-0.5*qm*qm);
				if (expm<tol) break;
				for(int pk=startk; pk>-1;pk++){
					qk=2.*M_PI*(k+pk*this->NPhi)/Lx;
					expk=exp(-0.5*qk*qk);
					if(expk<tol) break;
					out+=expm*expk*Hermite(qm+qk,i)*Hermite(qk-qm,i);
				}
				for(int pk=startk-1; pk<1;pk--){
					qk=2.*M_PI*(k+pk*this->NPhi)/Lx;
					expk=exp(-0.5*qk*qk);
					if(expk<tol) break;
					out+=expm*expk*Hermite(qm+qk,i)*Hermite(qk-qm,i);
				}
			}
			for(int pm=startm-1; pm<1;pm--){
				qm=2.*M_PI*(m+pm*this->NPhi)/Lx;
				expm=exp(-0.5*qm*qm);
				if (expm<tol) break;
				for(int pk=startk; pk>-1;pk++){
					qk=2.*M_PI*(k+pk*this->NPhi)/Lx;
					expk=exp(-0.5*qk*qk);
					if(expk<tol) break;
					out+=expm*expk*Hermite(qm+qk,i)*Hermite(qk-qm,i);
				}
				for(int pk=startk-1; pk<1;pk--){
					qk=2.*M_PI*(k+pk*this->NPhi)/Lx;
					expk=exp(-0.5*qk*qk);
					if(expk<tol) break;
					out+=expm*expk*Hermite(qm+qk,i)*Hermite(qk-qm,i);
				}
			}
			bigout+=HaldaneV[i]*out;
			out=0;
		}
	}		
	return 2*bigout;//here swapping c&d gives the same thing up to a - sign that accounts for fermion parity, so just 
}
template<class ART>
ART TorusSolver<ART>::four_body(int a, int b, int c, int d){ return four_body_haldane(a,b,c,d);}
//Hermite polynomials (using the 'probabalist' definition, see the wikipedia entry)
template<class ART>
double TorusSolver<ART>::Hermite(double x, int n){
	if(n==0)
		return 0;
	else if (n==1) return x;
	else if(n==2) return x*x-1;
	else if(n==3) return x*x*x-3*x;
	else{
		cout<<"asked for too large a Hermite polynomial, you should define more!"<<endl;
		exit(0);
	}
}

//loads single particle class that provides disorder potentials
template<class ART>
void TorusSolver<ART>::add_disorder(int seed, double V, int nLow, int nHigh){
	if(this->has_charge){
		cout<<"can't add disorder to a system with conserved charge!"<<endl;
		exit(0);
	}
	this->disorder=1;
	this->single=SingleSolver(this->oldNPhi,0);
	this->single.init(seed,V,nLow,nHigh);
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
