#ifndef TORUS_SOLVER_H
#define TORUS_SOLVER_H
#include "version.h"
#include "ManySolver.h"
#ifdef USE_CLUSTER
#include "matprod2.h"
#endif

extern"C"{
	complex<double> landau_coulomb_(int *k1, int *k2, int *k3, int *k4);
	void make_landau_coulomb_();
	void coulomb_setup_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void setup_z_function_table_();
}

template <class ART>
class TorusSolver:public ManySolver<ART>{

public:

	TorusSolver(int);
	double numerical_semidefinite_integral(double dx,double start,double tol,double z);
	void add_disorder(int seed, double V,int,int);
	void structure_factors(const vector<ART> &eigvec);
	void berry_phase();
	void run_groundstate();
	void run_finite_energy(int offset);
	double self_energy();

	void makeShrinker(int nx);
	Eigen::SparseMatrix<ART> density_operator(int mx, int my);
	bool arpack,haldane,duncan;
	Eigen::SparseMatrix<ART> shrinkMatrix;
private:

	double Lx,Ly,V;
	int count;
	double Misra_onehalf(double t,double z);

	ART two_body(int a,int b);//could make virtual
	ART four_body(int a,int b,int c,int d);//could make virtual
	ART four_body_coulomb(int a,int b, int c, int d);
	ART four_body_haldane(int a,int b, int c, int d);
	ART four_body_duncan(int a, int b, int c, int d);
	double Hermite(double x,int n);
	int get_charge(int state);//could make virtual

};

template <class ART>
TorusSolver<ART>::TorusSolver(int x):ManySolver<ART>(){
	//stuff unique to the torus
	double alpha=1;
	Ly=sqrt(2.*M_PI*this->NPhi*alpha);//aspect ratio is Lx/Ly=alpha
	Lx=Ly/alpha;
	if(this->disorder || this->project){
		this->single=SingleSolver(this->NPhi,0,Lx,Ly);	
		if(this->project && !this->disordered_projection) this->single.init_deltas_lattice(this->nHigh);
	}
	
	this->init(x);
	this->periodic=1;
	arpack=true;

}

template<class ART>
void TorusSolver<ART>::run_finite_energy(int offset){

#ifndef USE_CLUSTER
	cout<<"can't use this function if not on the cluster!"<<endl;
	exit(0);
#endif
	haldane=true;
	this->store_sparse=false;

	//which states to look at
	double minE,maxE;
	vector<double> windows;//which energies to look at
	for(int i=1;i<4;i++) windows.push_back(i/(1.*4));
	double jindex;
	clock_t t;

	//counters
	Eigen::VectorXd ee=Eigen::VectorXd::Zero(windows.size());
	Eigen::VectorXd ee2=Eigen::VectorXd::Zero(windows.size());
	Eigen::VectorXd kltot=Eigen::VectorXd::Zero(windows.size());
	Eigen::VectorXd rtot=Eigen::VectorXd::Zero(windows.size());
	Eigen::VectorXd oldrtot=Eigen::VectorXd::Zero(windows.size());

	vector<double> tempvec;
	
	//stuff for density of states calculation (increase these numbers if disorder strength>30)
	double startE=-30, endE=30, dE=0.01;
	int ngrid=floor((endE-startE)/dE);
	tempvec=vector<double>(ngrid,0);
	vector< vector<double> >DOS(windows.size(),tempvec);

	int stop=100; //how many eigenstates to look at in each energy window
	
	if(!arpack)
		tempvec=vector<double>(this->nStates,0);	
	else
		tempvec=vector<double>(stop*windows.size(),0);	

	vector< vector<double> > energies(this->NROD,tempvec);
	double temp;
	vector<double>::iterator low;
cout<<this->disorder<<endl;
	stringstream filename;
	ofstream energyout,eigvecout;

	vector < vector<double> > EE_levels_all(windows.size());
	vector < vector < vector<double> > > EE_levels_storage(windows.size());
	
#ifdef USE_CLUSTER
	MatrixWithProduct2 mat2(this->nStates);
	mat2.set_mode("superLU");
#endif
	for(int i=0;i<this->NROD;i++){
		//construct hamiltonian
		if(this->project && this->disordered_projection) this->single.init_deltas_random(i+this->random_offset,this->nLow,this->nHigh);
		if(this->disorder) this->single.init_whitenoise(i+offset,this->disorder_strength);

		t=clock();	
		this->make_Hnn();
		cout<<"time to make Hamiltonian "<<((float)(clock()-t))/(1.*CLOCKS_PER_SEC)<<endl;
//		for(int m=0;m<this->nStates;m++){
//			for(int n=0;n<this->nStates; n++){
//				matrixout<<real(this->EigenDense(n,m))<<" ";
//				imagmatrixout<<imag(this->EigenDense(n,m))<<" ";
//			}
//			matrixout<<endl;
//			imagmatrixout<<endl;
//		}
//		break;
		if(!arpack){
			this->EigenDenseEigs();	
			//compute windows, and get js from windows
			minE=this->eigvals[0];
			maxE=this->eigvals[this->nStates-1];
			//cout<<maxE<<" "<<minE<<endl;
			for(int w=0;w<this->nStates;w++) cout<<this->eigvals[w]<<endl;
			for(int w=0;w<windows.size();w++){
				low=lower_bound(this->eigvals.begin(),this->eigvals.end(),minE+windows[w]*(maxE-minE));
				jindex=low-this->eigvals.begin();
				if(jindex> this->nStates-stop) cout<<"window too close to the end"<<endl;

			//compute observables in each window
				temp=this->entanglement_entropy(this->eigvecs,this->states,jindex);
				ee(w)+=temp;
				ee2(w)+=temp*temp;
				oldrtot(w)+=stupid_spacings(this->eigvals,jindex,jindex+stop,w);	
			}
			//average the density of states
			energies[i]=this->eigvals; //annoyingly, need to save all the eigenvalues for later
		}else{
#ifdef USE_CLUSTER		//if you're not on the cluster you don't have access to the libraries for this
			t=clock();	
			if(this->store_sparse) mat2.CSR_from_Sparse(this->EigenSparse);
			else mat2.CSR_from_Dense(this->EigenDense);
			cout<<"time to translate to a sparse matrix"<<((float)(clock()-t))/(CLOCKS_PER_SEC)<<endl;
			t=clock();
			time_t timer1,timer2;
			time(&timer1);
			maxE=mat2.single_energy("LR");
			minE=mat2.single_energy("SR");
			time(&timer2);
	                cout<<"time to get upper and lower"<<((float)(clock()-t))/(1.*CLOCKS_PER_SEC)<<" "<<difftime(timer2,timer1)<<endl;
			cout<<maxE<<" "<<minE<<endl;

			double eps;
			for(int w=0;w<windows.size();w++){
				eps=windows[w]*(maxE-minE)+minE;
				cout<<"eps="<<eps<<endl;
				mat2.eigenvalues(stop,eps);
	//			for(int k=0;k<mat2.eigvals.size();k++) cout<<mat2.eigvals[k]<<endl;
				temp_oldr=stupid_spacings(mat2.eigvals,w);
				oldrtot(w)+=temp_oldr;
				EE_levels_all[w].insert(EE_levels_all[w].end(),mat2.eigvals.begin(),mat2.eigvals.end());
				EE_levels_storage[w].push_back(mat2.eigvals);
				//save the eigenvalues
				filename.str("");
				filename<<"/mnt/cmcomp1/geraedts/7/energy"<<w<<"_"<<offset<<"_"<<this->disorder_strength;
				energyout.open(filename.str().c_str(),ios::app);
				for(int x=0;x<(signed)mat2.eigvals.size();x++) energyout<<mat2.eigvals[x]<<endl;
				energyout.close();

				//save the eigenvectors
/*				filename.str("");
				filename<<"/mnt/cmcomp1/geraedts/7/eigvec"<<w<<"_"<<offset<<"_"<<this->disorder_strength;
				eigvecout.open(filename.str().c_str(),ios::out|ios::binary|ios::app);
				for(int x=0;x<(signed)mat2.eigvecs.size();x++){
					for(int y=0;y<(signed)mat2.eigvecs[x].size();y++){
						eigvecout.write((char*)&(mat2.eigvecs[x][y]),sizeof(complex<double>));
					}
				}
				eigvecout.close();*/
				mat2.release_after_LU();
			
			}
#endif			
		}//if arpack
	}//NROD

	vector<double> energy_grid,integrated_DOS,s,s_spacings;
	ofstream rout,Sout;
	for(int w=0;w<(signed)windows.size();w++){
		filename.str("");
		filename<<"rout"<<w<<"_"<<offset;
		rout.open(filename.str().c_str());
		filename.str("");
		filename<<"sout"<<w<<"_"<<offset;
		Sout.open(filename.str().c_str());

		sort(EE_levels_all[w].begin(),EE_levels_all[w].end());
		energy_grid=make_grid(EE_levels_all[w],200);
		integrated_DOS=make_DOS(EE_levels_all[w],energy_grid);
		
		for(int i=0;i<this->NROD;i++){
			s=make_S(EE_levels_storage[w][i],energy_grid,integrated_DOS);
			rout<<compute_r(s)<<endl;
			s_spacings=spacings(s);
			for(int j=0;j<(signed)s_spacings.size();j++)
				Sout<<s_spacings[j]<<endl;
		}
		Sout.close();
		rout.close();
	}

//	write_vector(oldrtot,"oldr",this->NROD);
}

template<class ART>
void TorusSolver<ART>::run_groundstate(){

	haldane=false;
	arpack=false;
	this->disorder=0;
	this->project=0;
	this->NROD=1;
	int stop=25;
	Eigen::VectorXd energy_sum;
	if(!arpack) energy_sum=Eigen::VectorXd::Zero(this->nStates);
	else energy_sum=Eigen::VectorXd::Zero(stop);

//	cout<<"make duncans stuff"<<endl;
//	complex<double> L1(Lx/sqrt(2),0), L2(0,Ly/sqrt(2));
//	set_l_(&(this->NPhi), &L1, &L2);
//	setup_z_function_table_();
//	coulomb_setup_();
//	make_landau_coulomb_();
//	cout<<"duncan's matrix elements"<<endl;
//	cout<<four_body_duncan(0,3,2,1)<<" "<<four_body_coulomb(0,3,2,1)<<endl;	
//	cout<<four_body_duncan(1,4,3,2)<<" "<<four_body_coulomb(1,4,3,2)<<endl;	
//	cout<<four_body_duncan(1,6,5,2)<<" "<<four_body_coulomb(1,6,5,2)<<endl;	

	for(int i=0;i<this->NROD;i++){
		//construct hamiltonian
		if(this->project && this->disordered_projection) this->single.init_deltas_random(i+this->random_offset,this->nLow,this->nHigh);
		if(this->disorder) this->single.init_whitenoise(i+this->random_offset,this->disorder_strength);

		this->make_Hnn();
		
		if(!arpack){
			this->EigenDenseEigs();	
		}else{
		//****A call to ARPACK++. The fastest of all methods		
			this->eigenvalues(stop);			
		}//if arpack
		
		for(int j=0;j<this->eigvals.size();j++) energy_sum(j)+=this->eigvals[j]/(1.*this->NROD);
			 
	}//NROD
//	write_vector(energy_sum,"energies");
	for(int i=0;i<10;i++) cout<<this->eigvals[i]<<endl;
	//" "<<self_energy()<<" "<<this->eigvals[0]/(1.*this->Ne)+self_energy()<<endl;
//	haldane=false;
//	this->make_Hnn();
//	Eigen::VectorXcd tempvec=Std_To_Eigen(this->eigvecs[0]);
//	complex<double> ee=tempvec.adjoint()*this->EigenDense*tempvec;
//	cout<<"energy: "<<real(ee)<<" "<<self_energy()<<" "<<real(ee)/(1.*this->Ne)+self_energy()<<endl;
//	structure_factors(this->eigvecs[0]);

//	duncan=true;
//	this->make_Hnn();
//	tempvec=Std_To_Eigen(this->eigvecs[0]);
//	ee=tempvec.adjoint()*this->EigenDense*tempvec;
//	cout<<"energy: "<<real(ee)<<" "<<self_energy()<<" "<<real(ee)/(1.*this->Ne)+self_energy()<<endl;
//	structure_factors(this->eigvecs[0]);

//	stringstream filename;
//	filename<<"eigenvectors"<<this->Ne;
//	ofstream eigout(filename.str().c_str());
//	for(unsigned int j=0;j<this->eigvecs[0].size();j++)
//		eigout<<real(this->eigvecs[0][j])<<" "<<imag(this->eigvecs[0][j])<<endl;
}

template<class ART>
void TorusSolver<ART>::berry_phase(){
	haldane=false;
	arpack=false;
	this->disorder=1;
	this->project=0;

	double A=0.01;
	double side=sqrt(A);
	//set up locations of the holes
	ofstream sumout("err_v_step",ios::app);
//	int steps=2;
	int step_array[]={5,10,15,20};
	vector<Eigen::MatrixXcd> overlaps;
	vector<Eigen::VectorXcd> psi0(3),psi1(3),psi2(3);
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es;
	Eigen::MatrixXcd total;
	for(int steps_c=0;steps_c<4;steps_c++){
		int steps=step_array[steps_c];
		double step=side/(1.*steps);
		vector<double> holes_x,holes_y;
		for(int x=0;x<steps;x++){
			holes_x.push_back(x*step);
			holes_y.push_back(0);
		}
		for(int y=0;y<this->NROD*steps;y++){
			holes_x.push_back(side);
			holes_y.push_back(y*step);
		}
		for(int x=steps;x>0;x--){
			holes_x.push_back(x*step);
			holes_y.push_back(this->NROD*side);
		}
		for(double y=this->NROD*steps;y>0;y--){
			holes_x.push_back(0);
			holes_y.push_back(y*step);
		}
		int nds=holes_x.size();
		overlaps=vector<Eigen::MatrixXcd>(nds,Eigen::MatrixXcd::Zero(3,3));	
	
		stringstream filename;
		cout<<this->nStates<<endl;	
		for(int b=0;b<nds;b++){
			this->single.init_hole(holes_x[b],holes_y[b]);
			this->make_Hnn();
			this->EigenDenseEigs();
			for(int i=0;i<3;i++) psi1[i]=Std_To_Eigen(this->eigvecs[i]);

			if(b>0){
				for(int gs1=0;gs1<3;gs1++){
					for(int gs2=0;gs2<3;gs2++){
						overlaps[b](gs1,gs2)=psi1[gs1].dot(psi2[gs2]);
					}
				}
			}
			if(b==nds-1){
				for(int gs1=0;gs1<3;gs1++){
					for(int gs2=0;gs2<3;gs2++){
						overlaps[0](gs1,gs2)=psi0[gs1].dot(psi1[gs2]);
					}
				}
			}			
			psi2=psi1;
			if(b==0) psi0=psi1;
			//print the energies to make sure that they make sense
	//		for(int i=0;i<5;i++) cout<<this->eigvals[i]<<" ";
	//		cout<<endl<<endl;

		}
		total=Eigen::MatrixXcd::Identity(3,3);
		double sum;
		vector<double> running_phase(3,0);
		for(int b=0;b<nds;b++){
			total=overlaps[b]*total;
	//		cout<<overlaps[b]<<endl;
			cout<<setprecision(15)<<holes_x[b]<<" "<<holes_y[b]<<endl;
			es.compute(overlaps[b]);
			sum=0;
			for(int i=0;i<3;i++){
				cout<<abs(es.eigenvalues()(i))<<" "<<arg(es.eigenvalues()(i))<<endl;
				sum+=arg(es.eigenvalues()(i));
				running_phase[i]+=arg(es.eigenvalues()(i));
			}
			cout<<"sum: "<<sum<<endl;
		}
		cout<<"total:"<<endl;
		cout<<total<<endl;
		es.compute(total);
		sum=0;
		double sum2=0;
		for(int i=0;i<3;i++){
			cout<<abs(es.eigenvalues()(i))<<" "<<arg(es.eigenvalues()(i))<<" "<<running_phase[i]<<endl;
			sum+=arg(es.eigenvalues()(i));
			sum2+=running_phase[i];
		}		
		cout<<"sum: "<<sum<<" "<<sum2<<endl;
		while(sum>M_PI) sum-=2*M_PI;
		while(sum<-M_PI) sum+=2*M_PI;
		sumout<<steps<<" "<<step<<" ";
		for(int i=0;i<3;i++) sumout<<abs(es.eigenvalues()(i))<<" "<<arg(es.eigenvalues()(i))<<" ";
		sumout<<total.norm();
		sumout<<endl;
	}
	sumout.close();
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
	
	for(int i=0;i<this->HaldaneV.size();i++){
		if(this->HaldaneV[i]!=0){
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
			bigout+=this->HaldaneV[i]*out;
			out=0;
		}
	}		
	return 2*bigout;//here swapping c&d gives the same thing up to a - sign that accounts for fermion parity, so just 
}
//uses duncan's library to generate coulomb matrix elements
//template<class ART>
//ART TorusSolver<ART>::four_body_duncan(int a, int b, int c, int d){
//	return landau_coulomb_(&a,&b,&c,&d);
//}

template<class ART>
ART TorusSolver<ART>::four_body(int a, int b, int c, int d){ 
	if(haldane==true) return four_body_haldane(a,b,c,d);
	//else if(duncan==true) return four_body_duncan(a,b,c,d);
	else return four_body_coulomb(a,b,c,d);	
}

//calculates structure factors
//right now only does qx=0 case since that is easiest
template<class ART>
void TorusSolver<ART>::structure_factors(const vector<ART> &eigvec){
	ofstream sqout("sq");
	vector<double> sk(this->NPhi), sk2(this->NPhi); 
	for(int i=0;i<this->nStates;i++){
		for(int k=0;k<this->NPhi;k++){
			if(bittest(this->states[i],k)) sk[k]+=norm(eigvec[i]);
			for(int n=0;n<this->NPhi;n++){
				if( bittest(this->states[i],n) && bittest(this->states[i],(n+k)%this->NPhi)) sk2[k]+=norm(eigvec[i]);
			}
		}
	}
	for(int i=0;i<this->NPhi;i++) sqout<<i<<" "<<sk[i]<<" "<<sk2[i]/(1.*this->NPhi)<<endl;
	sqout.close();		
}

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
