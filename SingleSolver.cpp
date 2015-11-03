#include "SingleSolver.h"

SingleSolver::SingleSolver(int _NPhi, int N_deltas, double _Lx, double _Ly, int _NROD, int _shift, double _Vstrength):NPhi(_NPhi),NROD(_NROD),shift(_shift),Vstrength(_Vstrength),Lx(_Lx), Ly(_Ly){
	if(Lx==0 || Ly==0){
		Lx=sqrt(2*M_PI*NPhi);//square aspect ratio, might want to change this later
		Ly=Lx;
	}
	qbounds=(CUTOFF+1)*NPhi-1;
	N_trunc=NPhi-2*N_deltas;
	energy_mesh=100+1;
	double meshstart=0.2/sqrt(1.*NPhi);
	energies=Eigen::VectorXd::LinSpaced(energy_mesh,-0.1*meshstart,0.1*meshstart);	
	spacings=Eigen::VectorXd::Zero(energy_mesh);
}

//the following two functions are designed to be called by a ManySolver
void SingleSolver::init(int seed,double V, int nLow, int nHigh){
	Vstrength=V;
	map <string, double> params;
	params["Vstrength"]=1.;
	params["nLow"]=nLow;
	params["nHigh"]=nHigh;
	Potential p_deltas("delta",qbounds,seed,Lx,Ly,params);
	make_Hamiltonian(p_deltas);
	Eigen::MatrixXcd old=Hnn;
	proj.compute(Hnn);
	params["Vstrength"]=1.;
	Potential p_gauss("gaussian",qbounds,seed,Lx,Ly,params);
	make_Hamiltonian(p_gauss);
//	cout<<Hnn<<endl;
//	cout<<proj.eigenvectors()<<endl;
//	cout<<proj.eigenvalues()<<endl;
}
void SingleSolver::printEnergy(int nLow, int nHigh){
	N_trunc=NPhi-nLow-nHigh;
	Eigen::MatrixXcd Hmm=Hnn_to_Hmm(proj, nLow);
//	cout<<Hmm<<endl;
	proj.compute(Hmm);
	cout<<proj.eigenvalues()(0)+proj.eigenvalues()(1)+proj.eigenvalues()(2)+proj.eigenvalues()(3)<<endl;
//	cout<<proj.eigenvalues()(0)+proj.eigenvalues()(3)<<endl;
	cout<<proj.eigenvalues()(0)+proj.eigenvalues()(2)<<endl;
	cout<<proj.eigenvalues()(1)+proj.eigenvalues()(2)<<endl;
}
complex<double> SingleSolver::getH(int a,int b){
	return Vstrength*Hnn(a,b);
}
void SingleSolver::switchVstrength(){ Vstrength=-Vstrength; }
complex<double> SingleSolver::evec(int a,int b){
	return proj.eigenvectors().col(a)(b);
}
void SingleSolver::visualizer(){
	map <string, double> params;
	params["Vstrength"]=1.;
	params["ndeltas"]=1.;
	Potential p_delta("delta",qbounds,0,Lx,Ly,params);
	p_delta.plot_potential();
	make_Hamiltonian(p_delta);
//	cout<<Hnn<<endl;
	Eigen::SelfAdjointEigenSolver<MatrixXcd> es(Hnn);
//	plot_density(es);
	density_near_x(es,p_delta.xloc[0],p_delta.yloc[0]);
//	cout<<es.eigenvalues()<<endl;
}
void SingleSolver::run(){

	Eigen::VectorXd full_sigma=Eigen::VectorXd::Zero(energy_mesh);//averaged of realizations and/or theta
	Eigen::VectorXd full_thouless=Eigen::VectorXd::Zero(energy_mesh);

	Eigen::MatrixXcd Hmm;
	Eigen::SelfAdjointEigenSolver<MatrixXcd> es_n(NPhi);
	Eigen::SelfAdjointEigenSolver<MatrixXcd> es_m(N_trunc);
	Eigen::VectorXd evals(N_trunc);
	map <string, double> params;
	params["ndeltas"]=2;//(NPhi-N_trunc)/2;

	for(int i=0;i<NROD;i++){
		params["Vstrength"]=Vstrength;
		
		Potential p_delta("delta",qbounds,i+shift*NROD,Lx,Ly,params);//make delta function potential
//		if (N_trunc<NPhi){
			make_Hamiltonian(p_delta);//make delta function hamiltonian
			es_n.compute(Hnn);//diagonalize hamiltonian
//		}
		for(int sign=0;sign<2;sign++){//zero for positive realization, 1 for negative
			if(sign==1) params["Vstrength"]=-Vstrength;
			Potential p_gauss("gaussian",qbounds,i,Lx,Ly,params);
			make_Hamiltonian(p_gauss);
//			if(N_trunc<NPhi) Hmm=Hnn_to_Hmm(es_n);
//			else Hmm=Hnn;
			cout<<es_n.eigenvectors()<<endl;
			Hmm=Hnn_to_Hmm(es_n,(NPhi-N_trunc)/2);

			es_m.compute(Hmm);

			//calculate Hall conductivity
			full_sigma-=conductivity(es_m);
			evals=es_m.eigenvalues();
			full_thouless+=ThoulessNumber(p_delta,p_gauss,evals);
		}//disorder sign loop
	}//ROD loop
	print(full_sigma,full_thouless);
}

//calculate conductivity for a realization of disorder
Eigen::VectorXd SingleSolver::conductivity(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es_m){
	Eigen::VectorXd sigma=Eigen::VectorXd::Zero(energy_mesh);
	complex<double> a1,a2;
	double temp=0;
	Eigen::VectorXd evals=es_m.eigenvalues();

	//calculate conductivities for each n,m ahead of time so we don't have to repeat these calculations
	Eigen::MatrixXd s_nm(N_trunc,N_trunc);
	for (int n=0;n<N_trunc;n++){
		for(int m=N_trunc-1;m>=0;m--){
			if (m==n) continue;
			a1=es_m.eigenvectors().col(m).adjoint()*dVdX*es_m.eigenvectors().col(n);
			a2=es_m.eigenvectors().col(n).adjoint()*dVdY*es_m.eigenvectors().col(m);
			s_nm(n,m)=-2*(a1.real()*a2.imag()+a1.imag()*a2.real())/pow(evals(m)-evals(n),2);
			
		}
	}

	for(int e=0;e<energy_mesh;e++){//compute conductivity only between states where one is less than the fermi energy and one is greater
		for (int n=0;n<N_trunc;n++){
			if(evals(n)>energies(e)) break;
			for(int m=N_trunc-1;m>=0;m--){
				if (m==n) continue;
				if(evals(m)<energies(e)) break;
				temp+=s_nm(n,m);			
			}
		sigma(e)+=temp-1;
		temp=0;
		}
	}//energy loop
//	cout<<sigma<<endl<<"-----------------"<<endl;
	return sigma;
}
void SingleSolver::print(Eigen::VectorXd &full_sigma,Eigen::VectorXd &full_thouless){
	full_sigma/=(2.*N_trunc*NROD);
	ofstream condout;
	condout.open("out");	
	for(int i=0;i<energy_mesh;i++)
		condout<<i<<" "<<energies(i)<<" "<<full_sigma(i)<<" "<<full_thouless(i)<<" "<<spacings(i)<<endl;
	condout.close();


}
Eigen::VectorXd SingleSolver::ThoulessNumber(Potential &deltas, Potential &gaussian, Eigen::VectorXd &oldEnergies){
	Eigen::MatrixXcd Hmm; //Hnn;
	Eigen::SelfAdjointEigenSolver<MatrixXcd> es_n(NPhi);
	Eigen::SelfAdjointEigenSolver<MatrixXcd> es_m(N_trunc);
	Eigen::VectorXd thouless(energy_mesh);
		
	if (N_trunc<NPhi){
		make_Hamiltonian(deltas,M_PI,0.);//make delta function hamiltonian
		//cout<<"Hnn print"<<endl<<Hnn<<endl;
		es_n.compute(Hnn);//diagonalize hamiltonian
	}
	make_Hamiltonian(gaussian,M_PI,0.);
	if(N_trunc<NPhi) Hmm=Hnn_to_Hmm(es_n,(NPhi-N_trunc)/2);
	else Hmm=Hnn;

	es_m.compute(Hmm);
	Eigen::VectorXd newEnergies=es_m.eigenvalues();
	for(int e=0;e<energy_mesh;e++){
		int n1=0;
		while(n1<N_trunc){
			if(oldEnergies(n1)>energies(e)) break;
			n1++;
		}
		if(n1==0 ||n1==N_trunc){
			thouless(e)=0;
			continue;
		}
//		while(n2<N_trunc){
//			if(newEnergies(n2)>energies(e)) break;
//			n2++;
//		}
//		if(n2==0 ||n2==N_trunc){
//			thouless(e)=0;
//			continue;
//		}
		thouless(e)=pow(newEnergies(n1)-oldEnergies(n1),2)+pow(newEnergies(n1-1)-oldEnergies(n1-1),2);
		spacings(e)+=pow(newEnergies(n1)-newEnergies(n1-1),2)+pow(oldEnergies(n1)-oldEnergies(n1-1),2);
	}
	return thouless;		

}
///make Hamiltonian from Vqs
//right now this only makes an lower-triangular Hamiltonian
//also computes derivatives of the potential for use when calculating Hall conductivities
void SingleSolver::make_Hamiltonian(Potential &pot,double theta_x,double theta_y){
//	Eigen::MatrixXcd Hnn;
	Hnn.resize(NPhi,NPhi); dVdX.resize(NPhi,NPhi); dVdY.resize(NPhi,NPhi);
	double qx,qy,arg, p_arg, exp_qx,exp_qy;
	complex<double> tempH,tempX,tempY,temp;
	double tol=1e-16;
	int mx;
	for(int n1=0;n1<NPhi;n1++){
		for(int n2=n1;n2<NPhi;n2++){
			for(int psign=0;psign<2;psign++){
				for(int p=psign;p>-1;p++){
					if(psign==0){
						mx=n2-n1+p*NPhi;	
						p_arg=p*theta_y;
					}else{
						mx=n2-n1-p*NPhi;
						p_arg=-p*theta_y;
					}
					qx=2*M_PI*mx/(1.*Lx);
					exp_qx=exp(-0.25*qx*qx);
					if(exp_qx<tol) break;

					for(int my=0;my>-1;my++){
						qy=2*M_PI*my/(1.*Ly);
						exp_qy=exp(-0.25*qy*qy);
						if(exp_qx*exp_qy<tol) break;
						if(my!=0){
							arg=M_PI*(1.*my*(n1+n2))/(1.*NPhi) +qy*theta_x/Lx+M_PI*my*p;
							for(int msign=-1;msign<2;msign+=2){
								temp=exp_qx * exp_qy * (pot.get_potential(mx,msign*my) * polar(1.0,msign*arg-p_arg));
								tempH+=temp;
								tempX+=qx*complex<double>(-temp.imag(),temp.real());
								tempY+=msign*qy*complex<double>(-temp.imag(),temp.real());
							}
						}else{
							temp=exp_qx*exp_qy*(pot.get_potential(mx,my)*polar(1.0,-p_arg));
							//temp/=(2.*M_PI*NPhi);
							tempH+=temp;
							tempX+=qx*complex<double>(-temp.imag(),temp.real());
						}
					}
				}
			}
			//cout<<"-----------------------"<<endl;
			Hnn(n2,n1)=tempH;
			Hnn(n1,n2)=conj(tempH);
			dVdX(n2,n1)=tempX; dVdX(n1,n2)=conj(tempX);
			dVdY(n2,n1)=tempY; dVdY(n1,n2)=conj(tempY);
			tempH=complex<double>(0,0); tempX=complex<double>(0,0); tempY=complex<double>(0,0);
		}
	}
//	return Hnn;
}
Eigen::MatrixXcd SingleSolver::Hnn_to_Hmm(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, int nLow){
	Eigen::MatrixXcd Hmm, tempX,tempY;
	Hmm.resize(N_trunc,N_trunc);
	tempX.resize(N_trunc,N_trunc);
	tempY.resize(N_trunc,N_trunc);
	Eigen::MatrixXcd rect=es.eigenvectors().middleCols(nLow,N_trunc);
	Hmm=rect.transpose().conjugate()*Hnn*rect;
	tempX=rect.transpose().conjugate()*dVdX*rect;
	tempY=rect.transpose().conjugate()*dVdY*rect;
	dVdX.resize(N_trunc,N_trunc);
	dVdY.resize(N_trunc,N_trunc);
	dVdX=tempX; dVdY=tempY;
	return Hmm;
}
void SingleSolver::plot_density(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, double theta_x,double theta_y){
	int gridsize=100;
	double x,y,arg;
	complex<double> temp;
	stringstream filename;
	ofstream overlapout;
	//initialize a vector of vectors of 2x2 matrices (I wish I had iTensor)
	MatrixXcd temp1(gridsize,gridsize);
	vector<MatrixXcd> temp2(NPhi,temp1);
	vector<vector<MatrixXcd> > overlaps(NPhi,temp2);

	for(int n1=0;n1<NPhi;n1++){
		for(int n2=n1;n2<NPhi;n2++){
			for(int ix=0;ix<gridsize;ix++){
				x=Lx*ix/(1.*gridsize);
				for(int iy=0;iy<gridsize;iy++){
					y=Ly*iy/(1.*gridsize);
					for(int p1=-CUTOFF;p1<=CUTOFF;p1++){
						for(int p2=-CUTOFF;p2<=CUTOFF;p2++){
							arg=(p2-p1)*theta_y+2.*M_PI*x/Lx*(n2-n1+NPhi*(p2-p1));
							temp+=exp(-0.5*pow( y-(2.*M_PI*(n2+p2*NPhi)+theta_x)/Lx,2) )*exp(-0.5*pow( y-(2.*M_PI*(n1+p1*NPhi)+theta_x)/Lx,2) )*polar(1.,arg);
						}
					}
					overlaps[n1][n2](ix,iy)=temp/(Lx*sqrt(M_PI));
					temp=complex<double> (0,0);	
				}
			}
//			filename.str("");
//			filename<<"visualizer/"<<n1<<"_"<<n2;
//			overlapout.open(filename.str().c_str());
//			overlapout<<overlaps[n1][n2].imag()<<endl;
//			overlapout.close();
		}
	}
	vector<MatrixXcd> probs(NPhi,temp1);
	Eigen::MatrixXcd evec=es.eigenvectors();
	for(int m=0;m<NPhi;m++){
		for(int n1=0;n1<NPhi;n1++){
			for(int n2=n1;n2<NPhi;n2++){
				if(n1==n2) probs[m]+=conj(evec(n1,m))*evec(n2,m)*overlaps[n1][n2];
				else probs[m]+=2.*conj(evec(n1,m))*evec(n2,m)*overlaps[n1][n2];
			}
		}
		filename.str("");
		filename<<"visualizer/"<<m;
		overlapout.open(filename.str().c_str());
		overlapout<<probs[m].real()<<endl;
		overlapout.close();
	}
				
}
void SingleSolver::density_near_x(Eigen::SelfAdjointEigenSolver<MatrixXcd> &es, double xloc,double yloc,double theta_x,double theta_y){
	int gridsize=500;
//	cout<<Lx<<" "<<Ly<<" "<<xloc<<" "<<yloc<<endl;
	double x,y,arg;
	complex<double> temp;
	stringstream filename;
	ofstream overlapout;
	//initialize a vector of vectors of 2x2 matrices (I wish I had iTensor)
	Eigen::VectorXcd temp1=Eigen::VectorXcd::Zero(gridsize+1);
	vector<Eigen::VectorXcd> temp2(NPhi,temp1);
	vector<vector<Eigen::VectorXcd> > overlaps(NPhi,temp2);

	for(int n1=0;n1<NPhi;n1++){
		for(int n2=n1;n2<NPhi;n2++){
			for(int ix=0;ix<=gridsize;ix++){
				y=Ly+yloc-Ly*ix/(1.*gridsize);
				x=xloc;////location of delta functions
				temp=complex<double> (0,0);	
				for(int p1=-CUTOFF;p1<=CUTOFF;p1++){
					for(int p2=-CUTOFF;p2<=CUTOFF;p2++){
						arg=(p2-p1)*theta_y+2.*M_PI*x/Lx*(n2-n1+NPhi*(p2-p1));
						temp+=exp(-0.5*pow( y-(2.*M_PI*(n2+p2*NPhi)+theta_x)/Lx,2) )*exp(-0.5*pow( y-(2.*M_PI*(n1+p1*NPhi)+theta_x)/Lx,2) )*polar(1.,arg);
					}
				}
				overlaps[n1][n2](ix)=temp/(Lx*sqrt(M_PI));
				if(n1==0&&n2==0&&y==0) cout<<temp<<" "<<overlaps[0][0](0)<<endl;
			}
//			filename.str("");
//			filename<<"visualizer/"<<n1<<"_"<<n2;
//			overlapout.open(filename.str().c_str());
//			overlapout<<overlaps[n1][n2].imag()<<endl;
//			overlapout.close();
		}
	}
	vector<Eigen::VectorXcd> probs(NPhi,temp1);
	Eigen::MatrixXcd evec=es.eigenvectors();
	//evec=Eigen::MatrixXcd::Identity(NPhi,NPhi);
	for(int m=0;m<NPhi;m++){
		for(int n1=0;n1<NPhi;n1++){
			for(int n2=n1;n2<NPhi;n2++){
				if(n1==n2) probs[m]+=conj(evec(n2,m))*evec(n1,m)*overlaps[n1][n2];
				else probs[m]+=2.*conj(evec(n1,m))*evec(n2,m)*overlaps[n1][n2];
			}
		}
		filename.str("");
		filename<<"visualizer/"<<m<<"delta";
		overlapout.open(filename.str().c_str());
		for(int i=0;i<=gridsize;i++)
			overlapout<<Lx-Lx*i/(gridsize)<<" "<<probs[m](i).real()<<endl;
		overlapout.close();
	}
}
//template<class ART>
//SingleSolver<ART>& SingleSolver<ART>::operator=(const SingleSolver<ART>& other)
//{

//  if (this != &other) { // Stroustrup suggestion.
//    this->ClearMem();
//    Copy(other);
//  }
//  return *this;

//} // operator=.
