#include "SingleSolver.h"
Potential::Potential(string _type, int _qbounds, int seed, double _Lx, double _Ly, map <string, double> &params):qbounds(_qbounds),type(_type),Lx(_Lx),Ly(_Ly){
	ran.seed(seed);
	V.resize(2*qbounds+1,qbounds+1);

	if (type=="gaussian") make_potential_gaussian(params);
	else if(type=="delta") make_potential_delta(params);
	else cout<<"type not recognized"<<endl;
}
void Potential::make_potential_gaussian(map <string, double> &params){
	double Vstrength=params["Vstrength"];
	Eigen::MatrixXcd temp_Vq;
	temp_Vq.resize(2*qbounds+1,qbounds+1);
	for(int i=0;i<2*qbounds+1;i++){
		for(int j=0;j<=qbounds;j++){
			if (i==qbounds && j==0) V(i,j)=complex<double>(ran.randNorm(0,Vstrength),0);
			else if(j==0 && i>qbounds) V(i,j)=conj(V(2*qbounds-i,j));
			else V(i,j)=complex<float>(ran.randNorm(0,Vstrength),ran.randNorm(0,Vstrength));
		}
	}
}
void Potential::make_potential_delta(map <string, double> &params){
	double qx,qy;
	complex<double> temp(0,0);
	int ndeltas=floor(params["ndeltas"]+0.1);
	double Vstrength=params["Vstrength"];
	for(int i=0;i<ndeltas;i++){
		xloc.push_back(ran.rand(Lx)); yloc.push_back(ran.rand(Ly));sign.push_back(1);//list of locations and signs of delta functions
		xloc.push_back(ran.rand(Lx)); yloc.push_back(ran.rand(Ly)); sign.push_back(-1);
	}
//			xloc.push_back(ran.rand(0)); yloc.push_back(ran.rand(0));sign.push_back(1);//list of locations and signs of delta functions
	for(int i=0;i<2*qbounds+1;i++){
		for(int j=0;j<=qbounds;j++){
			qx=2.*M_PI*(i-qbounds)/(1.*Lx); qy=2.*M_PI*j/(1.*Ly);
			temp=complex<double>(0,0);
			for(unsigned int k=0;k<xloc.size();k++)
				temp+=sign[k]*polar(Vstrength,-xloc[k]*qx-yloc[k]*qy);
			V(i,j)=temp;	
		}
	}
}
complex<double> Potential::get_potential(int mx,int my){
	//note that Vgauss only stores elements for qy>=0, if qy<0 you must take negative and complex conjugate
	//also in the x direction the index 0 stores the element mx=-CUTOFF
	if(mx<-qbounds || mx>qbounds || my<-qbounds || my>qbounds){
		cout<<"error getting potential"<<endl;
		exit(0);
	}
	if (my>=0) return V(mx+qbounds,my);
	else return conj(V(-mx+qbounds,-my));
}
int Potential::get_qbounds(){ return qbounds;}
void Potential::plot_potential(){
	if (type!="delta") cout<<"wrong type for plot potential"<<endl;
	ofstream outfile;
	outfile.open("visualizer/potential");
	outfile<<Lx<<" "<<Ly<<" "<<floor(Lx*Lx/2./M_PI+0.1)<<endl;
	for (unsigned int i=0;i<xloc.size();i+=2) outfile<<xloc[i]<<" "<<yloc[i]<<" "<<sign[i]<<endl;
	for (unsigned int i=1;i<xloc.size();i+=2) outfile<<xloc[i]<<" "<<yloc[i]<<" "<<sign[i]<<endl;
	outfile.close();
}
