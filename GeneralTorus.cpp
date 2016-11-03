#include "GeneralTorus.h"

GeneralTorus::GeneralTorus(int tcharge):ManySolver<complex<double> >(){
	ifstream params("genparams");
	double alpha, theta; //alpha is |L1|/|L2|, L1.L2=|L1||L2|cos(theta), theta is in degrees
	params>>alpha;
	params>>theta;
	
	theta=theta*2*M_PI/(360.);
	
	Lx=sqrt(2*M_PI*NPhi * alpha/sin(theta));
	Ly=Lx/alpha*sin(theta);
	LDelta=Lx/alpha*cos(theta);

//	cout<<"Ls:"<<Lx<<" "<<Ly<<" "<<LDelta<<endl;
	
	tau=complex<double>(LDelta/Lx,Ly/Lx);
//	cout<<"tau: "<<tau<<endl;
	cache=1;
	store_sparse=false;
	disorder=0;
	project=0;
	init(tcharge);
	//interaction_cache();
	periodic=1;
	
}

GeneralTorus::GeneralTorus(int tcharge, double alpha, double theta):ManySolver<complex<double> >(){
	//alpha is |L1|/|L2|, L1.L2=|L1||L2|cos(theta), theta is in degrees
	theta=theta*2*M_PI/(360.);
	
	Lx=sqrt(2*M_PI*NPhi * alpha/sin(theta));
	Ly=Lx/alpha*sin(theta);
	LDelta=Lx/alpha*cos(theta);

//	cout<<"Ls:"<<Lx<<" "<<Ly<<" "<<LDelta<<endl;
	
	tau=complex<double>(LDelta/Lx,Ly/Lx);
//	cout<<"tau: "<<tau<<endl;
	cache=1;
	store_sparse=false;
	disorder=0;
	project=0;
	init(tcharge);
	//interaction_cache();
	periodic=1;
	
}

void GeneralTorus::ground_state(int tcharge){
	if(charge!=tcharge) init(tcharge);
	make_Hnn();
//	cout<<EigenDense<<endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(EigenDense);
	for(int i=0;i<10;i++) cout<<es.eigenvalues()(i)<<" "<<charge<<endl;
//	return es.eigenvectors()(0);
}
int GeneralTorus::get_charge(int s){
	int out=0;
	for(unsigned int i=0;i<(unsigned)NPhi;i++){
		if (s & 1<<i) out+=i;
	}
	return out%NPhi;
}

void GeneralTorus::make_Vmk(){
	Vmk=vector< vector <complex<double> > > (NPhi,vector<complex<double> >(NPhi));

	complex<double> bigout,qyphase;
	int startm,startq;
	double expqy,expm,tol=1e-10;

	for(int m=0;m<NPhi;m++){
		
		for(int k=0;k<NPhi;k++){
			
			bigout=0.;

			for(int pmsign=1;pmsign>-2;pmsign-=2){
				if(pmsign==1) startm=0;
				else startm=-1;
				
				for(int pm=startm;true;pm+=pmsign){
					expm=exp(-M_PI*pow( (m+pm*NPhi)*imag(tau),2)/(imag(tau)*NPhi));
					if(expm<tol) break;
					
					for(int nqsign=1;nqsign>=-1;nqsign-=2){
						if(nqsign==1) startq=0;
						else startq=-1;
						
						for(int nqy=startq;true;nqy+=nqsign){
							if(m==0 && nqy==0 && pm==0) continue;
							expqy=exp(-M_PI*pow(nqy+real(tau)*(m+pm*NPhi),2)/(imag(tau)*NPhi));
							if (expqy<tol) break;
							qyphase=polar(1.,2*M_PI*k*nqy/(1.*NPhi));
				
							bigout+=expqy*qyphase*expm*V_Coulomb(2*M_PI*(m+pm*NPhi)/Lx,2*M_PI*(nqy+real(tau)*(m+pm*NPhi))/Ly);
							//here I for some reason attempted pseudopotentials
							//bigout+=expqy*qyphase*expm*( (k*k-(m+pm*NPhi)*(m+pm*NPhi))*pow(imag(tau)*2*M_PI/Lx,2) );


						}
					}
				}
			}
			Vmk[m][k]=Ly*Ly/(2.*M_PI*imag(tau)*NPhi)*bigout;
//			cout<<m<<" "<<k<<" "<<Vmk[m][k]<<" "<<Ly*Ly/(2.*M_PI*imag(tau)*NPhi)<<" "<<expm<<endl;
		}//k
	}//m
}
//this function replaces the one in manysolver with a better version which only works to calculate the Vmk, and then copies them into many positions
//the relation between a,b,c,d is taken from Mong, Zalatel, Pollmann
void GeneralTorus::interaction_cache(){
	int d,m,k;
	complex<double> temp;
	four_body_cache=vector<complex<double> > (pow(NPhi,4),complex<double>(0,0));
	make_Vmk();
	
	for(int a=0;a<NPhi;a++){
		for(int b=a+1;b<NPhi;b++){

			for(int c=0;c<NPhi;c++){
				d=supermod(a+b-c,NPhi);
				if (d>=c) continue;
				k=supermod(c-a,NPhi);
				m=supermod(d-a,NPhi);
 				//cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<m<<" "<<k<<" "<<Vmk[m][k]<<endl;
				four_body_cache[four_array_map(a,b,c,d)]+=Vmk[m][k];
				k=supermod(d-a,NPhi);
				m=supermod(c-a,NPhi);
 				//cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<m<<" "<<k<<" "<<Vmk[m][k]<<endl<<endl;
				four_body_cache[four_array_map(a,b,c,d)]-=Vmk[m][k];
			}	
		}
	}
}
complex<double> GeneralTorus::two_body(int a, int b){
	//actually don't need this function, only in here because I made twobody pure virtual in ManySolver
	cout<<"you called two body, but you shouldn't!"<<endl;
	return 0.;
}
complex<double> GeneralTorus::four_body(int a, int b, int c, int d){
	//actually don't need this function, only in here because I made twobody pure virtual in ManySolver
	cout<<"you called two body, but you shouldn't!"<<endl;
	return 0.;
}

//make shrinking matrix
void GeneralTorus::makeShrinker(int nx){
	vector<Eigen::Triplet< complex<double> > > triplets;
	vector<Eigen::Triplet< complex<double> > > temptrips;
	vector<int> found_states; 

	int temp,index;
	int invNu=NPhi/Ne;
	vector<int>::iterator it;
	int col=0, phase,sign;

//	for(int i=0;i<nStates;i++) cout<<(bitset<12>)states[i]<<endl;
//	cout<<endl;
	for(int i=0;i<nStates;i++){
		if(find(found_states.begin(),found_states.end(),i)!=found_states.end()) continue;
		
		it=states.begin()+i;
		temp=states[i];
		phase=0;
		temptrips.clear();
		sign=1;
		while(true){
			index=it-states.begin();
			//hopefully the only thing we need to do for a non-rectangular torus is add some constants to this line
			temptrips.push_back( Eigen::Triplet< complex<double> >( col,index,polar(1.*sign,phase*2*M_PI/(1.*Ne)*nx ) ) );
			found_states.push_back(index);
//			cout<<(bitset<12>)states[index]<<" "<<phase*nx%Ne<<" "<<sign<<endl;
			temp=cycle_M(temp,NPhi,invNu,sign);
			if(temp==states[i]) break;
			it=lower_bound(states.begin(),states.end(),temp);
			phase++;
			
		}
		//some cases only exist in certain momentum sectors (eg 0101) 
		if( (Ne%2!=0 && nx%(Ne/temptrips.size()) ) || (Ne%2==0 && (nx-Ne/2)%(Ne/temptrips.size()) ) ) {
//			cout<<"cancelled"<<endl;
			 continue;
		}

		for(unsigned int j=0;j<temptrips.size();j++) 
			temptrips[j]=Eigen::Triplet< complex<double> >(temptrips[j].col(), temptrips[j].row(), temptrips[j].value()/sqrt(temptrips.size()));
		triplets.insert(triplets.end(),temptrips.begin(),temptrips.end());
//		cout<<endl;
		col++;
	}
	shrinkMatrix=Eigen::SparseMatrix<complex<double> >(nStates,col);
//	cout<<"triplet size "<<triplets.size()<<endl;
//	for(unsigned int j=0;j<triplets.size();j++)
//		cout<<triplets[j].row()<<" "<<triplets[j].col()<<" "<<triplets[j].value()<<" "<<(bitset<10>)states[triplets[j].row()]<<endl;
	shrinkMatrix.setFromTriplets(triplets.begin(),triplets.end());
	shrinkMatrix=shrinkMatrix.adjoint();
}
//construct denstiy operator rho(kx,ky)
Eigen::SparseMatrix< complex<double> > GeneralTorus::density_operator(int mx, int my){
	//this operator doesn't conserve charge, and so it maps between different bases
	//this part constructs the basis to map to, it works just like make_states
	//the new basis has its charge INCREASED by kx
	int verbose=0;
	vector<int> new_states;
	if(mx==0) new_states=states;
	else{
		for(int i=0;i<intpow(2,NPhi);i++){
			if(count_bits(i)==Ne && get_charge(i)==((charge+mx)%NPhi) )
				new_states.push_back(i);
		}
	}
	int newNStates=new_states.size();

	if(verbose>0){	
		for(int i=0;i<(signed)states.size();i++) cout<<(bitset<16>)states[i]<<endl;
		cout<<endl;
		for(int i=0;i<(signed)new_states.size();i++) cout<<(bitset<16>)new_states[i]<<endl;
	}
	Eigen::SparseMatrix<complex<double> > rho(newNStates,nStates);
	vector<Eigen::Triplet<complex<double> > > triplets;

	//for every state, find all the states which you get from moving one electron kx, and add those matrix elements to rho
	complex<double>  prefactor;
	double sign,ky;
	int lmx;
	if(mx>NPhi/2) lmx=(mx-NPhi);
	else lmx=mx;
	if(my>NPhi/2) ky=2*M_PI*(my-NPhi);
	else ky=2*M_PI*my;
	vector<int>::iterator it;
	int newstate;
	prefactor=Ly/sqrt(2*M_PI*imag(tau)*NPhi);
	prefactor*=exp(-pow(ky,2)/(8*M_PI*imag(tau)*NPhi)); //qy^2
	prefactor*=exp(-M_PI*norm(tau)*pow(lmx,2)/(2.*imag(tau)*NPhi)); //qx^2
	prefactor*=polar(exp(-lmx*ky*real(tau)/(2.*imag(tau)*NPhi)), lmx*ky/(2.*NPhi)); //qxqy

	for(int i1=0; i1<nStates; i1++){
		if(verbose>1) cout<<(bitset<6>)states[i1]<<endl;
		for(int x=0; x<NPhi; x++){
			try{
				newstate=move_bit(states[i1],NPhi,x,mx);
			}catch(int e){
				continue;
			}
			if(verbose>1) cout<<x<<" "<<(bitset<6>)newstate<<endl;
			it=find(new_states.begin(),new_states.end(),newstate);
			if(it==new_states.end()) continue;
			//a minus sign from normal ordering if there's an even number of electrons
			if( newstate<states[i1] && Ne%2==0 ) sign=-1.;
			else sign=1.;
			//another minus sign for every electron this electron hops over
			if (x+mx>NPhi){
				for(int y=0;y<(x+mx)%NPhi;y++)
					if(states[i1] & 1<<y) sign*=-1;
				for(int y=x+1;y<NPhi;y++)
					if(states[i1] & 1<<y) sign*=-1;
			}else{
				for(int y=x+1;y<x+mx;y++) 
					if(states[i1] & 1<<y) sign*=-1;
			}
			triplets.push_back(Eigen::Triplet<complex<double> >( it-new_states.begin(),i1,sign*prefactor*polar(1.,ky*x/(1.*NPhi) ) ) );
		}
	}
	rho.setFromTriplets(triplets.begin(),triplets.end());
	if(verbose>1) cout<<rho<<endl;
	return rho;				
}
