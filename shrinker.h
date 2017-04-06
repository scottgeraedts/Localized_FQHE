#include "TorusSolver.h"

//unsigned int cycle_M(unsigned int in, int NPhi, int M, int &sign){
//	int out,old=in;
//	for(int i=0;i<M;i++){
//		out=cycle_bits(old,NPhi);
//		if (out<old && NPhi%4==0) sign*=-1;
//		old=out;
//	}
//	return out;
//}

//make shrinking matrix
template<class ART>
void TorusSolver<ART>::makeShrinker(int nx){
	vector<Eigen::Triplet< complex<double> > > triplets;
	vector<Eigen::Triplet< complex<double> > > temptrips;
	vector<int> found_states;

	int temp,index;
	int invNu=this->NPhi/this->Ne;
	vector<int>::iterator it;
	int col=0, phase,sign;

//	for(int i=0;i<this->nStates;i++) cout<<(bitset<12>)this->states[i]<<endl;
//	cout<<endl;
	for(int i=0;i<this->nStates;i++){
		if(find(found_states.begin(),found_states.end(),i)!=found_states.end()) continue;
		
		it=this->states.begin()+i;
		temp=this->states[i];
		phase=0;
		temptrips.clear();
		sign=1;
		while(true){
			//shift by invnu, adding an element each time, until you get back to the start
			index=it-this->states.begin();
			temptrips.push_back( Eigen::Triplet< complex<double> >( col,index,polar(1.*sign,phase*2*M_PI/(1.*this->Ne)*nx ) ) );
			found_states.push_back(index);
//			cout<<(bitset<12>)this->states[index]<<" "<<phase*nx%this->Ne<<" "<<sign<<endl;
			temp=cycle_M(temp,this->NPhi,invNu,sign);
			if(temp==this->states[i]) break;
			it=lower_bound(this->states.begin(),this->states.end(),temp);
			phase++;
			
		}
		//some cases only exist in certain momentum sectors (eg 0101) 
		if( (this->Ne%2!=0 && nx%(this->Ne/temptrips.size()) ) || (this->Ne%2==0 && (nx-this->Ne/2)%(this->Ne/temptrips.size()) ) ) {
//			cout<<"cancelled"<<endl;
			 continue;
		}

		for(unsigned int j=0;j<temptrips.size();j++) 
			temptrips[j]=Eigen::Triplet< complex<double> >(temptrips[j].col(), temptrips[j].row(), temptrips[j].value()/sqrt(temptrips.size()));
		triplets.insert(triplets.end(),temptrips.begin(),temptrips.end());
//		cout<<endl;
		col++;
	}
	shrinkMatrix=Eigen::SparseMatrix<ART>(this->nStates,col);
//	cout<<"triplet size "<<triplets.size()<<endl;
//	for(unsigned int j=0;j<triplets.size();j++)
//		cout<<triplets[j].row()<<" "<<triplets[j].col()<<" "<<triplets[j].value()<<" "<<(bitset<10>)this->states[triplets[j].row()]<<endl;
	shrinkMatrix.setFromTriplets(triplets.begin(),triplets.end());
	shrinkMatrix=shrinkMatrix.adjoint();
}
//construct denstiy operator rho(kx,ky)
template<class ART>
Eigen::SparseMatrix<ART> TorusSolver<ART>::density_operator(int mx, int my){
	//this operator doesn't conserve charge, and so it maps between different bases
	//this part constructs the basis to map to, it works just like make_states
	//the new basis has its charge INCREASED by kx
	int verbose=0;
	vector<int> new_states;
	if(mx==0) new_states=this->states;
	else{
		for(int i=0;i<intpow(2,this->NPhi);i++){
			if(count_bits(i)==this->Ne && get_charge(i)==((this->charge+mx)%this->NPhi) )
				new_states.push_back(i);
		}
	}
	int newNStates=new_states.size();

	if(verbose>0){	
		for(int i=0;i<this->states.size();i++) cout<<(bitset<16>)this->states[i]<<endl;
		cout<<endl;
		for(int i=0;i<new_states.size();i++) cout<<(bitset<16>)new_states[i]<<endl;
	}
	Eigen::SparseMatrix<ART> rho(newNStates,this->nStates);
	vector<Eigen::Triplet<ART> > triplets;

	//for every state, find all the states which you get from moving one electron kx, and add those matrix elements to rho
	ART prefactor;
	double sign,kx,ky;
	if(mx>this->NPhi/2) kx=2*M_PI/Lx*(mx-this->NPhi);
	else kx=2*M_PI/Lx*mx;
	if(my>this->NPhi/2) ky=2*M_PI/Ly*(my-this->NPhi);
	else ky=2*M_PI/Ly*my;
	vector<int>::iterator it;
	int newstate;
	prefactor=polar(exp(-0.25*(pow(kx,2)+pow(ky,2))), 0.5*kx*ky);
	for(int i1=0; i1<this->nStates; i1++){
		if(verbose>1) cout<<(bitset<6>)this->states[i1]<<endl;
		for(int x=0; x<this->NPhi; x++){
			try{
				newstate=move_bit(this->states[i1],this->NPhi,x,mx);
			}catch(int e){
				continue;
			}
			if(verbose>1) cout<<x<<" "<<(bitset<6>)newstate<<endl;
			it=find(new_states.begin(),new_states.end(),newstate);
			if(it==new_states.end()) continue;
			//a minus sign from normal ordering if there's an even number of electrons
			if( newstate<this->states[i1] && this->Ne%2==0 ) sign=-1.;
			else sign=1.;
			//another minus sign for every electron this electron hops over
			if (x+mx>this->NPhi){
				for(int y=0;y<(x+mx)%this->NPhi;y++)
					if(this->states[i1] & 1<<y) sign*=-1;
				for(int y=x+1;y<this->NPhi;y++)
					if(this->states[i1] & 1<<y) sign*=-1;
			}else{
				for(int y=x+1;y<x+mx;y++) 
					if(this->states[i1] & 1<<y) sign*=-1;
			}
			triplets.push_back(Eigen::Triplet<ART>( it-new_states.begin(),i1,sign*prefactor*polar(1.,ky*2*M_PI/Lx*x) ) );
		}
	}
	rho.setFromTriplets(triplets.begin(),triplets.end());
	if(verbose>1) cout<<rho<<endl;
	cout<<rho<<endl;
	return rho;				
}
