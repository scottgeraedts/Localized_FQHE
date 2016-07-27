#include "matprod2.h"

void MatrixWithProduct2::set_mode(string mode){ LU_mode=mode; }

void MatrixWithProduct2::CSR_from_Dense(Eigen::MatrixXcd EigenDense){

	oldE=0.;

	clock_t CPUtime=clock();
	time_t walltime=time(NULL);	
	if(first_time){
		nonzero=0;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++)
				if(abs(EigenDense(i,j))>1e-12) nonzero++;
		}
	   	if(verbose>0) cout<<"nonzero entries in eigendense: "<<nonzero<<endl;
		ia=new int[n+1];
		ja=new int[nonzero];
	 	vals=new complex<double>[nonzero];
		if(LU_mode=="arpack++") raw_evals=new complex<double>[n];
	}
	complex<double> *dense=new complex<double>[n*n];
	Eigen::Map<Eigen::MatrixXcd>(dense,n,n)=EigenDense;

	int sym_flag;
	if(LU_mode=="pardiso" || LU_mode=="none") sym_flag=1;
	else sym_flag=2;
	int job[6]={0,0,0,sym_flag,nonzero,1};

	int info=0;
	mkl_zdnscsr_(job,&n,&n,dense,&n,vals,ja,ia,&info);
	nonzero=ia[n];
	delete [] dense;

	walltime=time(NULL)-walltime;
	CPUtime=clock()-CPUtime;
	if(verbose>0) cout<<"copying the sparse matrix to dense took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;
	if(first_time){
		if(LU_mode=="pardiso" || LU_mode=="superLU") analyze_pattern();
		first_time=false;
	}
}
void MatrixWithProduct2::CSR_from_Sparse(Eigen::SparseMatrix<complex<double> > EigenSparse){
	nonzero=EigenSparse.nonZeros();
	if(first_time){
		vals=new complex<double>[nonzero];
		ia=new int[n+1];
		ja=new int[nonzero];
	}
	for(int i=0;i<nonzero;i++){
		vals[i]=*(EigenSparse.valuePtr()+i);
		ja[i]=*(EigenSparse.innerIndexPtr()+i);
	}
	for(int i=0;i<n+1;i++)
		ia[i]=*(EigenSparse.outerIndexPtr()+i);
	if(first_time){
		if(LU_mode=="pardiso" || LU_mode=="superLU") analyze_pattern();
		first_time=false;
	}
	
}
//does the analysis of the shape of the matrix, to speed up future solutions
void MatrixWithProduct2::analyze_pattern(){
//call to pardiso
	if(LU_mode=="pardiso"){	
		for(int i=0;i<64;i++) pt[i]=0.;

		maxfct=1;
		mnum=1;
		mtype=-4;
		nrhs=1;

		pardisoinit_(pt,&mtype,iparm);

		iparm[0]=1;//use custom parameters
		iparm[1]=2;
		iparm[2]=1;
		iparm[3]=0;
		iparm[4]=0;
		iparm[5]=0;
		iparm[6]=0;
		iparm[7]=2;//number of iterations
		iparm[8]=0;
		iparm[9]=8;//perturbation of pivots
		iparm[10]=0;
		iparm[11]=0;
		iparm[12]=0;	
		iparm[13]=0;	
		iparm[14]=0;	
		iparm[15]=0;	
		iparm[16]=0;	
		iparm[17]=-1;	
		iparm[18]=-1;	
		iparm[19]=0;
		iparm[20]=1;
		iparm[26]=0; //error checking
		iparm[34]=1; //arrays start at 0

		msglvl=0;
		int error=0;
		int phase=11;
		pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, vals, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if(error) cout<<"error analyzing with pardiso! "<<error<<endl;
	}else if (LU_mode=="superLU"){
		panel_size = sp_ienv(1);
		relax = sp_ienv(2);
		StatInit(&stat);
		set_default_options(&options);
		options.ColPerm=MMD_ATA;
		perm_c=new int[n];
		perm_r=new int[n];
		zCreate_CompCol_Matrix(&A,n,n,nonzero,(ldcomplex*)vals,ja,ia,SLU_NC,SLU_Z,SLU_GE);
		get_perm_c(options.ColPerm, &A, perm_c);
	}	
}

void MatrixWithProduct2::MultInv(complex<double> *v, complex<double> *w){
	int error=0;
	if(LU_mode=="pardiso"){
		int phase=33;
		pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, vals, ia, ja, &idum, &nrhs, iparm, &msglvl, v, w, &error);
		if(error) cout<<"problem in pardiso solution! "<<error<<endl;
	}else if(LU_mode=="superLU"){
		trans_t trans = NOTRANS;
		SuperMatrix B;
		zCreate_Dense_Matrix(&B,n,1,(ldcomplex*)v,n,SLU_DN,SLU_Z,SLU_GE);

//		zPrint_Dense_Matrix("B",&B);
//		zPrint_SuperNode_Matrix("L",&L);
		cout<<"call to MultInv"<<L.nrow<<" "<<L.ncol<<" "<<L.Stype<<" "<<L.Dtype<<" "<<L.Mtype<<endl;
		zgstrs(trans, &L, &U, perm_c, perm_r, &B, &stat, &error); 
		Destroy_SuperMatrix_Store(&B);
		for(int i=0;i<n;i++) w[i]=v[i];
	}
}

void MatrixWithProduct2::MultMv(complex<double> *v, complex<double> *w){
	char uplo='u';
	mkl_cspblas_zcsrsymv_(&uplo, &n, vals, ia, ja, v, w);
}

void MatrixWithProduct2::print(){
	int row=0, col=0, count=0;
	while(row<n){
		if(col==ja[count]){
			cout<<vals[count]<<" ";	
			count++;
			col++;
		}else{
			cout<<"0        ";
			col++;
		}
		if(col==n){	
			row++;
			col=0;
			cout<<endl;
		}
	}
	cout<<endl;
	for(int i=0;i<nonzero;i++) cout<<ja[i]<<" ";
	cout<<endl;
	for(int i=0;i<n+1;i++) cout<<ia[i]<<" ";
	cout<<endl;

	if(LU_mode=="pardiso") for(int i=0;i<64;i++) cout<<pt[i]<<"\t\t"<<iparm[i]<<endl;
}
void MatrixWithProduct2::ShiftInvert(double E){
//subtract a constant from the diagonal elements
	int row=0;
	for(int i=0;i<nonzero;i++){
		if(ja[i]==row) vals[i]-=(E-oldE);
		if(ia[row+1]==i+1) row++;
	}
	oldE=E;

	int error=0;
	if(LU_mode=="pardiso"){
		int phase=12;
		pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, vals, ia, ja, 
			&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if(error) cout<<"error factorizing with pardiso! "<<error<<endl;
	}else if(LU_mode=="superLU"){
		//free the memory in L and U so we can reuse it
		int * etree=new int[n];
		int lwork=0;
		double drop_tol=0.;
		SuperMatrix AC;
		sp_preorder(&options, &A, perm_c, etree, &AC);
		zgstrf(&options, &AC, drop_tol, relax, panel_size, etree, 
			NULL, lwork, perm_c, perm_r, &L, &U, &stat, &error);
		delete [] etree;
		Destroy_CompCol_Permuted(&AC);
	}
}
double MatrixWithProduct2::single_energy(string type){
	ARCompStdEig<double, MatrixWithProduct2 >  dprob(ncols(), 5, this, &MatrixWithProduct2::MultMv,type,(int)0, 1e-10,1e6);
	dprob.FindEigenvalues();
	return real(dprob.Eigenvalue(4));
}

int MatrixWithProduct2::eigenvalues(int stop, double E){
	vector<double>temp(n,0);
	int Nconverged;
	if (E==-100){
		if(verbose>0) cout<<"not doing shift invert"<<endl;
		ARCompStdEig<double, MatrixWithProduct2>  dprob(ncols(), stop, this, &MatrixWithProduct2::MultMv,"SR",(int)0, 1e-10,1e6);
		dprob.FindEigenvalues();
		dprob.FindEigenvectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),vector<complex<double> >(n,0));
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=real(dprob.Eigenvalue(k));
	//		eigvecs[k]=*(dprob.StlEigenvector(k));
		}
	}else{
		if(LU_mode=="none") cout<<"error! you didn't set the LU mode!"<<endl;
		time_t walltime=time(NULL);
		clock_t CPUtime=clock();
		if(verbose>0) cout<<"making sparse"<<endl;
		if(LU_mode=="pardiso" || LU_mode=="superLU") ShiftInvert(E);
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the LU decomposition took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;

		complex<double> shift(E,0);
		walltime=time(NULL);
		CPUtime=clock();

/*		complex<double> *v=new complex<double>[n];
		complex<double> *w=new complex<double>[n];
		for(int k=0;k<n;k++) v[k]=complex<double>(1.,0.5);
		MultInv(v,w);
		delete [] v;
		delete [] w;
*/
		ARCompStdEig<double, MatrixWithProduct2> dprob(n,stop, this, &MatrixWithProduct2::MultInv,"LM");
		if(LU_mode=="pardiso" || LU_mode=="superLU"){
			dprob.FindEigenvalues();
		}else{
			AREig(raw_evals,n,nonzero,vals,ja,ia,shift,stop);
		}
		walltime=time(NULL)-walltime;
		CPUtime=clock()-CPUtime;
		if(verbose>0) cout<<"the eigensolving took "<<(float)CPUtime/CLOCKS_PER_SEC<<" CPU time and "<<walltime<<" walltime"<<endl;
		
		eigvals=vector<double>(stop,0);
		//eigvecs=vector<vector< complex<double> > >(dprob.ConvergedEigenvalues(),vector<complex<double> > (n,0));
		for(int k=0;k<stop;k++){
			
			if(LU_mode=="pardiso" || LU_mode=="superLU") eigvals[k]=1./real(dprob.Eigenvalue(k))+E;
			else eigvals[k]=real(raw_evals[k]);
//			cout<<eigvals[k]<<endl;
//			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
		Nconverged=stop;
		sort(eigvals.begin(),eigvals.end());
	}
	//lowlevpos=sort_indexes(eigvals);
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}
//this should be called every time a different shift is required
void MatrixWithProduct2::release_after_LU(){

	if(LU_mode=="pardiso"){	
		mkl_free_buffers_();
	}else if(LU_mode=="superLU"){
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
	}
}

MatrixWithProduct2::~MatrixWithProduct2(){
	delete [] vals;
	delete [] ia;
	delete [] ja;
	if(LU_mode=="pardiso"){
		int phase=0;
		int error=0;
		msglvl=0;
		if(verbose>0) cout<<"about to start releasing pardisos memory"<<endl;
		pardiso_(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, &idum, &idum, &idum, &nrhs, &idum, &msglvl, &ddum, &ddum, &error);
		if(error) cout<<"error releasing pardiso's memory"<<endl;
	}else if(LU_mode=="superLU"){
		delete [] perm_c;
		delete [] perm_r;
		StatFree(&stat);
	}else if(LU_mode=="arpack++")
		delete [] raw_evals;
	if(verbose>0) cout<<"deallocated successfully"<<endl;
}

