/*
this can do all kinds of things to a matrix, all it needs is a matvec that inherits from it
*/

#ifndef MATPROD_H
#define MATPROD_H

//#include "lapacke.h"
#include <iostream>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/Core>
#include "arscomp.h"

using namespace std;

template<class ART>
class MatrixWithProduct {

 private:

	int m, n; // Number of rows and columns.
	double E1,E2;
	ART *dense;
	Eigen::SparseMatrix<ART> sparse;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<ART> > sparseLU_solver;
	int *ipiv;
	vector<double> eigvals;
	vector< vector<ART> > eigvecs;
	vector<int> lowlevpos;

 public:

  int nrows() { return m; }
  int ncols() { return n; }
	void setrows(int x){ m=x; n=x;}
	
	double getE(int a){return eigvals[lowlevpos[a]];}
	vector<ART> getEV(int a){return eigvecs[lowlevpos[a]];}

  virtual void MultMv(ART* v, ART* w) = 0;
  virtual Eigen::Matrix<ART, Eigen::Dynamic, 1> MultEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1>) = 0;//used so we can communicate with my lanczos, which is based on eigen

  void MultM2v(ART* v, ART* w);

  void MultInvDense(ART*v, ART* w);
  void makeDense();
  void denseLU();
  void denseSolve();
  void printDense();

  void MultInvSparse(ART*v, ART* w);
  void makeSparse();
  void sparseSolve();

  int eigenvalues(int k, int q, double E);
  double calcVarEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1> v);
    
  ~MatrixWithProduct();
  
  MatrixWithProduct(int nrows, double _E1=0, double _E2=0, int ncols = 0 )
  // Constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
    E1=_E1; E2=_E2;
    dense=NULL;
    ipiv=NULL;
  } // Constructor.
  MatrixWithProduct()
  // Constructor.
  {
    m = 1;
    n = 1;
    E1=0.; E2=0.;
    dense=NULL;
    ipiv=NULL;
  } // Constructor.
  
  vector<int> sort_indexes(const vector<double> &v);

}; // MatrixWithProduct

template<class ART>
void MatrixWithProduct<ART>::MultM2v(ART* v, ART* w)
{
	if(E1==0 && E2==0){
		cout<<"energies were never set so can't compute the square thingy"<<endl;
		exit(0);
	}
	double *w1=new double[this->ncols()];
	double *w2=new double[this->ncols()];
	MultMv(v,w1);
	for(int i=0;i<this->ncols();i++)
		w1[i]=w1[i]-E1*v[i];
	
	MultMv(w1,w);
	for(int i=0;i<this->ncols();i++)
		w[i]=-w[i]+E2*w1[i];
	delete [] w1;
	delete [] w2;

} //  MultM2v.

//turn the matvec into a dense matrix
template<class ART>
void MatrixWithProduct<ART>::makeDense(){
	dense=new ART[n*n];
	ART *v=new ART[n];
	ART *w=new ART[n];
	for(int i=0;i<n;i++){
		for(int j=0; j<n; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		MultMv(v,w);
		for(int j=0; j<n; j++){
			dense[i+j*n]=w[j];
		}
	}
	delete [] v;
	delete [] w;
	
}

template<class ART>
void MatrixWithProduct<ART>::printDense(){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++) cout<<dense[i+j*n]<<" ";
		cout<<endl;
	}
}
/*//compute the dense LU decomposition
template<class ART>
void MatrixWithProduct<ART>::denseLU(){
	ipiv=new int[n];
	LAPACKE_dsytrf(LAPACK_COL_MAJOR,'U',n,dense,n,ipiv);
}

//compute solution to linear system (i.e. multiply by the inverse)
template<class ART>
void MatrixWithProduct<ART>::MultInvDense(ART *v, ART *w){
	for(int i=0;i<n;i++) w[i]=v[i];
	LAPACKE_dsytrs(LAPACK_COL_MAJOR,'U',n,1,dense,n,ipiv,w,n);
}
//calls lapack on the dense matrix to get the solution
//since we can't call lapack with a template this needs to be changed depending on what ART is
template<class ART>
void MatrixWithProduct<ART>::denseSolve(){
	ART *w=new double[n];
	LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,dense,n,w); 

//	for(int i=0;i<n;i++){
//		if ( abs(w[i])<1)	cout<<w[i]<<endl;
//	}

	for(int i=0;i<n;i++) cout<<w[i]<<endl;
	delete [] w;
}
*/
template<class ART>
void MatrixWithProduct<ART>::makeSparse(){
	sparse.resize(n,n);
	vector<Eigen::Triplet<ART> > coeff;
	Eigen::Triplet<ART> temp;
	dense=new ART[n*n];
	ART *v=new ART[n];
	ART *w=new ART[n];
	for(int i=0;i<n;i++){
		for(int j=0; j<n; j++){
			if(i==j) v[j]=1;
			else v[j]=0;
			w[j]=0;
		}
		MultMv(v,w);
		for(int j=0; j<n; j++){
		//	if (w[j]!=0) temp=Eigen::Triplet<ART>(j,i,w[j]);
			if (abs(w[j])>1e-16) coeff.push_back( Eigen::Triplet<ART>(j,i,w[j]) );
		}
	}
	sparse.setFromTriplets(coeff.begin(), coeff.end() );
	delete [] v;
	delete [] w;
	//cout<<sparse<<endl;
//}
//template <class ART>
//void  MatrixWithProduct<ART>::sparseLU(){
	sparseLU_solver.compute(sparse);
	if(sparseLU_solver.info()!=0) {
	  // decomposition failed
	  cout<<"decomposition failed! "<<sparseLU_solver.info()<<endl;
	}	
}

template<class ART>
void MatrixWithProduct<ART>::MultInvSparse(ART *v, ART *w){

	Eigen::Matrix<ART, Eigen::Dynamic, 1> out(n);
	Eigen::Map <Eigen::Matrix<ART, Eigen::Dynamic, 1> > mapped_v(v,n);
	out=sparseLU_solver.solve(mapped_v);
	//for(int i=0;i<n;i++) w[i]=out(i);	
	w=out.data();
}


//this doesn't do what you think it does! It converts the sparse matrix to a dense one and solves the dense one
template <class ART>
void MatrixWithProduct<ART>::sparseSolve(){
	Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> dMat=Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>(sparse);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> > es(dMat);
	cout<<es.eigenvalues()<<endl;
}

template <class ART>
double MatrixWithProduct<ART>::calcVarEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1> evec){
	Eigen::Matrix<ART, Eigen::Dynamic, 1> w(n);
	double var,norm,eval;

	norm=evec.norm();
	evec/=norm;
	w=MultEigen(evec);
	eval=w.dot(evec); //should I also return this eigenvalue?
//	cout<<eval<<endl;
	w=w-eval*evec;
	return w.norm();
}

template <class ART> 
int MatrixWithProduct<ART>::eigenvalues(int stop, int q=1, double E=-100){
	vector<ART>temp(n,0);
	int Nconverged;
	if (E==-100){
		ARCompStdEig<double, MatrixWithProduct<ART> >  dprob(ncols(), stop, this, &MatrixWithProduct<ART>::MultMv,"SR",(int)0, 1e-10,1e6);
		dprob.FindEigenvectors();
		Nconverged=dprob.ConvergedEigenvalues();

		eigvals=vector<double>(Nconverged,0);
		eigvecs=vector<vector<ART> >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=dprob.Eigenvalue(k).real();
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
	}else{
		makeSparse();
		ARCompStdEig<double, MatrixWithProduct<ART> >  dprob(ncols(), stop, this, &MatrixWithProduct<ART>::MultMv,E,"SR",(int)0, 1e-10,1e6);
		dprob.FindEigenvalues();

		eigvals=vector<double>(dprob.ConvergedEigenvalues(),0);
		eigvecs=vector<vector<ART> >(dprob.ConvergedEigenvalues(),temp);
		for(int k=0;k<dprob.ConvergedEigenvalues();k++){
			eigvals[k]=1./dprob.Eigenvalue(k).real()+E;
			eigvecs[k]=*(dprob.StlEigenvector(k));
		}
		Nconverged=dprob.ConvergedEigenvalues();
	}
	
	lowlevpos=sort_indexes(eigvals);	
	
	//get the positions of the q smallest  eigenvectors
//	lowlevpos=vector<int>(q,Nconverged-1 );
//	for(int k=0;k<Nconverged;k++){
//		for(int j1=0;j1<q;j1++){
//			if(eigvals[k]<eigvals[lowlevpos[j1]]){
//				for(int j2=q-1;j2>j1;j2--) lowlevpos[j2]=lowlevpos[j2-1];
//				lowlevpos[j1]=k;
//				break;
//			}
//		}
//	}
	return Nconverged;
//	for(int i=0;i<Nconverged;i++) cout<<dprob.Eigenvalue(i)<<endl;
	
}
template <class ART>
vector<int> MatrixWithProduct<ART>::sort_indexes(const vector<double> &v) {

  // initialize original index locations
  vector<int> idx(v.size());
  for (int i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}	
//destructor, delete the dense matrices
template<class ART>
MatrixWithProduct<ART>::~MatrixWithProduct(){
	delete [] dense;
	delete [] ipiv;
}
#endif // MATPROD_H

