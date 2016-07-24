/*
 * uses MKLs sparse interface to do the matrix diagonalization
*/

#ifndef MATPROD2_H
#define MATPROD2_H
#define SUPERLU_INC_

#include <iostream>
#include <algorithm>
#include "utils.h"
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include "slu_zdefs.h"

using namespace std;

#include "areig.h"
#include "arscomp.h"
//#include "arssym.h"


extern"C"{
	//converter from dense to sparse
	void mkl_zdnscsr_(int *job, int *rows, int *cols, complex<double> *a, int *lda, complex<double> *dat, int *ja, int *ia, int *info);	
	//matrix multiplication
	void zgemm_(char *transa, char *transb, int *rows, int *cols, int *k, double *alpha, complex<double> *a, int *lda, complex<double> *b, int *ldb, double *beta, complex<double> *c, int *ldc);
	//sparse matrix-vector
	void mkl_cspblas_zcsrsymv_(char *uplo, int *rows, complex<double> *vals, int *ia, int *ja, complex<double> *v, complex<double> *w);
	//pardiso
	void pardisoinit_(int *pt, int *mtype, int *iparm);
	void pardiso_(int *pt, int *maxfct, int *mnum, int *mtype, int *phase, int *rows, complex<double> *a, int *ia, int *ja, int *perm, int *nrhs, int *iparm, int *msglvl, complex<double> *b, complex<double> *x, int *error);

	void mkl_free_buffers_();
	void mkl_mem_stat_();
	void mkl_disable_fast_mm_();
	void MKL_FreeBuffers();
}

class MatrixWithProduct2 {

 private:

	int n; // Number of rows and columns.
	int nonzero;
	double E1,E2;
	
	complex<double> *vals,*raw_evals;
	int *ia, *ja;

	//pardiso stuff
	int idum; complex<double> ddum;
	int maxfct, mnum, mtype, msglvl,nrhs;
	int pt[64];
	int iparm[64];

 //used to 'undo' shifts on the diagonal, so we can use the same matrix for all shift-inverts
	double oldE;

	bool first_time;

	//stuff for the calls to superlu
	int *perm_c, *perm_r;
	int relax, panel_size;
	SuperMatrix A,L,U;
	superlu_options_t options;
	SuperLUStat_t stat;
	string LU_mode;
	
 public:

	int verbose;
  int ncols() { return n; }
	void setrows(int x){ 
		n=x;
	}
	
	vector<double> eigvals;
	vector< vector< complex<double> > > eigvecs;
//	vector<int> lowlevpos;
//	double getE(int a){return eigvals[lowlevpos[a]];} //ARPACK sometimes returns eigenvalues in the wrong order, these functions correct that

	void set_mode(string mode);
	void CSR_from_Dense(Eigen::MatrixXcd);
	void CSR_from_Sparse(Eigen::SparseMatrix< complex<double> >);
	void analyze_pattern();	
	void ShiftInvert(double E);
	void print();
	
  void MultMv(complex<double> *v, complex<double> *w); //original Matvec
  void MultInv(complex<double> *v, complex<double> *w); //uses precomputed LU decomposition to solve a system when doing shift-invert

	Eigen::MatrixXd EigenDense;

  int eigenvalues(int k, double E=-100); //computes eigenvalues and eigenvectors using ARPACK, E should be -100 if you want the ground state. If E is not -100, does shift-invert around E
  double single_energy(string whichp); //finds only the highest or lowest state, used for estimating energy bounds before doing shift-invert to target a section of the spectrum

	void release_after_LU();    
  ~MatrixWithProduct2();
  MatrixWithProduct2(int nrows )
  // Constructor.
  {
    n = nrows;
	verbose=0;
	oldE=0.;
	first_time=true;
	if(verbose>0) cout<<"hitting the right constructor"<<endl;
	mkl_disable_fast_mm_();
	LU_mode="superLU";
  } // Constructor.
  MatrixWithProduct2()
  // Constructor.
  {
	cout<<"hitting the wrong constructor"<<endl;
    n = 1;
  } // Constructor.
  
	
}; // MatrixWithProduct
#endif
