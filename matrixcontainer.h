//** This object gives the full functionality of matprod (i.e. the ability to solve sparse and dense systems, including using shift-invert) 
//** to any Eigen Matrix object
//** Scott Geraedts, Oct 22, 2015
#ifndef MATRIXCONTAINER_H
#define MATRIXCONTAINER_H

#include "matprod.h"
#include "cblas.h"
#include <Eigen/Core>
using namespace std;

template<class ART>
class MatrixContainer: public MatrixWithProduct<ART> {

 private:
	Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> * H;
	Eigen::Matrix<ART, Eigen::Dynamic, 1> out;
	
 public:

	void MultMv(ART* v, ART* w);
	Eigen::Matrix<ART, Eigen::Dynamic, 1> MultEigen(Eigen::Matrix <ART, Eigen::Dynamic, 1>);

  MatrixContainer(int nv,Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic> &tH, double E1=0.,double E2=0.);

}; 

template<class ART>
inline MatrixContainer<ART>::MatrixContainer(int n, Eigen::Matrix<ART, Eigen::Dynamic, Eigen::Dynamic>& tH, double _E1, double _E2): MatrixWithProduct<ART>(n,_E1,_E2)
// Constructor
{
	H=&tH;
} // Constructor.

//  Matrix-vector multiplication w <- M*v.
template<class ART>
void MatrixContainer<ART>::MultMv(ART* v, ART* w){
//	Eigen::Map< Eigen::Matrix<ART, Eigen::Dynamic, 1> > wv(v, this->ncols(), 1);
//	for(int i=0;i<this->ncols();i++) cout<<wv(i)<<" ";
//	cout<<endl;
//	out=(*H)*wv;
//	for(int i=0;i<this->ncols();i++) cout<<out(i)<<" ";
//	cout<<endl;
//	w=out.data();
//	for(int i=0;i<this->ncols();i++) cout<<w[i]<<" ";
//	cout<<endl;
	
	double alpha=1., beta=0.;
	cblas_zgemv(CblasColMajor, CblasNoTrans, this->ncols(), this->ncols(), &alpha, (*H).data(), this->ncols(), v, 1, &beta, w, 1);
//	for(int i=0;i<this->ncols();i++) w[i]=0;
//	
//	for(int i=0;i<this->ncols(); i++){
//		for(int j=0;j<this->ncols(); j++)
//			w[i]+=(*H)(i,j)*v[j];
//	}	
} //  MultMv.

template<class ART>
Eigen::Matrix<ART, Eigen::Dynamic, 1> MatrixContainer<ART>::MultEigen(Eigen::Matrix<ART, Eigen::Dynamic, 1> v){
	Eigen::Matrix<ART, Eigen::Dynamic, 1> w=Eigen::Matrix<ART, Eigen::Dynamic, 1>::Zero(this->ncols());
	w=(*H)*v;
	return w;
} //  MultMv.
#endif 


