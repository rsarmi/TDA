#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <cmath>
#include <algorithm>
#include <string> 
#include <sstream>
#include <iterator>

using namespace Rcpp;
using namespace RcppParallel;


// ----------------  AUXILIAR PARA RBFDOT -----------------------
template <typename InputIterator1, typename InputIterator2>
inline double rbfdotAux(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2, const double sigma){
  // value to return
  double rval = 0;
  double varianza = pow(sigma, 2.0);
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  // for each input item
  while (it1 != end1) {
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    // accumulate if appropirate
    rval += pow(d1-d2, 2.0);
  }
  rval = exp(-rval/(2*varianza));
  return rval;
}
// ---------------------------------------------------------------



// -----------------------------------------------------  
// -----------------------------------------------------
//              LIBERÍA  UN SÓLO NÚCLEO
// -----------------------------------------------------
// -----------------------------------------------------

// -------------  FUNCIÓN BÁSICA DE CÁLCULO KERNEL -------------------
// [[Rcpp::export]]
double rbfdot(NumericVector x, NumericVector y, double sigma) {
  int n = x.size();
  double val = 0;
  for(int i = 0; i < n; ++i) {
    val += pow(x[i] - y[i], 2.0);
  }
  val = exp(-val/(2*pow(sigma, 2.0)));
  return val;  
}

// ---------------- RBFDOT NUEVO VECTOR A DATOS------------------
// [[Rcpp::export]]
NumericVector rbfdotProd(NumericVector x, NumericMatrix mat, double sigma) {
  NumericVector rmat(mat.nrow());
  for (int i = 0; i < mat.nrow(); i++) {
       rmat[i] = rbfdot(x, mat(i,_), sigma);
  }
  return rmat;
}

// ---------------- MATRIX DE KERNEL GAUSSIANO ----------------------
// [[Rcpp::export]]
NumericMatrix rbfdotMatrix(NumericMatrix mat, double sigma) {
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  for (int i = 0; i < mat.nrow(); i++) {
    for(int j = i; j < mat.nrow(); j++){
      rmat(i, j) = rbfdot(mat(j,_), mat(i,_), sigma);
      rmat(j, i) = rmat(i, j);
    }
  }
  return rmat;
}
// ------------------------------------------------------------


// ---------------   RBFDOT PRODUCTO VECTOR ----------------------
// [[Rcpp::export]]
NumericVector rbfdotCProd(NumericVector x, NumericMatrix dta, NumericMatrix kmatrix, double sigma) {
  // outpuvalue
  NumericVector kprodsC;
  // j-ésima columna de kmatrix
  NumericVector kprods(dta.nrow());
  NumericVector dtarow(dta.ncol());
  double rowMean;
  double nrows=kmatrix.nrow();
  double kprodsMean = std::accumulate(kprods.begin(), kprods.end(), 0)/nrows;
  double kmatrixSum = std::accumulate(kmatrix.begin(), kmatrix.end(), 0)/pow(nrows, 2.0);
  // aquí guardaremos el valor
  for(int k=0; k < dta.nrow(); k++){
    dtarow = dta.row(k);
    kprods[k] = rbfdotAux(x.begin(), x.end(), dtarow.begin(), sigma);
  }
  // necesitamos matriz de unos par acentrar
  for(int i=0; i<dta.nrow(); i++){
    dtarow = dta.row(i);    
    rowMean = std::accumulate(dtarow.begin(), dtarow.end(), 0)/nrows;
    kprodsC[i] = kprods[i] - kprodsMean - rowMean + kmatrixSum;
  }
  return kprodsC;
}
// ------------------------------------------------------------
  
// inline NumericVector nbhdDist(int npts, NumericMatrix dist){
//   NumericVector vals(dist.nrow());
//   NumericVector distrow(dist.nrow());
//   double nrows=dist.nrow();
//   for(int i=0; i < nrows; i++){
//     distrow = dist.row(i);
//     std::sort(distrow.begin(), distrow.end());
//     vals[i] = std::accumulate(distrow.begin(), distrow.begin() + npts, 0.0)/nrows;
//   }
//   return vals;
// }

  
// -------------------------------------------------------------------------------------------------------------------------
  
  
// -----------------------------------------------------  
// -----------------------------------------------------
//              LIBERÍA PARA EL PARALELO
// -----------------------------------------------------
// -----------------------------------------------------


// -----------------------------------------------------
//            RBFDOTMATRIX PARALLEL
// -----------------------------------------------------
struct KernelMatrix : public Worker {
  // input matrix to read from
  const RMatrix<double> mat;
  const double sigma;
  // output matrix to write to
  RMatrix<double> rmat;
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  KernelMatrix(const NumericMatrix mat, const double sigma, NumericMatrix rmat)
    : mat(mat), sigma(sigma), rmat(rmat) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {
        // rows we will operate on
        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2 = mat.row(j);
        // write to output matrix
        rmat(i,j) = rbfdotAux(row1.begin(), row1.end(), row2.begin(), sigma);
        rmat(j,i) = rmat(i,j);
      }
      rmat(i,i) = 1;
    }
  }
};
// [[Rcpp::export]]
NumericMatrix rbfdotMatrixParallel(NumericMatrix mat, const double sigma) {
  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  // create the worker
  KernelMatrix kernelMatrix(mat, sigma, rmat);
  // call it with parallelFor
  parallelFor(0, mat.nrow(), kernelMatrix);
  return rmat;
}
// -----------------------------------------------------


// -----------------------------------------------------
//             RBFDOTPROD   - PARALLEL
// -----------------------------------------------------
struct RbfdotProdParallel : public Worker {
  // input to read from
  const RVector<double> x;
  const RMatrix<double> dta;
  const double sigma;
  // output value
  RVector<double> prod;
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  RbfdotProdParallel(const NumericVector x, const NumericMatrix dta, double sigma, NumericVector prod)
    : x(x), dta(dta), sigma(sigma), prod(prod) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      RMatrix<double>::Row rowj = dta.row(j);
      prod[j] = rbfdotAux(rowj.begin(), rowj.end(), x.begin(), sigma);
    }
  }
};
// [[Rcpp::export]]
NumericVector rbfdotProdParallel(NumericVector x, NumericMatrix dta, const double sigma) {
  // allocate the matrix we will return
  NumericVector prod(dta.nrow());
  // create the worker
  RbfdotProdParallel parallelWorker(x, dta, sigma, prod);
  // call it with parallelFor
  parallelFor(0, dta.nrow(), parallelWorker);
  return prod;
}
// ------------------------------------------------------------------------

// -----------------------------------------------------
//                  RBFDOT CENTERED PARALLEL
// ------------------------------------------------------
NumericVector rbfdotProcCParallel(NumericVector x, NumericMatrix dta, NumericMatrix kmatrix, double sigma){
  NumericVector kprods = rbfdotProdParallel(x, dta, sigma);
  double nrows=dta.nrow();
  double kprodsMean = std::accumulate(kprods.begin(), kprods.end(), 0)/nrows;
  double kmatrixSum = std::accumulate(kmatrix.begin(), kmatrix.end(), 0)/pow(nrows, 2.0);
  double rowMean;
  NumericVector kprodsC(dta.nrow());
  NumericVector dtarow(dta.ncol());
  for(int i=0; i < dta.nrow(); i++){
    dtarow = dta.row(i);
    rowMean = std::accumulate(dtarow.begin(), dtarow.end(), 0)/nrows;
    kprodsC[i] = kprods[i] - kprodsMean - rowMean + kmatrixSum;
  }
  return kprodsC;
}


// -----------------------------------------------------

// -----------------------------------------------------
//           EUCLID DISTANCE PARALLEL
// -----------------------------------------------------
template <typename InputIterator1, typename InputIterator2>
inline double euclidAux(InputIterator1 begin1, InputIterator1 end1, InputIterator2 begin2){
  // value to return
  double rval = 0;
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  // for each input item
  while (it1 != end1) {
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    // accumulate if appropirate
    rval += pow(d1-d2, 2.0);
  }
  rval = sqrt(rval);
  return rval;
}
struct EuclidDistance : public Worker {
  // input matrix to read from
  const RMatrix<double> dta;
  // output matrix to write to
  RMatrix<double> rmat;
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  EuclidDistance(const NumericMatrix dta, NumericMatrix rmat)
    : dta(dta), rmat(rmat) {}
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {
        // rows we will operate on
        RMatrix<double>::Row row1 = dta.row(i);
        RMatrix<double>::Row row2 = dta.row(j);
        // write to output matrix
        rmat(i,j) = euclidAux(row1.begin(), row1.end(), row2.begin());
        rmat(j,i) = rmat(i,j);
      }
    }
  }
};
// [[Rcpp::export]]
NumericMatrix euclidDistParallel(NumericMatrix dta) {
  // allocate the matrix we will return
  NumericMatrix rmat(dta.nrow(), dta.nrow());
  // create the worker
  EuclidDistance parallelWorker(dta, rmat);
  // call it with parallelFor
  parallelFor(0, dta.nrow(), parallelWorker);
  return rmat;
}

