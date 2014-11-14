#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector colSumsC(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector colSums(ncol);
  for (int i = 0; i < ncol; i++) {
    double total = 0;
    for (int j = 0; j < nrow; j++) {
      total += mat(j, i);
    }
    colSums[i] = total;
  }
  return(colSums);
}

// [[Rcpp::export]]
NumericVector rowSumsC(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  NumericVector rowSums(nrow);
  for (int i = 0; i < nrow; i++) {
    double total = 0;
    for (int j = 0; j < ncol; j++) {
      total += mat(i, j);
    }
    rowSums[i] = total;
  }
  return(rowSums);
}

// [[Rcpp::export]]
double colMeanC(NumericMatrix data, int columnToAverage) {
  int nrow = data.nrow();
  double sum = 0;
  for(int i = 0; i < nrow; i++){
    sum += data(i, columnToAverage);
  }
  double mean = sum/nrow;
  return(mean);
}

// [[Rcpp::export]]
double rowMeanC(NumericMatrix data, int rowToAverage) {
  int ncol = data.ncol();
  double sum = 0;
  for(int i = 0; i < ncol; i++){
    sum += data(rowToAverage, i);
  }
  double mean = sum/ncol;
  return(mean);
}
