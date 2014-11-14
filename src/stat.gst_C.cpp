#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double gst_C(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowMeanC(NumericMatrix, int);

// gst calculation
  double Hs = rowMeanC(hets, 1);
  double Ht = rowMeanC(hets, 2);
  double est = 1 - (Hs /Ht);
  return(est);  
}
