#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double gstPrimeHedrick_C(NumericVector strata, NumericMatrix hets, double gstestimate){
//function declarations
  double rowMeanC(NumericMatrix, int);
  
// unique strata calculations 
  NumericVector uniqueS = unique(strata);
  int sizeuniqueS = uniqueS.size();  //numstrata
  
// gst calculation
  double Hs = rowMeanC(hets, 1);
  double gstmax = (((sizeuniqueS - 1) * (1 - Hs)) / (sizeuniqueS - 1 + Hs));
  double est = gstestimate/gstmax;

  return(est);
}
