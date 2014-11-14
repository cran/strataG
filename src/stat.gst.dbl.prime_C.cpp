#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double gstDblPrime_C(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowMeanC(NumericMatrix, int);
  double gstPrimeNei_C(NumericVector, NumericMatrix);
    
// Calculate primenei
  double primenei = gstPrimeNei_C(strata, hets);
      
// double prime calculation
  double Hs = rowMeanC(hets, 1);
  double est = (primenei / (1 - Hs));
    
  return(est);
}
