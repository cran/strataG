#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double gstPrimeNei_C(NumericVector strata, NumericMatrix hets){
//function declarations
  double rowMeanC(NumericMatrix, int);
  
// unique strata calculations 
  NumericVector uniqueS = unique(strata);
  uniqueS.sort();
  int sizeuniqueS = uniqueS.size();  //numstrata

// count strata frequency 
  int numindividuals = strata.size();
  NumericVector StrataFreq(sizeuniqueS);
  for(int i = 0; i < numindividuals; i++){
        StrataFreq(strata[i]) += 1;
    }
  
// gst calculation
  double Hs = rowMeanC(hets, 1);
  double Ht = rowMeanC(hets, 2);

  double est = ((sizeuniqueS * (Ht - Hs)) / ((sizeuniqueS * Ht) - Hs));
  
  return(est);
}
