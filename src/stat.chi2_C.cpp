#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]] 
double chi2H_C(NumericVector locusdata, NumericVector strata) {

// function declarations
  double calcchi2(NumericMatrix, NumericVector, NumericVector, int);
  NumericVector colSumsC(NumericMatrix);
  NumericVector rowSumsC(NumericMatrix);

// create new vectors containing the unique options
  NumericVector uniqueLD = unique(locusdata); 
  NumericVector uniqueS = unique(strata);
  uniqueLD.sort();
  uniqueS.sort();
  
// find size of each vector
  int sizeuniqueLD = uniqueLD.size();
  int sizeuniqueS = uniqueS.size();
  int sizedata = locusdata.size(); // number of individuals

// create blank table "obsfreq" for chi-square test
  NumericMatrix obsfreq(sizeuniqueLD, sizeuniqueS);
    
// counting obsfreq
      for (int i = 0; i < sizedata; i++){
          obsfreq(locusdata[i],strata[i]) += 1; 
      }    
       
// rowSums&colSums of obsfreq 
    NumericVector rowSums = rowSumsC(obsfreq);
    NumericVector colSums = colSumsC(obsfreq);

// Chi2 function call
    double chi2 = calcchi2(obsfreq, rowSums, colSums, sizedata);
 
   return(chi2); 
}

// [[Rcpp::export]]
double chi2D_C(NumericMatrix locusdata, NumericVector stratadata){

// function declarations
  double calcchi2(NumericMatrix, NumericVector, NumericVector, int);
  NumericVector colSumsC(NumericMatrix);
  NumericVector rowSumsC(NumericMatrix);
 
// Strata calculations
    int n = stratadata.size();
    NumericVector strata(n*2); // new vector used in frequency calculations
    for (int i = 0; i < n; i++){
      strata(i*2)     = stratadata(i);
      strata((i*2)+1) = stratadata(i);
    }

// Calculation for Number of Unique Strata
    NumericVector uniqueS = unique(stratadata);
    uniqueS.sort();
       
// Declared variables for loop  
  double chi2 = 0;
  int nrow = locusdata.nrow(), ncol = locusdata.ncol();

// Loop through columns (Every two columns of locusdata matrix is a different gene) 
  for (int column = 0; column < ncol; column+=2) {
    NumericVector templocusdata(nrow*2);
      for (int i = 0; i < nrow; i++){
        templocusdata(i*2)       = locusdata(i, column); 
        templocusdata((i*2) + 1) = locusdata(i, column+1);
      }

// Calculation for Number of Unique Alleles 
    NumericVector uniqueLD = unique(templocusdata);
    uniqueLD.sort();
    
// Remove NA from unique alleles
    for (int i = 0; i < uniqueLD.size(); i++){
      if(uniqueLD[i]== -1){
        uniqueLD.erase(i);
        break;
      }
    }

// size of each vector
      int sizeuniqueLD = uniqueLD.size();
      int sizeuniqueS = uniqueS.size();
      int sizedata = templocusdata.size();
    
// create blank table "obsfreq" for chi-square test
      NumericMatrix obsfreq(sizeuniqueLD, sizeuniqueS);
 
// counting obsfreq
      for (int i = 0; i < sizedata; i++){
         if(templocusdata[i]== -1){
            continue;
          }
          obsfreq(templocusdata[i], strata[i]) += 1;               
      } 
      
// rowSums&ColSums of obsfreq 
      NumericVector rowSums = rowSumsC(obsfreq);
      NumericVector colSums = colSumsC(obsfreq);

      int sum = 0;
      n = colSums.size();
      for (int i = 0; i< n; i++){
      sum += colSums[i];
      }
      
// Chi2 Calculation
       chi2 += calcchi2(obsfreq, rowSums, colSums, sum);
 
 }
  return(chi2);
}

// [[Rcpp::export]]
double calcchi2(NumericMatrix obsfreq, NumericVector rowSums, NumericVector colSums, int sizedata){
  int nrow = obsfreq.nrow();
  int ncol = obsfreq.ncol();
  double chi2 = 0;
  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      double expfreq = rowSums[i] * colSums[j] / sizedata;
      double temp = obsfreq(i,j) - expfreq;
      chi2 += temp*temp / expfreq;
    }
  }
  return(chi2);
  
}
