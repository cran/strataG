#include <Rcpp.h>
#include <vector>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double phist_C(NumericVector locusdata, NumericVector strata, NumericMatrix hapdist) {

// basic calculations of Strata, Locusdata, Unique, Size of vectors
  int numindividuals = locusdata.size();
  NumericVector uniqueS = unique(strata);
  uniqueS.sort();
  int sizeuniqueS = uniqueS.size();
  NumericVector uniqueLD = unique(locusdata);
  uniqueLD.sort();
  int sizeuniqueLD = uniqueLD.size();

// counting number in each strata (strataFreq), and number from each haplotype in each strata (strataHapFreq)
  NumericVector strataFreq(sizeuniqueS);
  NumericMatrix strataHapFreq(sizeuniqueLD, sizeuniqueS);
  for(int i = 0; i < numindividuals; i++) {
    strataFreq(strata(i))++;
    strataHapFreq(locusdata[i], strata[i])++;
  }
  
// Calculate sums of squares within strata (Eqn 8a)  
  double ssdwp = 0; 
  for(int column = 0; column < sizeuniqueS; column++) {
    NumericVector hapfreqS(sizeuniqueLD);
    NumericVector uniqueHap(sizeuniqueLD); 
    
    for(int i = 0; i < sizeuniqueLD; i ++) {
      if(strataHapFreq(i, column) > 0) {
        hapfreqS(i) = strataHapFreq(i,column);
        uniqueHap(i) = uniqueLD(i);    
      }
      else {
        uniqueHap(i) = -1;
      }
    }

    for (int i = uniqueHap.size(); i >= 0; i--) {
      if(uniqueHap[i] == -1) {
        uniqueHap.erase(i);
        hapfreqS.erase(i);
      }
    }    
    
    double sum = 0;
    int sizeUhaps = uniqueHap.size();
    for(int h1 = 0; h1 < sizeUhaps; h1++) {
      for(int h2 = 0; h2 < sizeUhaps; h2++) {
        sum += (hapdist(uniqueHap(h1), uniqueHap(h2)) * hapfreqS(h1) * hapfreqS(h2)); 
      }
    }
    sum = sum / (2 * strataFreq(column));
    ssdwp += sum;
  }

// Calculate sums of squares amoung strata (Eqn 8a)  
  double ssdap = 0; 
  for(int column = 0; column < sizeuniqueS; column++) {
    for(int column2 = 0; column2 < sizeuniqueS; column2++) {
      NumericVector hapfreqS1(sizeuniqueLD);
      NumericVector uniqueHap1(sizeuniqueLD);       
      NumericVector hapfreqS2(sizeuniqueLD);
      NumericVector uniqueHap2(sizeuniqueLD); 
      
// fill from first option
      for(int i = 0; i < sizeuniqueLD; i ++) {
        if(strataHapFreq(i, column) > 0) {
          hapfreqS1(i) = strataHapFreq(i,column);
          uniqueHap1(i) = uniqueLD(i);          
        }
        else {
          uniqueHap1(i) = -1;
        }
      }
      
// fill from second option
      for(int i = 0; i < sizeuniqueLD; i ++) {
        if(strataHapFreq(i, column2) > 0) {
          hapfreqS2(i) = strataHapFreq(i,column2);
          uniqueHap2(i) = uniqueLD(i);          
        }
        else {
          uniqueHap2(i) = -1;
        }
      }
      
// remove -1 from first and second options
      for (int i = uniqueHap1.size(); i >= 0; i--) {
        if(uniqueHap1[i] == -1) {
          uniqueHap1.erase(i);
          hapfreqS1.erase(i);
        }
      }    
      for (int i = uniqueHap2.size(); i >= 0; i--) {
        if(uniqueHap2[i] == -1) {
          uniqueHap2.erase(i);
          hapfreqS2.erase(i);
        }
      }  
      
// calculate ssdap
      double sum = 0;
      int sizeUhaps1 = uniqueHap1.size();
      int sizeUhaps2 = uniqueHap2.size();
      for(int h1 = 0; h1 < sizeUhaps1; h1++) {
        for(int h2 = 0; h2 < sizeUhaps2; h2++) {
          sum += (hapdist(uniqueHap1(h1), uniqueHap2(h2)) * hapfreqS1(h1) * hapfreqS2(h2)); 
        }
      }
      sum = sum / (2 * numindividuals);
      ssdap += sum;
    }
  }
  ssdap = ssdap - ssdwp;
   
//  # Calculate average sample size correction for among strata variance 
//  #  Eqn 9a in paper, but modified as in Table 8.2.1.1 from Arlequin v3.5.1 manual
//  #  (denominator is sum{I} - 1)
  double strataVariance = 0;
  for(int i = 0; i < sizeuniqueS; i++) {
    strataVariance += (strataFreq(i))*strataFreq(i) / numindividuals / (sizeuniqueS - 1);
  }
  double n = (numindividuals - strataVariance);
  
//  # Calculate variance components (Table 1)
//  #   Set MSD (SSD / df) equal to expected MSD
  double Vc = ssdwp / (numindividuals - sizeuniqueS);
  double Vb = ((ssdap / (sizeuniqueS - 1)) - Vc) / n;
  
  double est = Vb / (Vb + Vc);

  return(est);
}