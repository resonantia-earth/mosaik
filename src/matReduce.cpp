#include <Rcpp.h>
using namespace Rcpp;

//' Reduce a list of value vectors (c++)
//'
//' C++ function that uses a binary function to successively combine a list
//' of flat value vectors (all of equal length).
//' @param lVals [list(numeric)][list]\cr list of numeric vectors with equal
//'   length.
//' @param f [function(1)][function]\cr binary function to combine the vectors.
//' @return a single numeric vector of the same length as the input vectors.
//' @noRd
// [[Rcpp::export]]
NumericVector reduceCpp(List lVals, Function f) {
  int n = lVals.size();
  NumericVector out = clone(as<NumericVector>(lVals[0]));
  int len = out.size();

  for(int i = 1; i < n; i++){
    NumericVector vec = as<NumericVector>(lVals[i]);

    for(int j = 0; j < len; j++){
      NumericVector toAdd;
      toAdd.push_back(out[j]);
      toAdd.push_back(vec[j]);

      if(is_true(all(is_na(toAdd)))){
        out[j] = NA_REAL;
      } else{
        NumericVector theValue = f(na_omit(toAdd));
        out[j] = theValue[0];
      }
    }
  }

  return(out);
}
