#include <Rcpp.h>
using namespace Rcpp;

//' Binarise cell values by a threshold (c++)
//'
//' C++ function that converts a numeric vector to binary (0/1) values based on
//' a threshold. Values >= thresh become 1, values < thresh become 0, NA stays
//' NA.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values.
//' @param thresh [numeric(1)][numeric]\cr threshold value.
//' @return A numeric vector of the same length with binary values.
//' @noRd
// [[Rcpp::export]]
NumericVector binariseCpp(NumericVector &vals, double thresh) {
  int n = vals.size();
  NumericVector out(n);

  for(int i = 0; i < n; i++){
    if(NumericVector::is_na(vals[i])){
      out[i] = NA_REAL;
    } else {
      out[i] = (vals[i] >= thresh) ? 1.0 : 0.0;
    }
  }

  return out;
}
