#include <Rcpp.h>
using namespace Rcpp;

//' Check whether cell values are binary (c++)
//'
//' C++ function that checks if all non-NA values in a numeric vector are
//' either 0 or 1.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values.
//' @return A single logical value.
//' @export
// [[Rcpp::export]]
bool isBinaryCpp(NumericVector &vals) {
  int n = vals.size();

  for(int i = 0; i < n; i++){
    if(NumericVector::is_na(vals[i])) continue;
    if(vals[i] != 0.0 && vals[i] != 1.0){
      return false;
    }
  }

  return true;
}
