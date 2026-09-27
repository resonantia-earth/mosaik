#include <Rcpp.h>
using namespace Rcpp;

//' Get Matrix cell values (c++)
//'
//' C++ function that extracts the values of a matrix
//' @param mat [matrix(numeric)][matrix]\cr the matrix from which to extract the values.
//' @family extractor functions
//' @return numeric vector of the values in \code{mat}
//' @export
// [[Rcpp::export]]
NumericVector getValuesCpp(NumericMatrix &mat) {
  int mRows = mat.nrow(), mCols = mat.ncol();
  NumericVector values;

  for(int y = 0; y < mRows; y++){
    for(int x = 0; x < mCols; x++){
      values.push_back(mat(y, x));
    }
  }
  return(values);
}