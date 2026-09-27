#include <Rcpp.h>
using namespace Rcpp;

//' Scale values in a matrix (c++)
//'
//' C++ function that scales the values in a matrix between a provided range
//' @param mat [matrix(numeric)][matrix]\cr object in which the values are
//'   scaled.
//' @param range [numeric(2)][numeric]\cr two values (minimum and maximum)
//'   between which the values in \code{mat} shall be scaled.
//' @family matrix modify functions
//' @return a numeric matrix
//' @export
// [[Rcpp::export]]
NumericMatrix scaleMatrixCpp(NumericMatrix mat, NumericVector range){
  int mRows = mat.nrow(), mCols = mat.ncol();
  NumericMatrix out = clone(mat);
  double minVal = min(mat);
  double maxVal = max(mat);

  if(minVal == maxVal){
    return(out);
  }

  for(int y = 0; y < mRows; y++){
    for(int x = 0; x < mCols; x++){
      out(y, x) = (mat(y, x) - minVal) * (range[1] - range[0]) / (maxVal - minVal) + range[0];
    }
  }

  return(out);
}