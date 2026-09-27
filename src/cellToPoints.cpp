#include <Rcpp.h>
using namespace Rcpp;

//' Cell coordinates to data.frame (c++)
//'
//' C++ function that converts the cell coordinates of a matrix to a data.frame.
//' @param mat [matrix(numeric)][matrix]\cr object of which the cell coordinates
//'   are converted to a data.frame.
//' @return a data.frame of the cell coordinates
//' @export
// [[Rcpp::export]]
DataFrame cellToPointsCpp(NumericMatrix mat) {
  int mRows = mat.nrow(), mCols = mat.ncol();
  NumericVector xVal, yVal;
  IntegerVector values;

  for(int y = 0; y < mRows; y++){
    for(int x = 0; x < mCols; x++){

      xVal.push_back(x+0.5);
      yVal.push_back(mRows-(y+1)+0.5);
      values.push_back(mat(y, x));

    }
  }

  DataFrame out = DataFrame::create(Named("x")=xVal,
                                    Named("y")=yVal,
                                    Named("value")=values);

  return(out);
}
