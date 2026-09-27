#include <Rcpp.h>
using namespace Rcpp;

//' Rescale a gridded object by a factor (c++)
//'
//' C++ function that rescales a grid using nearest-neighbor resampling.
//' Operates on a flat vector (row-major) without requiring an R matrix.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values
//'   (row-major order).
//' @param nrow [integer(1)][integer]\cr number of rows in the input.
//' @param ncol [integer(1)][integer]\cr number of columns in the input.
//' @param factorRow [numeric(1)][numeric]\cr scaling factor for rows.
//' @param factorCol [numeric(1)][numeric]\cr scaling factor for columns.
//' @return A numeric vector of the rescaled grid in row-major order.
//'   The output dimensions are floor(nrow * factorRow) x floor(ncol *
//'   factorCol).
//' @details Uses nearest-neighbor interpolation. Each output cell maps back
//'   to the closest input cell. This is equivalent to
//'   \code{mmand::rescale()} with a box kernel for integer scaling factors.
//' @noRd
// [[Rcpp::export]]
List rescaleGridCpp(NumericVector &vals, int nrow, int ncol,
                    double factorRow, double factorCol) {

  int newRows = (int)(nrow * factorRow);
  int newCols = (int)(ncol * factorCol);
  int newN = newRows * newCols;

  NumericVector out(newN);

  for(int newRow = 0; newRow < newRows; newRow++){
    // map back to source row
    int srcRow = (int)(newRow / factorRow);
    if(srcRow >= nrow) srcRow = nrow - 1;

    for(int newCol = 0; newCol < newCols; newCol++){
      // map back to source column
      int srcCol = (int)(newCol / factorCol);
      if(srcCol >= ncol) srcCol = ncol - 1;

      int srcIdx = srcRow * ncol + srcCol;
      int dstIdx = newRow * newCols + newCol;

      out[dstIdx] = vals[srcIdx];
    }
  }

  return List::create(
    Named("vals") = out,
    Named("nrow") = newRows,
    Named("ncol") = newCols
  );
}
