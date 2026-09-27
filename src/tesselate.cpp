#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Assign each cell to the nearest seed point (c++)
//'
//' For each cell centroid, find the closest seed point (by squared Euclidean
//' distance) and assign the corresponding label.
//' @param xmin [numeric(1)][numeric]\cr left edge of the grid extent.
//' @param ymax [numeric(1)][numeric]\cr top edge of the grid extent.
//' @param res_x [numeric(1)][numeric]\cr cell width.
//' @param res_y [numeric(1)][numeric]\cr cell height.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param seed_x [numeric][numeric]\cr x coordinates of seed points.
//' @param seed_y [numeric][numeric]\cr y coordinates of seed points.
//' @param seed_val [integer][integer]\cr label for each seed point.
//' @return An integer vector of length \code{ncol * nrow} with the label of the
//'   nearest seed point per cell.
//' @noRd
// [[Rcpp::export]]
IntegerVector tesselateCpp(double xmin, double ymax,
                           double res_x, double res_y,
                           int ncol, int nrow,
                           NumericVector seed_x, NumericVector seed_y,
                           IntegerVector seed_val) {

  int n_cells = ncol * nrow;
  int n_seeds = seed_x.size();
  IntegerVector out(n_cells);

  for (int r = 0; r < nrow; r++) {
    double cy = ymax - (r + 0.5) * res_y;
    for (int c = 0; c < ncol; c++) {
      double cx = xmin + (c + 0.5) * res_x;
      int idx = r * ncol + c;

      double best_d = R_PosInf;
      int best_v = NA_INTEGER;
      for (int s = 0; s < n_seeds; s++) {
        double dx = cx - seed_x[s];
        double dy = cy - seed_y[s];
        double d = dx * dx + dy * dy;  // no sqrt needed for comparison
        if (d < best_d) {
          best_d = d;
          best_v = seed_val[s];
        }
      }
      out[idx] = best_v;
    }
  }

  return out;
}
