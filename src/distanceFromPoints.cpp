#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Distance from exact point coordinates to cell centroids (c++)
//'
//' For each cell in the grid, compute the minimum Euclidean distance to any of
//' the supplied source points. Coordinates are in CRS space (not cell indices).
//' @param xmin [numeric(1)][numeric]\cr left edge of the grid extent.
//' @param ymax [numeric(1)][numeric]\cr top edge of the grid extent.
//' @param res_x [numeric(1)][numeric]\cr cell width.
//' @param res_y [numeric(1)][numeric]\cr cell height.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param coords [matrix(numeric)][matrix]\cr two-column matrix of source
//'   point coordinates (x, y).
//' @return A numeric vector of length \code{ncol * nrow} with Euclidean
//'   distances in CRS units.
//' @noRd
// [[Rcpp::export]]
NumericVector distanceFromPointsCpp(double xmin, double ymax,
                                    double res_x, double res_y,
                                    int ncol, int nrow,
                                    NumericMatrix coords) {

  int n_cells = ncol * nrow;
  int n_pts = coords.nrow();
  NumericVector out(n_cells);

  for (int r = 0; r < nrow; r++) {
    double cy = ymax - (r + 0.5) * res_y;
    for (int c = 0; c < ncol; c++) {
      double cx = xmin + (c + 0.5) * res_x;
      int idx = r * ncol + c;

      double min_dist = R_PosInf;
      for (int p = 0; p < n_pts; p++) {
        double dx = cx - coords(p, 0);
        double dy = cy - coords(p, 1);
        double d = std::sqrt(dx * dx + dy * dy);
        if (d < min_dist) min_dist = d;
      }
      out[idx] = min_dist;
    }
  }

  return out;
}
