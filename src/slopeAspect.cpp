#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Slope and aspect from a DEM via Horn's (1981) 3x3 finite-difference method,
// the same estimator GDAL and GRASS use. Input elevation is a flat vector in
// row-major order (row r, column c at index r*ncol + c), matching the rest of
// the mosaik C++ backend.
//
// Slope is returned in radians (angle from horizontal). Aspect is returned in
// radians measured clockwise from north (0 = north-facing, pi/2 = east, ...),
// following the geographic convention; flat cells (zero gradient) get aspect
// -1 as a sentinel so callers can treat them as direction-less.
//
// xres / yres are the cell sizes in the SAME linear unit as the elevation
// values (e.g. metres). For geographic (degree) extents the slope is not
// physically meaningful; that is a user-data concern documented in R.
//
// [[Rcpp::export]]
List slopeAspectCpp(NumericVector elev, int nrow, int ncol,
                    double xres, double yres) {
  int n = nrow * ncol;
  NumericVector slope(n);
  NumericVector aspect(n);

  // Horn weights the 8 neighbours; edges clamp to the nearest in-grid cell
  // (replicate-edge), so border cells get a one-sided estimate rather than NA.
  for (int r = 0; r < nrow; r++) {
    for (int c = 0; c < ncol; c++) {
      int rm = (r > 0) ? r - 1 : 0;
      int rp = (r < nrow - 1) ? r + 1 : nrow - 1;
      int cm = (c > 0) ? c - 1 : 0;
      int cp = (c < ncol - 1) ? c + 1 : ncol - 1;

      // 3x3 window (a b c / d e f / g h i), row-major indexing
      double a = elev[rm * ncol + cm];
      double b = elev[rm * ncol + c ];
      double cc= elev[rm * ncol + cp];
      double d = elev[r  * ncol + cm];
      double f = elev[r  * ncol + cp];
      double g = elev[rp * ncol + cm];
      double h = elev[rp * ncol + c ];
      double i = elev[rp * ncol + cp];

      int idx = r * ncol + c;

      if (NumericVector::is_na(elev[idx])) {
        slope[idx]  = NA_REAL;
        aspect[idx] = NA_REAL;
        continue;
      }

      // dz/dx grows eastward, dz/dy southward (rows run north to south)
      double dzdx = ((cc + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * xres);
      double dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + cc)) / (8.0 * yres);

      double mag = std::sqrt(dzdx * dzdx + dzdy * dzdy);
      slope[idx] = std::atan(mag);

      if (mag == 0.0) {
        aspect[idx] = -1.0;  // flat: no defined aspect
      } else {
        // the gradient points uphill; the slope faces the opposite way
        double az = std::atan2(dzdx, dzdy);
        az = az + M_PI;
        if (az < 0.0)        az += 2.0 * M_PI;
        if (az >= 2.0 * M_PI) az -= 2.0 * M_PI;
        aspect[idx] = az;
      }
    }
  }

  return List::create(Named("slope") = slope, Named("aspect") = aspect);
}
