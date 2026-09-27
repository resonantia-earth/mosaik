#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Generate a height map (c++)
//'
//' C++ function that generates a height map based on the diamond-square
//' algorithm
//' @param mat [matrix(numeric)][matrix]\cr the template matrix into which
//'   height values are written.
//' @param stepSize [numeric(.)][numeric]\cr the side length of the squares at
//'   each iteration, largest first, each half the one before; the first is the
//'   side of \code{mat} minus one.
//' @param roughness [numeric(1)][numeric]\cr the Hurst exponent (fBm) or its
//'   equivalent, the roughness factor (DSa). Bounded between 0 and 1.
//' @param startDev [numeric(1)][numeric]\cr standard deviation of the random
//'   offset at the first iteration; it shrinks by \code{0.5^roughness} at each
//'   following one.
//' @family matrix modify functions
//' @return A numeric matrix of the same dimensions as \code{mat} into which the
//'   values of a height map have been written.
//' @export
// [[Rcpp::export]]
NumericMatrix diamondSquareCpp(NumericMatrix mat, NumericVector stepSize, double roughness, double startDev){

  const int n = stepSize.size();

  for(int z = 0; z < n; ++z){
    int halfStep = stepSize[z]/2;
    int aStep = stepSize[z];

    // diamond step
    // In this step the point in the center of a square is computed. This is
    // always a cell that exists and we can calculate the simple mean of the
    // corner values plus deviation.
    // a     b
    //
    //    X
    //
    // c     d
    for(int x = 0; x < stepSize[0]; x += aStep){
      for(int y = 0; y < stepSize[0]; y += aStep){
        float a = mat(y, x);
        float b = mat(y, x + aStep);
        float c = mat(y + aStep, x);
        float d = mat(y + aStep, x + aStep);
        NumericVector diam = (a + b + c + d)/4 + rnorm(1, 0, startDev);
        mat(y + halfStep, x + halfStep) = diam[0];
      }
    }

    // square step
    // In this step we calculate two points in the center of a diamond. The
    // base cell (a) is the same as in the previous step. e and f can be outsde
    // of the array. In these cases we need to assign them a value manually.
    // Here this is the average of the other three values with which e or f
    // would be averaged otherwise. Hence, the new value of e and f is
    // determined by their three respective neighbours only.
    //       e
    //
    //    a  X  b
    //
    // f  X  d
    //
    //    c
    //
    //
    for(int x = 0; x < stepSize[0]; x += aStep){
      for(int y = 0; y < stepSize[0]; y += aStep){
        float a = mat(y, x);
        float b = mat(y, x + aStep);
        float c = mat(y + aStep, x);

        float d = mat(y + halfStep, x + halfStep);
        float e = y - halfStep < 0 ? (a + b + d)/3 : mat(y - halfStep, x + halfStep);
        float f = x - halfStep < 0 ? (a + c + d)/3 : mat(y + halfStep, x - halfStep);

        NumericVector square1 = (a + b + d + e)/4 + rnorm(1, 0, startDev);
        NumericVector square2 = (a + c + d + f)/4 + rnorm(1, 0, startDev);
        mat(y, x + halfStep) = square1[0];
        mat(y + halfStep, x) = square2[0];
      }
    }
    startDev = startDev * std::pow(0.5, roughness);

  }

  return mat;
}
