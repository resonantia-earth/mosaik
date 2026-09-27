#include <Rcpp.h>
#include <functional>
using namespace Rcpp;


int edtFun(int x, int i, int gI){
  return(
    std::pow((x-i), 2) + std::pow(gI, 2)
  );
}

int mdtFun(int x, int i, int gI){
  return(
    std::abs(x-i) + gI
  );
}

int cdtFun(int x, int i, int gI){
  return(
    std::max(std::abs(x-i), gI)
  );
}

double edtSep(double i, double x, int gI, int gX, int maxVal){
  return(
    std::floor(
      (std::pow(x, 2) - std::pow(i, 2) + std::pow(gX, 2) - std::pow(gI, 2)) / (2*(x - i))
    )
  );
}

double mdtSep(double i, double x, int gI, int gX, int maxVal){
  if(gX >= (gI + x - i)){
    return(
      maxVal
    );
  } else if(gI > (gX + x - i)){
    return(
      -maxVal
    );
  } else{
    return(std::floor((gX - gI + x + i) /2));
  }
}

double cdtSep(double i, double x, int gI, int gX, int maxVal){
  if(gI <= gX){
    return(
      std::max(i + gX, floor((i + x)/2))
    );
  } else{
    return(
      std::min(x - gI, floor((i + x)/2))
    );
  }
}

// [[Rcpp::plugins(cpp11)]]


//' Calculate a distance map (c++)
//'
//' C++ function that derives a distance map of foreground patches in a flat
//' grid vector.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values
//'   (row-major order). Foreground cells have value 1.
//' @param nrow [integer(1)][integer]\cr number of rows in the grid.
//' @param ncol [integer(1)][integer]\cr number of columns in the grid.
//' @param method [character(1)][character]\cr which of the methods
//'   \code{"euclidean"}, \code{"manhattan"} or \code{"chessboard"} to use.
//' @family grid modify functions
//' @return A numeric vector of the same length with distance values.
//' @export
// [[Rcpp::export]]
NumericVector distanceCpp(NumericVector &vals, int nrow, int ncol, String method){

  int maxVal = ncol + nrow;
  std::vector<double> temp(nrow * ncol);
  NumericVector out(nrow * ncol);

  std::function<int(int, int, int)> fun;
  std::function<double(double, double, int, int, int)> sep;

  if(method == "euclidean"){
    fun = edtFun;
    sep = edtSep;
  } else if(method == "manhattan"){
    fun = mdtFun;
    sep = mdtSep;
  } else if(method == "chessboard"){
    fun = cdtFun;
    sep = cdtSep;
  }

  // first phase: scan columns
  // row-major: index = row * ncol + col
  for(int x = 0; x < ncol; x++){

    if(vals[0 * ncol + x] == 1){
      temp[0 * ncol + x] = 0;
    } else {
      temp[0 * ncol + x] = maxVal;
    }

    // scan 1 (top to bottom)
    for(int y = 1; y <= nrow - 1; y++){
      if(vals[y * ncol + x] == 1){
        temp[y * ncol + x] = 0;
      } else{
        temp[y * ncol + x] = 1 + temp[(y - 1) * ncol + x];
      }
    }

    // scan 2 (bottom to top)
    for(int y = nrow - 2; y >= 0; y--){
      if(temp[(y + 1) * ncol + x] < temp[y * ncol + x]){
        temp[y * ncol + x] = 1 + temp[(y + 1) * ncol + x];
      }
    }
  }

  // second phase: scan rows
  for(int y = 0; y <= nrow - 1; y++){

    std::vector<double> s(ncol);
    std::vector<double> t(ncol);

    int q = 0;
    s[0] = 0;
    t[0] = 0;

    // scan 3
    for(int x = 1; x <= ncol - 1; x++){

      int Fi = fun(t[q], s[q], temp[y * ncol + (int)s[q]]);
      int Fu = fun(t[q], x, temp[y * ncol + x]);

      while((q >= 1) & (Fi > Fu)){
        q--;
        Fi = fun(t[q], s[q], temp[y * ncol + (int)s[q]]);
        Fu = fun(t[q], x, temp[y * ncol + x]);
      }
      if(q < 0){
        q = 0;
        s[q] = x;
      } else{

        int w = 1 + sep(s[q], x, temp[y * ncol + (int)s[q]], temp[y * ncol + x], maxVal);
        if(w < ncol){
          q++;
          s[q] = x;
          t[q] = w;
        }
      }
    }

    // scan 4
    for(int x = ncol - 1; x >= 0; x--){
      int dis = fun(x, s[q], temp[y * ncol + (int)s[q]]);
      out[y * ncol + x] = dis;
      if(x == t[q]){
        q--;
      }
    }
  }

  return(out);
}
