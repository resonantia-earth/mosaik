#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Count cell edges in a grid (c++)
//'
//' C++ function that counts the number of edges between different values in a
//' flat grid vector
//' @param vals [numeric][numeric]\cr the values to analyse (row-major flat
//'   vector)
//' @param nrow [integer(1)][integer]\cr number of rows in the grid
//' @param ncol [integer(1)][integer]\cr number of columns in the grid
//' @family count functions
//' @return A data.frame of the unique values and their edges in X and Y
//' @noRd
// [[Rcpp::export]]
DataFrame countCellEdgesCpp(NumericVector &vals, int nrow, int ncol) {
  int n = vals.size();

  // collect unique non-NA values
  std::vector<double> uvals;
  for(int i = 0; i < n; i++){
    if(NumericVector::is_na(vals[i])) continue;
    bool found = false;
    for(size_t j = 0; j < uvals.size(); j++){
      if(vals[i] == uvals[j]){ found = true; break; }
    }
    if(!found) uvals.push_back(vals[i]);
  }
  std::sort(uvals.begin(), uvals.end());

  int elements = uvals.size();
  IntegerVector outValues(elements);
  IntegerVector edgesX(elements), edgesY(elements);
  for(int j = 0; j < elements; j++){
    outValues[j] = (int)uvals[j];
  }

  // helper to find position via binary search
  auto findPos = [&](double v) -> int {
    int lo = 0, hi = elements - 1;
    while(lo <= hi){
      int mid = (lo + hi) / 2;
      if(uvals[mid] == v) return mid;
      else if(uvals[mid] < v) lo = mid + 1;
      else hi = mid - 1;
    }
    return -1;
  };

  // iterate grid in row-major order: y = row, x = col
  // row-major: index = y * ncol + x
  for(int y = 0; y < nrow; y++){
    for(int x = 0; x < ncol; x++){
      int idx = y * ncol + x;
      double focal = vals[idx];
      if(NumericVector::is_na(focal)) continue;

      int pos = findPos(focal);
      if(pos < 0) continue;

      // right neighbour (x+1)
      if(x < ncol - 1){
        int rIdx = y * ncol + (x + 1);
        double rVal = vals[rIdx];
        if(!NumericVector::is_na(rVal) && rVal != focal){
          edgesX[pos] += 1;
        }
      }

      // left neighbour (x-1)
      if(x > 0){
        int lIdx = y * ncol + (x - 1);
        double lVal = vals[lIdx];
        if(!NumericVector::is_na(lVal) && lVal != focal){
          edgesX[pos] += 1;
        }
      }

      // bottom neighbour (y+1)
      if(y < nrow - 1){
        int bIdx = (y + 1) * ncol + x;
        double bVal = vals[bIdx];
        if(!NumericVector::is_na(bVal) && bVal != focal){
          edgesY[pos] += 1;
        }
      }

      // top neighbour (y-1)
      if(y > 0){
        int tIdx = (y - 1) * ncol + x;
        double tVal = vals[tIdx];
        if(!NumericVector::is_na(tVal) && tVal != focal){
          edgesY[pos] += 1;
        }
      }
    }
  }

  DataFrame out = DataFrame::create(Named("value") = outValues,
                                    Named("edgesX") = edgesX,
                                    Named("edgesY") = edgesY);
  return(out);
}
