#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Count unique cell values in a grid (c++)
//'
//' C++ function that counts the number of unique values in a flat grid vector
//' @param vals [numeric][numeric]\cr the values to count (column-major flat
//'   vector)
//' @param nrow [integer(1)][integer]\cr number of rows in the grid
//' @param ncol [integer(1)][integer]\cr number of columns in the grid
//' @family count functions
//' @return A data.frame of the unique values and their number
//' @export
// [[Rcpp::export]]
DataFrame countCellValuesCpp(NumericVector &vals, int nrow, int ncol) {
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
  IntegerVector outCells(elements);
  for(int j = 0; j < elements; j++){
    outValues[j] = (int)uvals[j];
  }

  // count occurrences
  for(int i = 0; i < n; i++){
    if(NumericVector::is_na(vals[i])) continue;
    // binary search for position
    int lo = 0, hi = elements - 1, pos = -1;
    double v = vals[i];
    while(lo <= hi){
      int mid = (lo + hi) / 2;
      if(uvals[mid] == v){ pos = mid; break; }
      else if(uvals[mid] < v) lo = mid + 1;
      else hi = mid - 1;
    }
    if(pos >= 0) outCells[pos] += 1;
  }

  DataFrame out = DataFrame::create(Named("value") = outValues,
                                    Named("cells") = outCells);
  return(out);
}
