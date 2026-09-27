#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

//' Count cell adjacencies in a grid (c++)
//'
//' C++ function that counts the adjacencies between unique values in a flat
//' grid vector
//' @param vals [numeric][numeric]\cr the values to analyse (row-major flat
//'   vector)
//' @param nrow [integer(1)][integer]\cr number of rows in the grid
//' @param ncol [integer(1)][integer]\cr number of columns in the grid
//' @param doublecount [logical(1)][logical]\cr whether or not to tally up each
//'   adjacency for both neighbouring cells, or only once.
//' @param eightconn [logical(1)][logical]\cr whether to use 8-connectivity
//'   (queen: also count the four diagonal neighbours) or 4-connectivity (rook:
//'   orthogonal neighbours only).
//' @family count functions
//' @return A matrix of values*values with their adjacencies
//' @noRd
// [[Rcpp::export]]
NumericMatrix countCellAdjacenciesCpp(NumericVector &vals, int nrow, int ncol,
                                      bool doublecount, bool eightconn = false) {
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
  NumericMatrix out(elements, elements);

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

  // Neighbour offsets as (dx, dy). Counting the "forward" half (right, bottom,
  // and the two forward diagonals) visits every adjacency exactly once;
  // doublecount adds the mirrored half so each is tallied from both cells.
  std::vector<std::pair<int,int>> offsets = {{1, 0}, {0, 1}};   // right, bottom
  if(eightconn){
    offsets.push_back({1, 1});    // bottom-right
    offsets.push_back({-1, 1});   // bottom-left
  }
  if(doublecount){
    size_t half = offsets.size();
    for(size_t k = 0; k < half; k++){
      offsets.push_back({-offsets[k].first, -offsets[k].second});
    }
  }

  // iterate grid: row-major, index = y * ncol + x
  for(int y = 0; y < nrow; y++){
    for(int x = 0; x < ncol; x++){
      double focal = vals[y * ncol + x];
      if(NumericVector::is_na(focal)) continue;

      int posFocal = findPos(focal);
      if(posFocal < 0) continue;

      for(size_t k = 0; k < offsets.size(); k++){
        int nx = x + offsets[k].first;
        int ny = y + offsets[k].second;
        if(nx < 0 || nx >= ncol || ny < 0 || ny >= nrow) continue;
        double nVal = vals[ny * ncol + nx];
        if(NumericVector::is_na(nVal)) continue;
        int posN = findPos(nVal);
        if(posN >= 0) out(posFocal, posN) += 1;
      }
    }
  }

  return(out);
}
