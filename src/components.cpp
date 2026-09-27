#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// Union-Find (disjoint set) helpers
static int uf_find(std::vector<int> &parent, int i) {
  while(parent[i] != i){
    parent[i] = parent[parent[i]]; // path compression
    i = parent[i];
  }
  return i;
}

static void uf_union(std::vector<int> &parent, std::vector<int> &rank, int a, int b) {
  int ra = uf_find(parent, a);
  int rb = uf_find(parent, b);
  if(ra == rb) return;
  if(rank[ra] < rank[rb]) std::swap(ra, rb);
  parent[rb] = ra;
  if(rank[ra] == rank[rb]) rank[ra]++;
}

//' Connected-component labeling (c++)
//'
//' C++ function that labels connected components in a gridded object using a
//' two-pass union-find algorithm.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values (row-major
//'   order).
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param connectivity [integer(1)][integer]\cr neighbourhood connectivity;
//'   \code{4} for rook/diamond (up/down/left/right) or \code{8} for
//'   queen/box (including diagonals). Default is \code{4}.
//' @return An integer vector of the same length where each connected component
//'   of non-NA cells sharing the same value has a unique label. NA cells remain
//'   NA.
//' @noRd
// [[Rcpp::export]]
IntegerVector componentsCpp(NumericVector &vals, int nrow, int ncol,
                            int connectivity = 4) {
  int n = vals.size();
  IntegerVector labels(n, NA_INTEGER);

  // parent and rank arrays for union-find
  std::vector<int> parent(n + 1);
  std::vector<int> rnk(n + 1, 0);
  for(int i = 0; i <= n; i++) parent[i] = i;

  int nextLabel = 1;

  // Pass 1: assign provisional labels and record equivalences
  // row-major: index = row * ncol + col
  for(int row = 0; row < nrow; row++){
    for(int col = 0; col < ncol; col++){

      int idx = row * ncol + col;

      if(NumericVector::is_na(vals[idx])) continue;

      double val = vals[idx];
      int neighborLabel = 0;

      // check left neighbor (col - 1)
      if(col > 0){
        int left = row * ncol + (col - 1);
        if(!NumericVector::is_na(vals[left]) && vals[left] == val){
          if(neighborLabel == 0){
            neighborLabel = labels[left];
          } else {
            uf_union(parent, rnk, neighborLabel, labels[left]);
          }
        }
      }

      // check top neighbor (row - 1)
      if(row > 0){
        int top = (row - 1) * ncol + col;
        if(!NumericVector::is_na(vals[top]) && vals[top] == val){
          if(neighborLabel == 0){
            neighborLabel = labels[top];
          } else {
            uf_union(parent, rnk, neighborLabel, labels[top]);
          }
        }
      }

      // diagonal neighbors for 8-connectivity
      // Row-major scan (left-to-right, top-to-bottom): already-visited
      // backward neighbors are top-left, top-right, and left-col neighbors.
      if(connectivity == 8){
        // top-left (row - 1, col - 1)
        if(row > 0 && col > 0){
          int tl = (row - 1) * ncol + (col - 1);
          if(!NumericVector::is_na(vals[tl]) && vals[tl] == val){
            if(neighborLabel == 0){
              neighborLabel = labels[tl];
            } else {
              uf_union(parent, rnk, neighborLabel, labels[tl]);
            }
          }
        }

        // top-right (row - 1, col + 1)
        if(row > 0 && col < ncol - 1){
          int tr = (row - 1) * ncol + (col + 1);
          if(!NumericVector::is_na(vals[tr]) && vals[tr] == val){
            if(neighborLabel == 0){
              neighborLabel = labels[tr];
            } else {
              uf_union(parent, rnk, neighborLabel, labels[tr]);
            }
          }
        }
      }

      if(neighborLabel == 0){
        // new component
        labels[idx] = nextLabel;
        nextLabel++;
      } else {
        labels[idx] = uf_find(parent, neighborLabel);
      }
    }
  }

  // Pass 2: resolve all labels to their root
  // Also re-number sequentially (1, 2, 3, ...)
  std::vector<int> remap(nextLabel, 0);
  int finalLabel = 0;

  for(int i = 0; i < n; i++){
    if(IntegerVector::is_na(labels[i])) continue;

    int root = uf_find(parent, labels[i]);
    if(remap[root] == 0){
      finalLabel++;
      remap[root] = finalLabel;
    }
    labels[i] = remap[root];
  }

  return labels;
}
