#include <Rcpp.h>
#include <queue>
#include <vector>
#include <limits>
using namespace Rcpp;


//' Cost-distance and shortest path on a grid (C++)
//'
//' Dijkstra's algorithm on a regular grid with per-cell traversal costs.
//' Returns both the cumulative cost surface and a predecessor map for
//' path tracing.
//'
//' @param cost [numeric(.)][numeric]\cr flat vector of per-cell traversal
//'   costs (row-major order).  Values must be >= 0.  NA cells are
//'   impassable barriers.
//' @param from [integer(.)][integer]\cr 1-based cell indices of source
//'   cells.
//' @param nrow [integer(1)][integer]\cr number of rows in the grid.
//' @param ncol [integer(1)][integer]\cr number of columns in the grid.
//' @param diagonal [logical(1)][logical]\cr if TRUE, use 8-connectivity
//'   (diagonal moves cost sqrt(2) * average cost).  If FALSE, use
//'   4-connectivity (rook moves only).
//' @return A list with two elements:
//'   \describe{
//'     \item{dist}{Numeric vector of cumulative cost-distances from the
//'       nearest source.  Inf for unreachable cells.}
//'     \item{pred}{Integer vector of 1-based predecessor cell indices
//'       (-1 for source cells, 0 for unreachable cells).  Walk
//'       backwards from any target cell to reconstruct the path.}
//'   }
//' @family grid utility functions
//' @export
// [[Rcpp::export]]
List costDistanceCpp(NumericVector cost, IntegerVector from,
                     int nrow, int ncol, bool diagonal) {

  int n = nrow * ncol;

  // output vectors
  std::vector<double> dist(n, std::numeric_limits<double>::infinity());
  std::vector<int> pred(n, 0);  // 0 = unreachable

  // priority queue: (cost, cell_index)
  typedef std::pair<double, int> pdi;
  std::priority_queue<pdi, std::vector<pdi>, std::greater<pdi>> pq;

  // seed source cells

  for (int i = 0; i < from.size(); i++) {
    int idx = from[i] - 1;  // convert to 0-based
    if (idx >= 0 && idx < n && !NumericVector::is_na(cost[idx])) {
      dist[idx] = 0.0;
      pred[idx] = -1;  // source marker
      pq.push(std::make_pair(0.0, idx));
    }
  }

  // direction offsets: 4-connected (N, S, E, W)
  int dr4[] = {-1, 1, 0, 0};
  int dc4[] = {0, 0, 1, -1};
  // additional diagonal offsets (NE, NW, SE, SW)
  int dr8[] = {-1, 1, 0, 0, -1, -1, 1, 1};
  int dc8[] = {0, 0, 1, -1, 1, -1, 1, -1};

  int ndirs = diagonal ? 8 : 4;
  int* dr = diagonal ? dr8 : dr4;
  int* dc = diagonal ? dc8 : dc4;

  double sqrt2 = 1.41421356237;

  while (!pq.empty()) {
    pdi top = pq.top();
    pq.pop();
    double d = top.first;
    int u = top.second;

    // skip if we already found a shorter path
    if (d > dist[u]) continue;

    int r = u / ncol;
    int c = u % ncol;

    for (int i = 0; i < ndirs; i++) {
      int nr = r + dr[i];
      int nc_new = c + dc[i];

      if (nr < 0 || nr >= nrow || nc_new < 0 || nc_new >= ncol) continue;

      int v = nr * ncol + nc_new;
      if (NumericVector::is_na(cost[v])) continue;

      // edge cost = average of endpoint costs * distance factor
      double step = (cost[u] + cost[v]) / 2.0;
      if (i >= 4) step *= sqrt2;  // diagonal

      double new_dist = d + step;
      if (new_dist < dist[v]) {
        dist[v] = new_dist;
        pred[v] = u + 1;  // store as 1-based
        pq.push(std::make_pair(new_dist, v));
      }
    }
  }

  // convert dist to R numeric vector
  NumericVector rdist(n);
  IntegerVector rpred(n);
  for (int i = 0; i < n; i++) {
    rdist[i] = dist[i];
    rpred[i] = pred[i];
  }

  return List::create(
    Named("dist") = rdist,
    Named("pred") = rpred
  );
}
