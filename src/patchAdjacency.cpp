#include <Rcpp.h>
#include <vector>
#include <map>
#include <utility>
using namespace Rcpp;

// Union-Find over contact cells, used to count spatially distinct contact
// regions per patch pair. Local to this translation unit.
static int pa_find(std::vector<int> &parent, int i) {
  while(parent[i] != i){
    parent[i] = parent[parent[i]];
    i = parent[i];
  }
  return i;
}

static void pa_union(std::vector<int> &parent, int a, int b) {
  int ra = pa_find(parent, a);
  int rb = pa_find(parent, b);
  if(ra != rb) parent[rb] = ra;
}

//' Patch-adjacency and contact regions on a labelled grid (c++)
//'
//' Given a grid of patch labels (e.g. from \code{componentsCpp}), count for
//' every ordered pair of distinct patches (i) how many adjacent cell pairs
//' connect them (contact length) and (ii) in how many spatially distinct
//' places they touch (contact regions). The second is the loop-vs-branch
//' discriminator for MSPA and is not derivable from the first.
//'
//' @param labels [integer(.)][integer]\cr flat vector of patch labels
//'   (row-major). NA cells are ignored.
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param eightconn [logical(1)][logical]\cr 8-connectivity (queen) if TRUE,
//'   4-connectivity (rook) if FALSE, for what counts as contact between two
//'   patches.
//' @return A list with \code{ids} (sorted unique patch labels, the row/column
//'   order of both matrices), \code{adjacency} (double-counted cell-pair
//'   contacts), and \code{regions} (distinct contact regions per pair).
//' @details Contact-cell linkage for the region count is always 8-connected:
//'   two contact cells belonging to the same patch pair are the same region if
//'   they are queen-adjacent. This is independent of \code{eightconn}, which
//'   only governs whether two patches are considered in contact at all.
//' @noRd
// [[Rcpp::export]]
List patchAdjacencyCpp(IntegerVector &labels, int nrow, int ncol,
                       bool eightconn) {
  int n = labels.size();

  // sorted unique labels -> 0-based position
  std::vector<int> ids;
  for(int i = 0; i < n; i++){
    if(IntegerVector::is_na(labels[i])) continue;
    ids.push_back(labels[i]);
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  int k = ids.size();

  std::map<int,int> pos;
  for(int i = 0; i < k; i++) pos[ids[i]] = i;

  NumericMatrix adjacency(k, k);

  // Contact-counting offsets: forward half (right, bottom [, both forward
  // diagonals]); each undirected contact is then tallied in both directions.
  std::vector<std::pair<int,int>> offs = {{1,0},{0,1}};
  if(eightconn){ offs.push_back({1,1}); offs.push_back({-1,1}); }

  // For regions: every contact cell of an ordered pair (a,b) gets a running id,
  // unioned with same-(a,b) contact cells that are 8-adjacent. A "contact cell"
  // here is the focal cell of patch a that borders patch b. We record, per
  // (a,b), the contact-cell node id assigned to each grid cell (or -1).
  // Because a cell can border several different patches, we key the node store
  // by (pairKey, cellIndex).
  std::vector<int> ufParent;                      // union-find node store
  std::map<std::pair<long long,int>, int> nodeOf; // (pairKey, cellIdx) -> node

  auto pairKey = [&](int a, int b) -> long long {
    return (long long)a * (long long)k + (long long)b;
  };

  // 8-neighbour offsets for linking contact cells into regions
  const int rlx[8] = {-1,0,1,-1,1,-1,0,1};
  const int rly[8] = {-1,-1,-1,0,0,1,1,1};

  for(int y = 0; y < nrow; y++){
    for(int x = 0; x < ncol; x++){
      int a = labels[y * ncol + x];
      if(IntegerVector::is_na(a)) continue;
      int pa = pos[a];

      for(size_t o = 0; o < offs.size(); o++){
        int nx = x + offs[o].first;
        int ny = y + offs[o].second;
        if(nx < 0 || nx >= ncol || ny < 0 || ny >= nrow) continue;
        int b = labels[ny * ncol + nx];
        if(IntegerVector::is_na(b) || b == a) continue;
        int pb = pos[b];

        // contact length, both directions (double count)
        adjacency(pa, pb) += 1;
        adjacency(pb, pa) += 1;

        // register the two focal cells as contact-region nodes for their
        // ordered pairs: cell (x,y) borders b -> node for (a,b);
        // cell (nx,ny) borders a -> node for (b,a).
        int idxA = y * ncol + x;
        int idxB = ny * ncol + nx;

        auto getNode = [&](int p, int q, int cellIdx) -> int {
          std::pair<long long,int> key(pairKey(p, q), cellIdx);
          auto it = nodeOf.find(key);
          if(it != nodeOf.end()) return it->second;
          int id = ufParent.size();
          ufParent.push_back(id);
          nodeOf[key] = id;
          return id;
        };
        getNode(pa, pb, idxA);
        getNode(pb, pa, idxB);
      }
    }
  }

  // Link contact cells of the same ordered pair that are 8-adjacent.
  for(auto &kv : nodeOf){
    long long key = kv.first.first;
    int cellIdx = kv.first.second;
    int node = kv.second;
    int cy = cellIdx / ncol;
    int cx = cellIdx % ncol;
    for(int d = 0; d < 8; d++){
      int nx = cx + rlx[d];
      int ny = cy + rly[d];
      if(nx < 0 || nx >= ncol || ny < 0 || ny >= nrow) continue;
      std::pair<long long,int> nkey(key, ny * ncol + nx);
      auto it = nodeOf.find(nkey);
      if(it != nodeOf.end()) pa_union(ufParent, node, it->second);
    }
  }

  // Count distinct region roots per ordered pair.
  // Map (pairKey) -> set of roots seen.
  NumericMatrix regions(k, k);
  std::map<long long, std::vector<int>> rootsPerPair;
  for(auto &kv : nodeOf){
    long long key = kv.first.first;
    int root = pa_find(ufParent, kv.second);
    rootsPerPair[key].push_back(root);
  }
  for(auto &kv : rootsPerPair){
    long long key = kv.first;
    std::vector<int> &roots = kv.second;
    std::sort(roots.begin(), roots.end());
    roots.erase(std::unique(roots.begin(), roots.end()), roots.end());
    int a = key / k;
    int b = key % k;
    regions(a, b) = roots.size();
  }

  IntegerVector idsOut(k);
  for(int i = 0; i < k; i++) idsOut[i] = ids[i];

  return List::create(Named("ids") = idsOut,
                      Named("adjacency") = adjacency,
                      Named("regions") = regions);
}
