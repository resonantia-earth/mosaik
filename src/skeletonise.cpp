#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' Morphological skeletonisation (c++)
//'
//' C++ function that computes the morphological skeleton of a binary gridded
//' object using the hit-or-miss thinning algorithm. Iteratively removes
//' boundary pixels while preserving topology (connectivity).
//' @param vals [numeric(.)][numeric]\cr flat vector of binary cell values
//'   (row-major order). Must contain only 0, 1, and NA.
//' @param nrow [integer(1)][integer]\cr number of rows.
//' @param ncol [integer(1)][integer]\cr number of columns.
//' @param anchor [numeric(.)][numeric]\cr optional flat vector, same length as
//'   \code{vals}, marking anchor cells (non-zero) that are never removed during
//'   thinning. This yields an anchored skeleton whose branches stay connected
//'   to the anchor set and whose free ends rest on it. Pass an all-zero (or
//'   empty) vector for the ordinary unanchored skeleton.
//' @param homotopic [logical(1)][logical]\cr if \code{TRUE}, order-independent
//'   homotopic thinning (Ranwez & Soille 2002) instead of Zhang-Suen.
//' @return An integer vector of the same length with the skeleton (1) and
//'   background (0). NA cells remain NA.
//' @details Zhang-Suen iterates two directional sub-passes. Homotopic thinning
//'   deletes simple points -- cells whose removal keeps both foreground and
//'   background connected as before -- but only those that are also
//'   independent of every simple neighbour, since two adjacent simple pixels
//'   may each be removable alone and not together. That is Ronse's strong
//'   8-deletability, and it is what makes the result independent of the scan
//'   order. Ported from miallib's skel.c (stype 0: foreground 8-connected,
//'   background 4-connected).
//' @export
// [[Rcpp::export]]
IntegerVector skeletoniseCpp(NumericVector &vals, int nrow, int ncol,
                             NumericVector anchor = NumericVector::create(),
                             bool homotopic = false) {
  int n = vals.size();

  // anchor lookup: a cell is anchored when the anchor vector marks it non-zero
  bool useAnchor = (anchor.size() == n);
  auto isAnchor = [&](int idx) -> bool {
    return useAnchor && !NumericVector::is_na(anchor[idx]) && anchor[idx] != 0.0;
  };

  // work on an integer grid (1 = foreground, 0 = background, -1 = NA)
  std::vector<int> grid(n);
  for(int i = 0; i < n; i++){
    if(NumericVector::is_na(vals[i])){
      grid[i] = -1;
    } else {
      grid[i] = (vals[i] != 0.0) ? 1 : 0;
    }
  }

  // helper to get value at (row, col), treating out-of-bounds as 0
  auto getVal = [&](int r, int c) -> int {
    if(r < 0 || r >= nrow || c < 0 || c >= ncol) return 0;
    int idx = r * ncol + c;
    return (grid[idx] == 1) ? 1 : 0;
  };

  // --- order-independent anchored homotopic thinning -----------------------
  //
  // Ported from ec-jrc/jeolib-miallib core/c/skel.c (binOIthin, stype = 0,
  // i.e. foreground 8-connected and background 4-connected). Deleting simple
  // points one at a time is ORDER DEPENDENT: two adjacent simple pixels may
  // each be removable alone but not together. A pixel is therefore deleted
  // only when it is simple AND independent of every simple neighbour, which is
  // Ronse's criterion of strong 8-deletability. That independence test is what
  // makes the result invariant to the scan order.
  //
  // References: Ranwez & Soille (2002); Iwanowski & Soille (2007).
  if(homotopic){

    // skel.c getngbshift8 order: W, E, N, S, NW, SW, NE, SE
    static const int dr[8] = { 0,  0, -1,  1, -1,  1, -1,  1};
    static const int dc[8] = {-1,  1,  0,  0, -1, -1,  1,  1};

    // skel.c simple_pixel, table for stype 0. The neighbourhood is encoded as
    // W | E<<1 | N<<2 | S<<3 | NW<<4 | SW<<5 | NE<<6 | SE<<7.
    static const unsigned char htab[256] = {
    0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,
    1,1,0,0,1,1,1,1,0,1,0,1,0,1,1,0,
    1,1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,
    0,1,0,0,0,1,0,1,0,1,0,1,0,1,1,0,
    1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,0,
    0,0,0,0,1,1,1,1,0,0,0,1,0,1,1,0,
    0,0,0,0,0,1,0,1,0,0,1,1,0,1,1,0,
    0,0,0,0,0,1,0,1,0,0,0,1,0,1,1,0,
    1,0,1,0,0,0,1,1,1,1,1,1,0,1,1,0,
    0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,0,
    0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,0,
    0,0,0,0,0,0,0,1,0,1,0,1,0,1,1,0,
    0,0,1,0,0,0,1,1,0,0,1,1,0,1,1,0,
    0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,
    0,0,0,0,0,0,0,1,0,0,1,1,0,1,1,0,
    0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0
    };

    // skel.c simple_pair, tables for stype 0
    static const int simplePair2[4]  = {0,1,1,1};
    static const int simplePair4[16] = {0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0};

    auto inBounds = [&](int r, int c) -> bool {
      return r >= 0 && r < nrow && c >= 0 && c < ncol;
    };

    auto isSimple = [&](int r, int c) -> bool {
      if(!inBounds(r, c) || grid[r * ncol + c] != 1) return false;
      int code = 0;
      for(int k = 0; k < 8; k++){
        int rr = r + dr[k], cc = c + dc[k];
        if(inBounds(rr, cc) && grid[rr * ncol + cc] == 1) code |= (1 << k);
      }
      return htab[code] != 0;
    };

    // skel.c test_anchor: true when the cell is NOT an anchor
    auto notAnchor = [&](int r, int c) -> bool {
      if(!inBounds(r, c)) return true;
      return !isAnchor(r * ncol + c);
    };

    // skel.c testsimple: neighbour k is unanchored, simple, and equal-valued
    auto testSimple = [&](int r, int c, int k) -> bool {
      int rr = r + dr[k], cc = c + dc[k];
      if(!inBounds(rr, cc)) return false;
      if(!notAnchor(rr, cc)) return false;
      if(!isSimple(rr, cc)) return false;
      int a = grid[r * ncol + c] == 1 ? 1 : 0;
      int b = grid[rr * ncol + cc] == 1 ? 1 : 0;
      return a == b;
    };

    // value at (r,c) as the C code reads it, out of bounds being background
    auto val = [&](int r, int c) -> int {
      if(!inBounds(r, c)) return 0;
      return grid[r * ncol + c] == 1 ? 1 : 0;
    };

    // skel.c simple_pair: is p independent of its neighbour k? Also returns the
    // configuration type used by the triple/quadruple test below.
    auto simplePair = [&](int r, int c, int k, int &conftype) -> int {
      int self = val(r, c);
      if(k < 4){
        // non-diagonal neighbour: four cells a, b, c, d around the pair
        static const int q[4][4] = {{2,4,3,5}, {2,6,3,7}, {0,4,1,6}, {0,5,1,7}};
        int ia = q[k][0], ib = q[k][1], ic = q[k][2], id = q[k][3];
        int code = 0;
        code |= (val(r + dr[id], c + dc[id]) >= self);
        code |= (val(r + dr[ic], c + dc[ic]) >= self) << 1;
        code |= (val(r + dr[ib], c + dc[ib]) >= self) << 2;
        code |= (val(r + dr[ia], c + dc[ia]) >= self) << 3;
        conftype = code + 5;
        return simplePair4[code];
      } else {
        // diagonal neighbour: two cells a, b flanking the pair
        static const int q[4][2] = {{0,2}, {3,0}, {2,1}, {1,3}};
        int ia = q[k-4][0], ib = q[k-4][1];
        int code = 0;
        code |= (val(r + dr[ib], c + dc[ib]) >= self);
        code |= (val(r + dr[ia], c + dc[ia]) >= self) << 1;
        conftype = code + 1;
        return simplePair2[code];
      }
    };

    // skel.c num_sngb: how many of the first n neighbours are >= the centre
    auto numSngb = [&](int r, int c, int n) -> int {
      int self = val(r, c), sum = 0;
      for(int i = 0; i < n; i++) if(val(r + dr[i], c + dc[i]) >= self) sum++;
      return sum;
    };

    // skel.c indep_simple, stype 0
    auto indepSimple = [&](int r, int c) -> bool {
      int allindep = 1, ngbnum = 0, ct = 0, conftype = 0;
      for(int k = 0; k < 8 && allindep; k++){
        if(testSimple(r, c, k)){
          ngbnum++;
          if(simplePair(r, c, k, conftype) == 0) allindep = 0;
          ct = conftype;
        }
      }
      if(allindep){
        // p may still belong to a triple or quadruple of mutually simple
        // pixels that cannot all go at once
        bool susp = (ngbnum == 2 && (ct == 2 || ct == 3 || ct == 6 ||
                                     ct == 7 || ct == 9 || ct == 13)) ||
                    (ngbnum == 3 && (ct == 4 || ct == 8 || ct == 17));
        if(susp){
          int founddiff = 0;
          for(int k = 0; k < 8 && !founddiff; k++){
            if(testSimple(r, c, k)){
              if(numSngb(r + dr[k], c + dc[k], 8) != ngbnum) founddiff = 1;
            }
          }
          if(!founddiff) allindep = 0;
        }
      }
      return allindep != 0;
    };

    // detection pass then deletion pass, to a fixpoint. Because every marked
    // pixel is independent of the others, the whole set can be removed at once.
    std::vector<int> mark(n);
    while(true){
      std::fill(mark.begin(), mark.end(), 0);
      bool any = false;
      for(int r = 0; r < nrow; r++){
        for(int c = 0; c < ncol; c++){
          int idx = r * ncol + c;
          if(grid[idx] != 1) continue;
          if(isAnchor(idx)) continue;
          if(isSimple(r, c) && indepSimple(r, c)){ mark[idx] = 1; any = true; }
        }
      }
      if(!any) break;
      for(int i = 0; i < n; i++) if(mark[i]) grid[i] = 0;
    }

    IntegerVector outOI(n, NA_INTEGER);
    for(int i = 0; i < n; i++) outOI[i] = (grid[i] == -1) ? NA_INTEGER : grid[i];
    return outOI;
  }

  // Zhang-Suen thinning
  bool changed = true;
  while(changed){
    changed = false;

    // Sub-iteration 1
    std::vector<int> toRemove;
    for(int row = 0; row < nrow; row++){
      for(int col = 0; col < ncol; col++){
        int idx = row * ncol + col;
        if(grid[idx] != 1) continue;
        if(isAnchor(idx)) continue;   // anchor cells are never removed

        // 8-neighbors: P2=N, P3=NE, P4=E, P5=SE, P6=S, P7=SW, P8=W, P9=NW
        int p2 = getVal(row-1, col);
        int p3 = getVal(row-1, col+1);
        int p4 = getVal(row, col+1);
        int p5 = getVal(row+1, col+1);
        int p6 = getVal(row+1, col);
        int p7 = getVal(row+1, col-1);
        int p8 = getVal(row, col-1);
        int p9 = getVal(row-1, col-1);

        int B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
        if(B < 2 || B > 6) continue;

        // transitions 0->1 in the ordered sequence P2..P9,P2
        int A = 0;
        A += (p2 == 0 && p3 == 1) ? 1 : 0;
        A += (p3 == 0 && p4 == 1) ? 1 : 0;
        A += (p4 == 0 && p5 == 1) ? 1 : 0;
        A += (p5 == 0 && p6 == 1) ? 1 : 0;
        A += (p6 == 0 && p7 == 1) ? 1 : 0;
        A += (p7 == 0 && p8 == 1) ? 1 : 0;
        A += (p8 == 0 && p9 == 1) ? 1 : 0;
        A += (p9 == 0 && p2 == 1) ? 1 : 0;
        if(A != 1) continue;

        // conditions for sub-iteration 1
        if(p2 * p4 * p6 != 0) continue;
        if(p4 * p6 * p8 != 0) continue;

        toRemove.push_back(idx);
      }
    }
    for(int idx : toRemove){
      grid[idx] = 0;
      changed = true;
    }

    // Sub-iteration 2
    toRemove.clear();
    for(int row = 0; row < nrow; row++){
      for(int col = 0; col < ncol; col++){
        int idx = row * ncol + col;
        if(grid[idx] != 1) continue;
        if(isAnchor(idx)) continue;   // anchor cells are never removed

        int p2 = getVal(row-1, col);
        int p3 = getVal(row-1, col+1);
        int p4 = getVal(row, col+1);
        int p5 = getVal(row+1, col+1);
        int p6 = getVal(row+1, col);
        int p7 = getVal(row+1, col-1);
        int p8 = getVal(row, col-1);
        int p9 = getVal(row-1, col-1);

        int B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
        if(B < 2 || B > 6) continue;

        int A = 0;
        A += (p2 == 0 && p3 == 1) ? 1 : 0;
        A += (p3 == 0 && p4 == 1) ? 1 : 0;
        A += (p4 == 0 && p5 == 1) ? 1 : 0;
        A += (p5 == 0 && p6 == 1) ? 1 : 0;
        A += (p6 == 0 && p7 == 1) ? 1 : 0;
        A += (p7 == 0 && p8 == 1) ? 1 : 0;
        A += (p8 == 0 && p9 == 1) ? 1 : 0;
        A += (p9 == 0 && p2 == 1) ? 1 : 0;
        if(A != 1) continue;

        // conditions for sub-iteration 2
        if(p2 * p4 * p8 != 0) continue;
        if(p2 * p6 * p8 != 0) continue;

        toRemove.push_back(idx);
      }
    }
    for(int idx : toRemove){
      grid[idx] = 0;
      changed = true;
    }
  }

  // build output
  IntegerVector out(n, NA_INTEGER);
  for(int i = 0; i < n; i++){
    if(grid[i] == -1){
      out[i] = NA_INTEGER;
    } else {
      out[i] = grid[i];
    }
  }

  return out;
}
