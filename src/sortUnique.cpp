#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

//' Extract and sort unique elements (c++)
//'
//' C++ function that extract and returns a sorted vector of numeric values.
//' @param x [numeric(.)][numeric]\cr the vector of numeric values
//' @family extractor functions
//' @return a numeric vector of the unique and sorted values.
//' @export
// [[Rcpp::export]]
NumericVector sortUniqueCpp(NumericVector x) {
  // [[Rcpp::plugins(cpp11)]]
  std::unordered_set<double> seen;
  int n = x.size();
  std::vector<double> out;

  for (int i = 0; i < n; ++i) {
    if (seen.insert(x[i]).second) out.push_back(x[i]);
  }

  std::sort(out.begin(), out.end());

  return wrap(out);
}
