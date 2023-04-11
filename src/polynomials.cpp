#include <Rcpp.h>
using namespace Rcpp;

//' Evaluate Monomials.
//'
//' @param x A n x d matrix. The n input vectors, each of dimension d.
//' @param degrees A p x d matrix. The degrees of the p output monomial for each of the d dimensions.
//' @returns A n x p matrix. The values of the p monomials for each of the n inputs vectors.
//' @export
// [[Rcpp::export]]
NumericMatrix evaluateMonomials(NumericMatrix x, IntegerMatrix degrees) {

  int n = x.nrow();
  int d = x.ncol();
  int p = degrees.nrow();
  NumericMatrix out(n, p);
  double v;

  List powers(d);
  for (int j = 0; j < d; ++j) {
    int maxDegree = 0;
    for (int k = 0; k < p; ++k) {
      if (degrees(k, j) > maxDegree) {
        maxDegree = degrees(k, j);
      }
    }
    NumericMatrix powers1(n, maxDegree+1);
    for (int i = 0; i < n; ++i) {
      powers1(i, 0) = 1;
    }
    for (int l = 1; l <= maxDegree; ++l) {
      for (int i = 0; i < n; ++i) {
        powers1(i, l) = powers1(i, l-1) * x(i, j);
      }
    }
    powers[j] = powers1;
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < p; ++k) {
      v = 1;
      for (int j = 0; j < d; ++j) {
        NumericMatrix pwrs = powers[j];
        v *= pwrs(i, degrees(k, j));
      }
      out(i, k) = v;
    }
  }

  return out;
}
