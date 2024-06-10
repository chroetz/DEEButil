#include <Rcpp.h>
using namespace Rcpp;

//' Get the index of the closest row.
//'
//' @param target A n x p matrix.
//' @param query A vector of length p.
//' @returns The smallest index of a row of target that has the smallest
//'   Euclidean distance to query.
//' @export
// [[Rcpp::export]]
int whichMinDist(NumericMatrix target, NumericVector query) {
  int nrow = target.nrow(), ncol = target.ncol();
  double minDist = INFINITY;
  double minIdx = -1;
  double dst;
  double v;
  for (int i = 0; i < nrow; ++i) {
    dst = 0;
    for (int j = 0; j < ncol; ++j) {
      v = target(i, j) - query(j);
      dst += v*v;
    }
    if (dst < minDist) {
      minDist = dst;
      minIdx = i;
    }
  }
  return minIdx+1;
}

//' @export
// [[Rcpp::export]]
NumericVector distToVec(NumericMatrix x, NumericVector y) {
  int n = x.nrow(), d = x.ncol();
  NumericVector out(n);
  double dst;
  double v;
  for (int i = 0; i < n; ++i) {
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = x(i, j) - y(j);
      dst += v*v;
    }
    out[i] = sqrt(dst);
  }
  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector distSqrToVec(NumericMatrix target, NumericVector query) {
  int n = target.nrow(), d = target.ncol();
  NumericVector out(n);
  double dst;
  double v;
  for (int i = 0; i < n; ++i) {
    dst = 0;
    for (int j = 0; j < d; ++j) {
      v = target(i, j) - query(j);
      dst += v*v;
    }
    out[i] = dst;
  }
  return out;
}

