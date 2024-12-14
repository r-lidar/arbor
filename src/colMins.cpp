#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector colMins(NumericMatrix mat) {
  int ncol = mat.ncol();
  IntegerVector result(ncol, NA_INTEGER);

  for (int col = 0; col < ncol; ++col) {
    double minVal = R_PosInf;
    int minIndex = -1;
    bool allNA = true;

    for (int row = 0; row < mat.nrow(); ++row) {
      double value = mat(row, col);

      if (!NumericVector::is_na(value)) {
        allNA = false;
        if (value < minVal) {
          minVal = value;
          minIndex = row + 1; // Convert to 1-based index
        }
      }
    }

    if (!allNA) {
      result[col] = minIndex;
    }
  }

  return result;
}
