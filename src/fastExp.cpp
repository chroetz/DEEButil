#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector fast_exp_neg_sq(NumericVector x) {
    int n = x.size();
    NumericVector result(n);

    // Predefined points for piece-wise linear approximation
    const double exp0 = 1.0;                 // exp(-0^2)
    const double exp0_5 = 0.778800783;       // exp(-0.5^2)
    const double exp1 = 0.367879441;         // exp(-1^2)
    const double exp1_5 = 0.105399224;       // exp(-1.5^2)
    const double exp2 = 0.018315639;         // exp(-2^2)
    const double exp2_5 = 0.001930454;       // exp(-2.5^2)
    const double exp3 = 0.0001234098;        // exp(-3^2)
    const double exp3_5 = 0.000011109;       // exp(-3.5^2)
    const double exp4 = 0.0000003354626;     // exp(-4^2)
    const double exp4_5 = 0.000000005184705; // exp(-4.5^2)
    const double exp5 = 0.000000001388794;   // exp(-5^2)

    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        double xi2 = xi * xi;

        if (xi2 <= 0.25) {
            // Linear approximation in [0, 0.5^2]
            result[i] = exp0 + (xi2 - 0.0) * (exp0_5 - exp0) / (0.25 - 0.0);
        } else if (xi2 <= 1.0) {
            // Linear approximation in [0.5^2, 1^2]
            result[i] = exp0_5 + (xi2 - 0.25) * (exp1 - exp0_5) / (1.0 - 0.25);
        } else if (xi2 <= 2.25) {
            // Linear approximation in [1^2, 1.5^2]
            result[i] = exp1 + (xi2 - 1.0) * (exp1_5 - exp1) / (2.25 - 1.0);
        } else if (xi2 <= 4.0) {
            // Linear approximation in [1.5^2, 2^2]
            result[i] = exp1_5 + (xi2 - 2.25) * (exp2 - exp1_5) / (4.0 - 2.25);
        } else if (xi2 <= 6.25) {
            // Linear approximation in [2^2, 2.5^2]
            result[i] = exp2 + (xi2 - 4.0) * (exp2_5 - exp2) / (6.25 - 4.0);
        } else if (xi2 <= 9.0) {
            // Linear approximation in [2.5^2, 3^2]
            result[i] = exp2_5 + (xi2 - 6.25) * (exp3 - exp2_5) / (9.0 - 6.25);
        } else if (xi2 <= 12.25) {
            // Linear approximation in [3^2, 3.5^2]
            result[i] = exp3 + (xi2 - 9.0) * (exp3_5 - exp3) / (12.25 - 9.0);
        } else if (xi2 <= 16.0) {
            // Linear approximation in [3.5^2, 4^2]
            result[i] = exp3_5 + (xi2 - 12.25) * (exp4 - exp3_5) / (16.0 - 12.25);
        } else if (xi2 <= 20.25) {
            // Linear approximation in [4^2, 4.5^2]
            result[i] = exp4 + (xi2 - 16.0) * (exp4_5 - exp4) / (20.25 - 16.0);
        } else if (xi2 <= 25.0) {
            // Linear approximation in [4.5^2, 5^2]
            result[i] = exp4_5 + (xi2 - 20.25) * (exp5 - exp4_5) / (25.0 - 20.25);
        } else {
            // For large xi, approximate as 0
            result[i] = 0.0;
        }
    }

    return result;
}
