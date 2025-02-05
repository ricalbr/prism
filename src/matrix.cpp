#include "../include/prism/matrix.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace Eigen;

long long comb(int n, int k) {
    if (k > n || k < 0)
        return 0;
    if (k == 0 || k == n)
        return 1;
    std::vector<long long> c(k + 1, 0);
    c[0] = 1;
    for (int i = 1; i <= n; ++i) {
        for (int j = std::min(i, k); j > 0; --j) {
            c[j] = c[j] + c[j - 1];
        }
    }
    return c[k];
}

std::unordered_map<std::string, double> memo;
template <typename... Args> std::string make_key(Args... args) {
    std::ostringstream oss;
    ((oss << args << ","), ...);
    return oss.str();
}

double P(int D, double eta, double dcr, int N, int C, int K) {
    if (N < 0 || C < 0 || C > N)
        return 0.0;
    if (C == 0 && N == 0) {
        return comb(D, K) * pow(dcr, K) * pow(1 - dcr, D - K);
    }

    std::string key = make_key(D, eta, dcr, N, C, K);
    if (memo.find(key) != memo.end()) {
        return memo[key];
    }

    double term1 =
        ((1 - eta) + eta * (K + C) / D) * P(D, eta, dcr, N - 1, C, K);
    double term2 =
        (eta * (D - (K + C - 1)) / D) * P(D, eta, dcr, N - 1, C - 1, K);
    double result = term1 + term2;

    memo[key] = result;
    return result;
}

double analytic(int D, double eta, double dcr, int N, int L) {
    double sum = 0.0;
    for (int i = 0; i <= L; ++i) {
        sum += P(D, eta, dcr, N, i, L - i);
    }
    return sum;
}

void get_spad_matrix(MatrixXd &A, int num_det, double eta, double dcr,
                     double xtk) {

    for (int n = 0; n <= num_det; ++n) {
        for (int m = 0; m <= num_det; ++m) {
            A(n, m) = analytic(num_det, eta, dcr, m, n);
        }
    }
}
