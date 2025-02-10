#include "../include/prism/xtalk.hpp"
#include "../include/prism/matrix.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <unordered_map>

template <typename... Args> std::string make_key(Args... args) {
    std::ostringstream oss;
    ((oss << args << ","), ...);
    return oss.str();
}

double p4_xt(double eps, int events) {
    double p = 1 - std::sqrt(1 - eps);
    double q = 1 - p;

    if (events < 1) {
        return 0;
    } else if (events == 1) {
        return 1 - eps;
    } else if (events == 2) {
        return 4 * p * std::pow(q, 6);
    } else if (events == 3) {
        return 18 * std::pow(p, 2) * std::pow(q, 8);
    } else if (events == 4) {
        return 4 * std::pow(p, 3) * std::pow(q, 8) * (1 + 3 * q + 18 * q * q);
    } else if (events == 5) {
        return 5 * std::pow(p, 4) * std::pow(q, 10) * (8 + 24 * q + 55 * q * q);
    } else {
        double p5 = p4_xt(eps, 5);
        double sum = 0.0;
        for (int i = 1; i <= 4; i++) {
            sum += p4_xt(eps, i);
        }
        return p5 * std::pow((1 - (p5 / (1 - sum))), events - 5);
    }
}

double p2_xt(double eps, int events) {
    double p = 1 - std::sqrt(1 - eps);
    double q = 1 - p;

    if (events < 1) {
        return 0;
    } else if (events == 1) {
        return 1 - eps;
    } else if (events == 2) {
        return 2 * p * std::pow(q, 3);
    } else if (events == 3) {
        return 5 * std::pow(p, 2) * std::pow(q, 4);
    } else if (events == 4) {
        return 14 * std::pow(p, 3) * std::pow(q, 5);
    } else if (events == 5) {
        return 42 * std::pow(p, 4) * std::pow(q, 6);
    } else {
        double p5 = p2_xt(eps, 5);
        double sum = 0.0;
        for (int i = 1; i <= 4; i++) {
            sum += p2_xt(eps, i);
        }
        return p5 * std::pow((1 - (p5 / (1 - sum))), events - 5);
    }
}

std::unordered_map<std::string, double> p_m_cache;
double p_m(double p, int k, int m, int nn) {
    std::string key = make_key(k, m, nn);
    if (p_m_cache.find(key) != p_m_cache.end())
        return p_m_cache[key];
    if (nn == 2) {
        if (m < 1) {
            double result = p2_xt(p, k);
            p_m_cache[key] = result;
            return result;
        } else {
            double result = 0;
            for (int i = 1; i <= k - m; ++i) {
                result += p_m(p, k - i, m - 1, nn) * p2_xt(p, i);
            }
            p_m_cache[key] = result;
            return result;
        }
    } else if (nn == 4) {
        if (m < 1) {
            double result = p4_xt(p, k);
            p_m_cache[key] = result;
            return result;
        } else {
            double result = 0;
            for (int i = 1; i <= k - m; ++i) {
                result += p_m(p, k - i, m - 1, nn) * p4_xt(p, i);
            }
            p_m_cache[key] = result;
            return result;
        }
    } else {
        return 0;
    }
}

MatrixXd x_matrix_gallego(double eps, int num, int nn) {
    MatrixXd mat = MatrixXd::Identity(num + 1, num + 1);
    if (eps <= 1e-9) {
        return mat;
    }

    for (int i = 0; i <= num; ++i) {
        for (int j = 0; j < num; ++j) {
            mat(i, j + 1) = p_m(eps, i, j, nn);
        }
    }
    return mat;
}

MatrixXd x_matrix_simple(double eps, int num) {
    MatrixXd mat = MatrixXd::Identity(num + 1, num + 1);
    if (eps <= 1e-9) {
        return mat;
    }

    for (int i = 0; i < num + 1; i++) {
        for (int j = 0; j < num + 1; j++) {
            mat(i, j) =
                comb(j, i - j) * pow(eps, i - j) * pow(1 - eps, 2 * j - i);
        }
    }
    return mat;
}
