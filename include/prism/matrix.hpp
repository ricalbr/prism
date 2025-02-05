#ifndef PRISM_MATRIX_HPP
#define PRISM_MATRIX_HPP

#include <Eigen/Dense>

long long comb(int n, int k);

template <typename... Args> std::string make_key(Args... args);

double analytic(int D, double eta, double dcr, int N, int L);

void get_spad_matrix(Eigen::MatrixXd &A, int num_det, double eta, double dcr,
                     double xtk);

#endif // PRISM_MATRIX_HPP
