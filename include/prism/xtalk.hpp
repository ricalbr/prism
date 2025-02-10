#ifndef PRISM_XTALK_HPP
#define PRISM_XTALK_HPP

#include <Eigen/Dense>

using namespace Eigen;

double p4_xt(double eps, int events);

double p2_xt(double eps, int events);

double p_m(double p, int k, int m, int nn);

MatrixXd x_matrix_simple(double eps, int num);

MatrixXd x_matrix_gallego(double eps, int num, int nn);

#endif // PRISM_XTALK_HPP
