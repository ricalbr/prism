#ifndef PRISM_XTALK_HPP
#define PRISM_XTALK_HPP

#include <Eigen/Dense>

using namespace Eigen;

double p2_xt(double eps, int events);

double p_m(double p, int k, int m);

MatrixXd x_matrix_simple(double eps, int num);

MatrixXd x_matrix_std(double eps, int num);

MatrixXd x_matrix_gallego(double eps, int num);

#endif // PRISM_XTALK_HPP
