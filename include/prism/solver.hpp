#ifndef PRSIM_SOLVER_HPP
#define PRSIM_SOLVER_HPP

#include <Eigen/Dense>

using namespace Eigen;

VectorXd EMESolver(MatrixXd &A, VectorXd &b, double alpha = 1e-3, int iterations = (int)1e8, double eps = 1e-15,
                   int max_stagnation = (int)1e5);

#endif // !PRSIM_SOLVER_HPP
