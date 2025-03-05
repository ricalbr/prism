#include "../include/prism/solver.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace Eigen;

VectorXd EMESolver(MatrixXd &A, VectorXd &b, double alpha, int iterations, double eps, int max_stagnation) {
    VectorXd EME = VectorXd::Zero(A.cols());
    VectorXd pn = VectorXd::Ones(A.cols()) / A.cols();

    int curr_iter = 0;
    int stagnation_counter = 0;
    double start_dist = 1.0;
    double last_dist = std::numeric_limits<double>::infinity();

    if (start_dist <= eps) {
        throw std::invalid_argument("Initial `dist` must be greater than `epsilon`.");
    }

    while (curr_iter <= iterations) {

        // EME computation
        ArrayXd mat_pn = A * pn;
        RowVectorXd tmp = (b.array() / mat_pn).matrix();
        tmp *= A;
        VectorXd EM = tmp.transpose();
        EM = EM.array() * pn.array();

        VectorXd log_pn = pn.array().log();
        VectorXd pn_log_pn = pn.dot(log_pn) * VectorXd::Ones(log_pn.size());

        VectorXd E = alpha * (log_pn - pn_log_pn).array() * pn.array();
        E = (E.array().isFinite()).select(E, 0.0); // Replace -inf with 0
        EME = EM - E;

        // Check convergence
        double dist = (EME - pn).norm();
        if (dist <= eps) {
            std::cout << "\nConverged to desired precision.\n";
            break;
        } else if (dist >= last_dist) {
            stagnation_counter++;
            if (stagnation_counter >= max_stagnation) {
                std::cout << "\nStagnation detected: the distance has not "
                             "improved for "
                          << max_stagnation << " iterations.\n";
                break;
            } else {
                stagnation_counter = 0;
            }
        }
        pn = EME;
        last_dist = std::min(last_dist, dist);
        curr_iter++;
    }
    return pn;
}
