#ifndef PRISM_SIMULATION_HPP
#define PRISM_SIMULATION_HPP

#include "random.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <functional>

int simulate(int rows, int cols, int num_ph, double eta, const Eigen::VectorXd &dcr, double xtlk,
             Xoshiro256PlusPlus &rng);
Eigen::VectorXd get_clicks_array(int rows, int cols, double eta, const Eigen::VectorXd &dcr, double xtlk,
                                 int iterations, const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution,
                                 uint64_t seed);

#endif // PRISM_SIMULATION_HPP
