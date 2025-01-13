#ifndef PRISM_SIMULATION_HPP
#define PRISM_SIMULATION_HPP

#include "random.hpp"
#include <cstdint>
#include <functional>
#include <vector>

int simulate(int num_det, int num_ph, double eta, const std::vector<double> &dcr, double xtlk, Xoshiro256PlusPlus &rng);
std::vector<double> get_clicks_array(int num_det, double eta, const std::vector<double> &dcr, double xtlk,
                                     int iterations,
                                     const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution,
                                     uint64_t seed);

#endif // PRISM_SIMULATION_HPP
