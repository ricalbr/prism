#ifndef PRISM_CONFIG_HPP
#define PRISM_CONFIG_HPP

#include "random.hpp"
#include <functional>
#include <string>
#include <tuple>

struct SimulationConfig {
    int rows;
    int cols;
    double eta;
    double xtk;
    int iterations;
    double dcr_min;
    double dcr_max;
    double alpha;
};
std::tuple<SimulationConfig, std::function<int(Xoshiro256PlusPlus &)>> parse_yaml(const std::string &filename);

#endif // PRISM_CONFIG_HPP
