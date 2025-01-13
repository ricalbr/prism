#include "include/prism/config.hpp"
#include "include/prism/simulation.hpp"
#include "include/prism/utils.hpp"
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h>

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
        return 1;
    }

    auto [sim_config, phot_dist] = parse_yaml(argv[1]);

    // Create dark count rates using logspace
    std::vector<double> dcr(sim_config.num_det);
    double log_max = std::log10(sim_config.dcr_max), log_min = std::log10(sim_config.dcr_min);
    double step = (log_max - log_min) / (sim_config.num_det - 1);
    for (int i = 0; i < sim_config.num_det; ++i) {
        dcr[i] = std::pow(10.0, log_min + i * step);
    }
    auto frequencies = get_clicks_array(sim_config.num_det, sim_config.eta, dcr, sim_config.xtk, sim_config.iterations,
                                        phot_dist, 12345);

    write_array(frequencies);

    return 0;
}
