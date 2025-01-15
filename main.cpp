#include "include/prism/config.hpp"
#include "include/prism/matrix.hpp"
#include "include/prism/simulation.hpp"
#include "include/prism/solver.hpp"
#include "include/prism/utils.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <yaml-cpp/yaml.h>

using namespace Eigen;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
        return 1;
    }

    auto [sim_config, phot_dist] = parse_yaml(argv[1]);

    // VectorXd dcr = VectorXd::LinSpaced(sim_config.num_det,
    // std::log(sim_config.dcr_min), std::log10(sim_config.dcr_max));
    VectorXd exponents =
        VectorXd::LinSpaced(sim_config.num_det, std::log10(sim_config.dcr_min),
                            std::log10(sim_config.dcr_max));
    VectorXd dcr =
        exponents.unaryExpr([](double exp) { return std::pow(10.0, exp); });

    // Compute SPAD matrix
    std::cout << "Computing SPAD matrix\t";
    MatrixXd V = MatrixXd::Zero(sim_config.num_det + 1, sim_config.num_det + 1);
    get_spad_matrix(V, sim_config.num_det, sim_config.eta, dcr.mean(),
                    sim_config.xtk);
    std::cout << "[DONE]\n";

    // // Simulate clicks
    std::cout << "Simulating clicks\t";
    uint64_t seed = 123456789; // Seed for RNG
    VectorXd c = get_clicks_array(sim_config.num_det, sim_config.eta, dcr,
                                  sim_config.xtk, sim_config.iterations,
                                  phot_dist, seed);
    std::cout << "[DONE]\n";
    // write_array(c, "freq.txt");

    // // Statistics retrieval (EME)
    std::cout << "Retrieving photon statistics...\t";
    VectorXd p = EMESolver(V, c, sim_config.alpha);
    write_array(p, "p.txt");

    return 0;
}
