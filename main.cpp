#include "include/prism/config.hpp"
#include "include/prism/matrix.hpp"
#include "include/prism/simulation.hpp"
#include "include/prism/solver.hpp"
#include "include/prism/utils.hpp"
#include "include/prism/xtalk.hpp"
#include <Eigen/Dense>
#include <chrono>
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

    VectorXd exponents =
        VectorXd::LinSpaced(sim_config.num_det, std::log10(sim_config.dcr_min),
                            std::log10(sim_config.dcr_max));
    VectorXd dcr =
        exponents.unaryExpr([](double exp) { return std::pow(10.0, exp); });

    // Compute SPAD matrix
    std::cout << "Computing SPAD matrix\t\t";
    MatrixXd V = MatrixXd::Zero(sim_config.num_det + 1, sim_config.num_det + 1);
    get_spad_matrix(V, sim_config.num_det, sim_config.eta, dcr.mean(),
                    sim_config.xtk);
    MatrixXd XG = x_matrix_gallego(sim_config.xtk, sim_config.num_det);
    MatrixXd XS = x_matrix_simple(sim_config.xtk, sim_config.num_det);
    MatrixXd VG = XG * V;
    MatrixXd VS = XS * V;
    std::cout << "[DONE]\n";

    // // Simulate clicks
    std::cout << "Simulating clicks\t\t";

    uint64_t seed = 123456789; // Seed for RNG
    auto start_time = std::chrono::high_resolution_clock::now();
    VectorXd c = get_clicks_array(sim_config.num_det, sim_config.eta, dcr,
                                  sim_config.xtk, sim_config.iterations,
                                  phot_dist, seed);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    std::cout << "[DONE]\n";
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " s\n";
    // write_array(c, "freq.txt");

    // Statistics retrieval (EME)
    std::cout << "Retrieving photon statistics...\t" << std::endl;

    std::cout << "\t- Standard method...\t";
    VectorXd pS = EMESolver(VS, c, sim_config.alpha);
    write_array(pS, "pS.txt");

    std::cout << "\t- Gallego method...\t";
    VectorXd pG = EMESolver(VG, c, sim_config.alpha);
    write_array(pG, "pG.txt");

    return 0;
}
