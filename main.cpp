#include "include/prism/config.hpp"
#include "include/prism/matrix.hpp"
#include "include/prism/simulation.hpp"
#include "include/prism/solver.hpp"
#include "include/prism/utils.hpp"
#include "include/prism/xtalk.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <yaml-cpp/yaml.h>

using namespace Eigen;

void saveMatrix(const Eigen::MatrixXd &M, const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore: impossibile aprire il file " << filename << "\n";
        return;
    }

    Eigen::IOFormat format(Eigen::FullPrecision, Eigen::DontAlignCols, " ",
                           "\n");
    file << M.format(format);
    file.close();
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml>\n";
        return 1;
    }

    auto [sim_config, phot_dist] = parse_yaml(argv[1]);

    int rows = sim_config.rows, cols = sim_config.cols;
    int num_det = rows * cols;

    VectorXd exponents =
        VectorXd::LinSpaced(num_det, std::log10(sim_config.dcr_min),
                            std::log10(sim_config.dcr_max));
    VectorXd dcr =
        exponents.unaryExpr([](double exp) { return std::pow(10.0, exp); });

    // Compute SPAD matrix
    std::cout << "Computing SPAD matrix\t\t";
    MatrixXd V = MatrixXd::Zero(num_det + 1, num_det + 1);
    get_spad_matrix(V, num_det, sim_config.eta, dcr.mean(), sim_config.xtk);
    float on_border =
        2.0 * static_cast<float>(sim_config.rows + (sim_config.cols - 1)) /
        static_cast<float>(sim_config.rows * sim_config.cols);
    float avg_nn = 3 * on_border + 4 * (1 - on_border);
    std::cout << "avg_nn " << avg_nn << std::endl;
    float xtk = 1 - std::pow(1 - sim_config.xtk, avg_nn);
    MatrixXd XT = x_matrix_gallego(xtk, num_det, 4);
    MatrixXd VG = XT * V;
    std::cout << "[DONE]\n";

    // Simulate clicks
    std::cout << "Simulating clicks\t\t";
    uint64_t seed = 123456789; // Seed for RNG
    auto start_time = std::chrono::high_resolution_clock::now();

    VectorXd c =
        get_clicks_array(rows, cols, sim_config.eta, dcr, sim_config.xtk,
                         sim_config.iterations, phot_dist, seed);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;
    std::cout << "[DONE]\n";
    std::cout << "Elapsed time: " << elapsed_seconds.count() << " s\n";
    // write_array(c, "freq.txt");

    // Statistics retrieval (EME)
    std::cout << "Retrieving photon statistics...\t" << std::endl;
    VectorXd p_out = EMESolver(VG, c, sim_config.alpha);
    write_array(p_out, "p_out.txt");

    return 0;
}
