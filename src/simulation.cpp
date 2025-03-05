#include "../include/prism/simulation.hpp"
#include "../include/prism/random.hpp"
#include <Eigen/Dense>
#include <omp.h>
#include <queue>
#include <set>

#include <iostream>

using namespace Eigen;

int simulate(int rows, int cols, int num_ph, double eta, const Eigen::VectorXd &dcr, double xtk,
             Xoshiro256PlusPlus &rng) {

    Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> detector_mat(rows, cols);
    detector_mat.setConstant(false);

    std::queue<std::pair<int, int>> clicked, nn;
    int num_det = rows * cols;

    // Dark counts
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            if (rng.next_double() < dcr[r + c * rows]) {
                detector_mat(r, c) = true;
                clicked.push({r, c});
            }
        }
    }

    // Photon counts
    for (int i = 0; i < num_ph; ++i) {
        if (rng.next_double() < eta) {
            int r = rng.next() % rows, c = rng.next() % cols;
            detector_mat(r, c) = true;
            clicked.push({r, c});
        }
    }

    // Crosstalk
    if (xtk > 0.0) {
        std::vector<std::pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

        auto is_valid = [&](int a, int b) { return (0 <= a) && (a < rows) && (0 <= b) && (b < cols); };

        do {
            // find nearest neighbors
            while (!clicked.empty()) {
                auto &[r, c] = clicked.front();
                clicked.pop();
                for (auto &[drow, dcol] : directions) {
                    int nr = r + drow, nc = c + dcol;
                    if (is_valid(nr, nc) && !detector_mat(nr, nc)) {
                        nn.push({nr, nc});
                    }
                }
            }
            // MC step for evaluating clicks
            std::set<std::pair<int, int>> s;
            while (!nn.empty()) {
                if ((double)rand() / RAND_MAX < xtk) {
                    s.insert(nn.front());
                    detector_mat(nn.front().first, nn.front().second) = true;
                }
                nn.pop();
            }

            // a set ensure the uniqueness of the clicked sites, a (r,c) site can only click once even though it may be
            // evaluated more than once (4 times at max)
            for (std::pair<int, int> elem : s) {
                clicked.push(elem);
            }
        } while (!clicked.empty());
    }
    return (detector_mat.cast<int>().array() == 1).count();
}

VectorXd get_clicks_array(int rows, int cols, double eta, const Eigen::VectorXd &dcr, double xtk, int iterations,
                          const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution, uint64_t seed) {

    int num_det = rows * cols;
    VectorXd frequencies = VectorXd::Zero(num_det + 1);

#pragma omp parallel
    {
        VectorXd local_frequencies = VectorXd::Zero(num_det + 1);
        Xoshiro256PlusPlus rng(seed + omp_get_thread_num());

#pragma omp for schedule(static, 100)
        for (int i = 0; i < iterations; ++i) {
            int num_ph = photon_distribution(rng);
            int sum_clicks = simulate(rows, cols, num_ph, eta, dcr, xtk, rng);
            if (sum_clicks < num_det) { // Ensure valid index
                local_frequencies[sum_clicks] += 1.0;
            } else {
                std::cout << "Saturating... measured clicks: " << sum_clicks << "\n";
                local_frequencies[local_frequencies.size() - 1] += 1.0;
            }
        }

#pragma omp critical
        {
            for (size_t i = 0; i < frequencies.size(); ++i) {
                frequencies[i] += local_frequencies[i];
            }
        }
    }

    for (auto &freq : frequencies) {
        freq /= iterations;
    }

    return frequencies;
}
