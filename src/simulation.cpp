#include "../include/prism/simulation.hpp"
#include "../include/prism/random.hpp"
#include <Eigen/Dense>
#include <omp.h>
#include <queue>

#include <iostream>
#include <unordered_set>

// Funzione di hashing personalizzata per std::pair<int, int>
struct hash_pair {
    size_t operator()(const std::pair<int, int> &p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

using namespace Eigen;

int simulate(int rows, int cols, int num_ph, double eta,
             const Eigen::VectorXd &dcr, double xtk, Xoshiro256PlusPlus &rng) {
    std::unordered_set<std::pair<int, int>, hash_pair> detector_mat;
    std::queue<std::pair<int, int>> clicked;

    // Dark counts
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            if (rng.next_double() < dcr[c + r * rows]) {
                detector_mat.insert({r, c});
                clicked.push({r, c});
            }
        }
    }

    // Photon counts
    for (int i = 0; i < num_ph; ++i) {
        if (rng.next_double() < eta) {
            int r = rng.next() % rows, c = rng.next() % cols;
            if (detector_mat.insert({r, c})
                    .second) { // .second is true only if the element is
                               // inserted in the set (for the first time)
                clicked.push({r, c});
            }
        }
    }

    // Crosstalk
    if (xtk > 0.0) {
        static const std::vector<std::pair<int, int>> directions = {
            {-1, 0}, {1, 0}, {0, -1}, {0, 1}};

        auto is_valid = [&](int a, int b) {
            return (0 <= a) && (a < rows) && (0 <= b) && (b < cols);
        };

        std::queue<std::pair<int, int>> nn;
        do {
            while (!clicked.empty()) {
                // current clicked pixel
                auto [r, c] = clicked.front();
                clicked.pop();

                // find nearest neighbors
                for (const auto &[drow, dcol] : directions) {
                    int nr = r + drow, nc = c + dcol;
                    if (is_valid(nr, nc) &&
                        detector_mat.find({nr, nc}) == detector_mat.end()) {
                        nn.push({nr, nc});
                    }
                }
            }

            // MC step to evaluate clicks
            while (!nn.empty()) {
                auto [nr, nc] = nn.front();
                nn.pop();

                if (rng.next_double() < xtk) {
                    if (detector_mat.insert({nr, nc}).second) {
                        clicked.push({nr, nc});
                    }
                }
            }
        } while (!clicked.empty());
    }
    return detector_mat.size();
}
VectorXd get_clicks_array(
    int rows, int cols, double eta, const Eigen::VectorXd &dcr, double xtk,
    int iterations,
    const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution,
    uint64_t seed) {

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
                std::cout << "Saturating... measured clicks: " << sum_clicks
                          << "\n";
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
