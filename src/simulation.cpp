#include "../include/prism/simulation.hpp"
#include "../include/prism/random.hpp"
#include <Eigen/Dense>
#include <bitset>
#include <omp.h>

#include <iostream>

using namespace Eigen;

const int MAX_DET = 256;

int simulate(int num_det, int num_ph, double eta, const Eigen::VectorXd &dcr, double xtk, Xoshiro256PlusPlus &rng) {
    std::bitset<MAX_DET> click_space;

    // Dark counts
    for (int i = 0; i < num_det; ++i) {
        if (rng.next_double() < dcr[i]) {
            click_space.set(i);
        }
    }

    // Photon counts
    for (int i = 0; i < num_ph; ++i) {
        if (rng.next_double() < eta) {
            click_space.set(rng.next() % num_det);
        }
    }

    // Crosstalk
    if (xtk > 0.0) {
        std::bitset<MAX_DET> nbors = (click_space << 1) | (click_space >> 1);
        nbors &= ~click_space;

        for (int i = 0; i < num_det; ++i) {
            if (nbors[i] && rng.next_double() < xtk) {
                click_space.set(i);
            }
        }
    }

    return click_space.count();
}

VectorXd get_clicks_array(int num_det, double eta, const Eigen::VectorXd &dcr, double xtk, int iterations,
                          const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution, uint64_t seed) {
    VectorXd frequencies = VectorXd::Zero(num_det + 1);

#pragma omp parallel
    {
        VectorXd local_frequencies = VectorXd::Zero(num_det + 1);
        Xoshiro256PlusPlus rng(seed + omp_get_thread_num());

#pragma omp for schedule(static, 100)
        for (int i = 0; i < iterations; ++i) {
            int num_ph = photon_distribution(rng);
            int sum_clicks = simulate(num_det, num_ph, eta, dcr, xtk, rng);
            if (sum_clicks < local_frequencies.size()) { // Ensure valid index
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
