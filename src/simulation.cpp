#include "../include/prism/simulation.hpp"
#include "../include/prism/random.hpp"
#include <Eigen/Dense>
#include <bitset>
#include <omp.h>

#include <iostream>

using namespace Eigen;

const int MAX_DET = 256;

int simulate(int num_det, int num_ph, double eta, const Eigen::VectorXd &dcr,
             double xtk, Xoshiro256PlusPlus &rng) {
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
        std::bitset<MAX_DET> click_space_old = click_space;
        std::vector<int> mask(num_det, 0);
        std::vector<int> nbors(num_det, 0);
        for (int i = 0; i < num_det; ++i) {
            mask[i] = click_space[i];
        }
        while (true) {
            for (int i = 0; i < num_det; i++) {
                int left = 0, right = 0;
                if (i > 0) {
                    left = (click_space ^ click_space_old)[i - 1];
                }
                if (i < num_det - 1) {
                    right = (click_space ^ click_space_old)[i + 1];
                }
                nbors[i] = left + right;
            }

            for (int i = 0; i < num_det; i++) {
                if (mask[i] == 0) {
                    mask[i] = nbors[i];
                } else {
                    mask[i] = 0;
                }
            }

            click_space_old = click_space;
            for (int i = 0; i < num_det; i++) {
                for (int j = 0; j < mask[i]; j++) {
                    if (rng.next_double() < xtk) {
                        click_space.set(i);
                    }
                }
            }

            if ((click_space ^ click_space_old).none()) {
                break;
            }

            for (int i = 0; i < num_det; i++) {
                mask[i] = click_space_old[i] + click_space[i];
            }
        }
    }
    return click_space.count();
}

VectorXd get_clicks_array(
    int num_det, double eta, const Eigen::VectorXd &dcr, double xtk,
    int iterations,
    const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution,
    uint64_t seed) {
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
