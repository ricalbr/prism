#include "../include/prism/simulation.hpp"
#include "../include/prism/random.hpp"
#include <bitset>
#include <omp.h>

const int MAX_DET = 256;

int simulate(int num_det, int num_ph, double eta, const std::vector<double> &dcr, double xtk, Xoshiro256PlusPlus &rng) {
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

std::vector<double> get_clicks_array(int num_det, double eta, const std::vector<double> &dcr, double xtk,
                                     int iterations,
                                     const std::function<int(Xoshiro256PlusPlus &)> &photon_distribution,
                                     uint64_t seed) {
    std::vector<double> frequencies(num_det + 1, 0.0);

#pragma omp parallel
    {
        std::vector<double> local_frequencies(num_det + 1, 0.0);
        Xoshiro256PlusPlus rng(seed + omp_get_thread_num());

#pragma omp for schedule(static, 100)
        for (int i = 0; i < iterations; ++i) {
            int num_ph = photon_distribution(rng);
            int sum_clicks = simulate(num_det, num_ph, eta, dcr, xtk, rng);
            if (sum_clicks <= num_det) {
                local_frequencies[sum_clicks] += 1.0;
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
