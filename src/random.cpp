#include "../include/prism/random.hpp"
#include <cmath>
#include <cstdint>

Xoshiro256PlusPlus::Xoshiro256PlusPlus(uint64_t seed) {
    for (int i = 0; i < 4; ++i) {
        state[i] = seed = seed * 6364136223846793005ULL + 1;
    }
}

uint64_t Xoshiro256PlusPlus::rotl(const uint_fast64_t x, int k) const { return (x << k) | (x >> (64 - k)); }

uint_fast64_t Xoshiro256PlusPlus::next() {
    const uint64_t result = rotl(state[0] + state[3], 23) + state[0];

    const uint64_t t = state[1] << 17;

    state[2] ^= state[0];
    state[3] ^= state[1];
    state[1] ^= state[2];
    state[0] ^= state[3];

    state[2] ^= t;

    state[3] = rotl(state[3], 45);

    return result;
}

double Xoshiro256PlusPlus::next_double() { return (next() >> 11) * (1.0 / (1ULL << 53)); }

int generate_photons_poisson(double mean, Xoshiro256PlusPlus &rng) {
    double L = exp(-mean);
    int k = 0;
    double p = 1.0;

    do {
        ++k;
        p *= rng.next_double();
    } while (p > L);

    return k - 1;
}

int generate_photons_discrete(const std::vector<double> &probabilities, Xoshiro256PlusPlus &rng) {
    double random_value = rng.next_double();
    double cumulative = 0.0;
    for (size_t i = 0; i < probabilities.size(); ++i) {
        cumulative += probabilities[i];
        if (random_value < cumulative) {
            return i;
        }
    }
    return probabilities.size() - 1;
}
