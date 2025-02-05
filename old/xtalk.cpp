#include <Eigen/Dense>
#include <bitset>
#include <iostream>

const int D = 40;               // Dimensione dell'array (modificabile)
const double PROBABILITY = 0.2; // Probabilità iniziale per l'array
const double XTK = 0.05;        // Probabilità di aggiornamento

// xoshiro256++ random number generator implementation
struct Xoshiro256PlusPlus {
    uint64_t state[4];

    Xoshiro256PlusPlus(uint64_t seed) {
        for (int i = 0; i < 4; ++i) {
            state[i] = seed = seed * 6364136223846793005ULL + 1;
        }
    }

    uint64_t rotl(const uint64_t x, int k) const {
        return (x << k) | (x >> (64 - k));
    }

    uint64_t next() {
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

    double next_double() { return (next() >> 11) * (1.0 / (1ULL << 53)); }
};

using namespace Eigen;

int main() {
    std::bitset<D> click_space;
    std::array<int, D> mask = {};
    std::array<int, D> nbors = {};
    std::bitset<D> click_space_old = {};

    Xoshiro256PlusPlus rng(12345);

    for (int i = 0; i < D; ++i) {
        if (rng.next_double() < PROBABILITY) {
            click_space.set(i);
        }
        mask[i] = click_space[i];
    }

    while (true) {

        for (int i = 0; i < D; ++i) {
            int left = (i > 0) ? (click_space ^ click_space_old)[i - 1] : 0;
            int right =
                (i < D - 1) ? (click_space ^ click_space_old)[i + 1] : 0;
            nbors[i] = left + right;
        }

        for (int i = 0; i < D; ++i) {
            mask[i] = (mask[i] == 0) ? nbors[i] : 0;
        }

        click_space_old = click_space;
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < mask[i]; ++j) {
                if (rng.next_double() < XTK) {
                    click_space.set(i);
                }
            }
        }

        if ((click_space ^ click_space_old).none()) {
            break;
        }

        for (int i = 0; i < D; ++i) {
            mask[i] = click_space_old[i] + click_space[i];
        }
    }
}
