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
    int mask[D] = {};
    std::bitset<D> click_space_old = {};

    Xoshiro256PlusPlus rng(12345);

    for (int i = 0; i < D; ++i) {
        if (rng.next_double() < PROBABILITY) {
            click_space.set(i);
        }
        mask[i] = click_space[i];
    }

    std::cout << "init\t\t";
    for (int i = 0; i < D; i++) {
        if (click_space[i]) {
            std::cout << "*";
        } else {
            std::cout << "-";
        }
    }
    std::cout << std::endl;

    while (true) {
        int nbors[D];

        // Calcolo dei vicini (find_nbors integrata)
        for (int i = 0; i < D; ++i) {
            int left = (i > 0) ? (click_space ^ click_space_old)[i - 1] : 0;
            int right =
                (i < D - 1) ? (click_space ^ click_space_old)[i + 1] : 0;
            nbors[i] = left + right;
        }

        // Aggiorna il mascheramento
        for (int i = 0; i < D; ++i) {
            mask[i] = (mask[i] == 0) ? nbors[i] : 0;
        }

        std::cout << "nbors\t\t";
        for (int i = 0; i < D; i++) {
            if (mask[i] == 2) {
                std::cout << "x";
            } else if (mask[i] == 1) {
                std::cout << "v";
            } else {
                std::cout << " ";
            }
        }
        std::cout << std::endl;

        click_space_old = click_space;
        for (int i = 0; i < D; ++i) {
            for (int j = 0; j < mask[i]; ++j) {
                if (rng.next_double() < XTK) {
                    click_space.set(i);
                }
            }
        }

        std::cout << "xtalk\t\t";
        for (int i = 0; i < D; i++) {
            if (click_space_old[i] - click_space[i]) {
                std::cout << "*";
            } else if (click_space[i]) {
                std::cout << "*";
            } else {
                std::cout << "-";
            }
        }
        std::cout << std::endl;

        if ((click_space ^ click_space_old).count() == 0) {
            break;
        }

        // Aggiorna il mascheramento considerando il nuovo stato
        for (int i = 0; i < D; ++i) {
            mask[i] = click_space_old[i] + click_space[i];
        }

        // std::cout << "nbors\t\t";
        // for (int i = 0; i < D; i++) {
        //     if (mask[i] == 2) {
        //         std::cout << "x";
        //     } else if (mask[i] == 1) {
        //         std::cout << "v";
        //     } else {
        //         std::cout << " ";
        //     }
        // }
        // std::cout << std::endl;
        // std::cout << "\n\nnew state\t";
        // for (int i = 0; i < D; i++) {
        //     if (click_space[i]) {
        //         std::cout << "1";
        //     } else {
        //         std::cout << "-";
        //     }
        // }
        // std::cout << std::endl;
    }
    return 0;
}
