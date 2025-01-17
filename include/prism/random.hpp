#ifndef PRISM_RANDOM_HPP
#define PRISM_RANDOM_HPP

#include <cstdint>
#include <vector>

class Xoshiro256PlusPlus {
  public:
    explicit Xoshiro256PlusPlus(uint_fast64_t seed);
    uint_fast64_t next();
    double next_double();

  private:
    uint_fast64_t state[4];
    uint_fast64_t rotl(const uint_fast64_t x, int k) const;
};

int generate_photons_poisson(double mean, Xoshiro256PlusPlus &rng);
int generate_photons_discrete(const std::vector<double> &probabilities, Xoshiro256PlusPlus &rng);

#endif // PRISM_RANDOM_HPP
