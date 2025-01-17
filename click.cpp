#include <bits/stdc++.h>
#include <omp.h>
#include <yaml-cpp/yaml.h>

using namespace std;

// Adjust this constant to the maximum number of detectors
const int MAX_DET = 256;

// xoshiro256++ random number generator implementation
struct Xoshiro256PlusPlus {
    uint64_t state[4];

    Xoshiro256PlusPlus(uint64_t seed) {
        for (int i = 0; i < 4; ++i) {
            state[i] = seed = seed * 6364136223846793005ULL + 1;
        }
    }

    uint64_t rotl(const uint64_t x, int k) const { return (x << k) | (x >> (64 - k)); }

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

// Poissonian photon generator using xoshiro256++
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

// Discrete photon generator using a probability distribution
int generate_photons_discrete(const vector<double> &probabilities, Xoshiro256PlusPlus &rng) {
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

// Simulates the click space for a detector array
int simulate(int num_det, int num_ph, double eta, const vector<double> &dcr, double xtk, Xoshiro256PlusPlus &rng) {
    bitset<MAX_DET> click_space;

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
        bitset<MAX_DET> nbors = (click_space << 1) | (click_space >> 1);
        nbors &= ~click_space;

        for (int i = 0; i < num_det; ++i) {
            if (nbors[i] && rng.next_double() < xtk) {
                click_space.set(i);
            }
        }
    }

    // Count and return the number of detectors that clicked
    return click_space.count();
}

// Parallelized get_clicks_array function
template <typename PhotonGenerator>
vector<double> parallel_get_clicks_array(int num_det, double eta, const vector<double> &dcr, double xtk, int iterations,
                                         PhotonGenerator photon_distribution, uint64_t seed) {
    vector<double> frequencies(num_det + 1, 0.0);

#pragma omp parallel
    {
        // Thread-local storage for frequencies
        vector<double> local_frequencies(num_det + 1, 0.0);
        Xoshiro256PlusPlus rng(seed + omp_get_thread_num());

#pragma omp for schedule(static, 10)
        for (int i = 0; i < iterations; ++i) {
            int num_ph = photon_distribution(rng);
            int sum_clicks = simulate(num_det, num_ph, eta, dcr, xtk, rng);
            if (sum_clicks <= num_det) {
                local_frequencies[sum_clicks] += 1.0;
            }
        }

#pragma omp critical
        {
            for (int i = 0; i <= num_det; ++i) {
                frequencies[i] += local_frequencies[i];
            }
        }
    }

    for (auto &freq : frequencies) {
        freq /= iterations;
    }

    return frequencies;
}

void write_array(const vector<double> &arr, const string &filename = "c.txt") {
    ofstream file(filename);
    if (file.is_open()) {
        for (const auto &e : arr) {
            file << e << "\n";
        }
        file.close();
    } else {
        cerr << "Unable to open file\n";
    }
}

struct SimulationConfig {
    int num_det;
    double eta;
    double xtk;
    double iterations;
    double dcr_min;
    double dcr_max;
};

SimulationConfig parse_yaml(const YAML::Node &config) {
    SimulationConfig sim_config;
    sim_config.num_det = config["num_det"].as<int>();
    sim_config.eta = config["eta"].as<double>();
    sim_config.xtk = config["xtk"].as<double>();
    sim_config.iterations = config["iterations"].as<double>();
    sim_config.dcr_min = config["dcr_min"].as<double>();
    sim_config.dcr_max = config["dcr_max"].as<double>();
    return sim_config;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <config.yaml>\n";
        return 1;
    }

    // Parse YAML configuration
    YAML::Node config = YAML::LoadFile(argv[1]);
    SimulationConfig sim_config = parse_yaml(config);
    string dist_type = config["photon_distribution"]["type"].as<string>();

    // function per gestire la generazione dei fotoni
    function<int(Xoshiro256PlusPlus &)> phot_dist;

    if (dist_type == "poisson") {
        double mean = config["photon_distribution"]["mean"].as<double>();
        phot_dist = [mean](Xoshiro256PlusPlus &rng) { return generate_photons_poisson(mean, rng); };
    } else if (dist_type == "discrete") {
        vector<double> probs = config["photon_distribution"]["probabilities"].as<vector<double>>();
        phot_dist = [probs](Xoshiro256PlusPlus &rng) { return generate_photons_discrete(probs, rng); };
    } else {
        throw invalid_argument("Unsupported photon distribution type");
    }

    // Clock timing
    double t1 = omp_get_wtime();

    // Create dark count rates using logspace
    vector<double> dcr(sim_config.num_det);
    double log_min = log10(sim_config.dcr_min), log_max = log10(sim_config.dcr_max);
    double step = (log_max - log_min) / (sim_config.num_det - 1);
    for (int i = 0; i < sim_config.num_det; ++i) {
        dcr[i] = pow(10.0, log_min + i * step);
    }

    // Seed for RNG
    uint64_t seed = 123456789;

    // Run the simulation using template-based photon generator
    auto frequencies = parallel_get_clicks_array(sim_config.num_det, sim_config.eta, dcr, sim_config.xtk,
                                                 sim_config.iterations, phot_dist, seed);

    // Write results
    write_array(frequencies);

    double t2 = omp_get_wtime();
    printf("Elapsed time is %.2lf seconds.\n\n", t2 - t1);

    return 0;
}
