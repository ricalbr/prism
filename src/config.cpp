#include "../include/prism/config.hpp"
#include "../include/prism/random.hpp"
#include <functional>
#include <yaml-cpp/yaml.h>

std::tuple<SimulationConfig, std::function<int(Xoshiro256PlusPlus &)>>
parse_yaml(const std::string &filename) {

    YAML::Node config = YAML::LoadFile(filename);
    SimulationConfig sim_config;
    sim_config.rows = config["num_row"].as<int>();
    sim_config.cols = config["num_col"].as<int>();
    sim_config.eta = config["eta"].as<double>();
    sim_config.xtk = config["xtk"].as<double>();
    sim_config.iterations = config["iterations"].as<double>();
    sim_config.dcr_min = config["dcr_min"].as<double>();
    sim_config.dcr_max = config["dcr_max"].as<double>();
    sim_config.alpha = config["alpha"].as<double>();

    std::function<int(Xoshiro256PlusPlus &)> phot_dist;

    std::string dist_type =
        config["photon_distribution"]["type"].as<std::string>();
    if (dist_type == "poisson") {
        double mean = config["photon_distribution"]["mean"].as<double>();
        phot_dist = [mean](Xoshiro256PlusPlus &rng) {
            return generate_photons_poisson(mean, rng);
        };
    } else if (dist_type == "discrete") {
        std::vector<double> probs =
            config["photon_distribution"]["probabilities"]
                .as<std::vector<double>>();
        phot_dist = [probs](Xoshiro256PlusPlus &rng) {
            return generate_photons_discrete(probs, rng);
        };
    } else {
        throw std::invalid_argument("Unsupported photon distribution type");
    }
    return {sim_config, phot_dist};
}
