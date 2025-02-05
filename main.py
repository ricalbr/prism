from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import yaml


def load_yaml(filename):
    with open(filename, "r") as file:
        return yaml.safe_load(file)


def create_distribution(photon_distribution):
    dist_type = photon_distribution.get("type", "poisson")  # Default: Poisson

    if dist_type == "poisson":
        mean = photon_distribution.get("mean", 1)  # Default: 1
        return stats.poisson(mean)

    elif dist_type == "discrete":
        probabilities = photon_distribution.get("probabilities", [1.0])  # Default: one state
        values = np.arange(len(probabilities))  # Discrete index states
        return stats.rv_discrete(name="custom", values=(values, probabilities))

    else:
        raise ValueError(f"Photon statistics '{dist_type}' not supported!")


def get_fidelity(filename, distribution, num_det, plot=False):
    p_n = np.loadtxt(filename)
    p_space = np.array(range(num_det + 1))
    p_stat = distribution.pmf(p_space)
    fid = 1 - np.linalg.norm(p_n - p_stat)

    # PLOT
    plt.figure(1)
    plt.bar(p_space, p_n, alpha=0.6, color="red", label="Analytic + EME")
    plt.plot(p_stat, "o--k", label="Theoretical p(n)")

    plt.legend()
    plt.xlabel("n")
    plt.ylabel("p(n)")
    plt.ylim([0, 1])
    if plot:
        plt.show()
    return fid


def main():

    # Load data from config.yaml
    yaml_data = load_yaml("config.yaml")
    num_det = yaml_data.get("num_det")
    photon_distribution = yaml_data.get("photon_distribution", {})
    distribution = create_distribution(photon_distribution)

    # Standard
    fid = get_fidelity('pS.txt', distribution, num_det)
    print(f"Fidelity of the reconstruction for standard method: {fid:.3f}.")

    # Gallego
    fid = get_fidelity('pG.txt', distribution, num_det)
    print(f"Fidelity of the reconstruction for Gallego method: {fid:.3f}.")


if __name__ == "__main__":
    main()
