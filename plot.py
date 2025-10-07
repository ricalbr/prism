from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
from matplotlib import rc

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Nimbus Sans'], 'size': 14})
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


def get_fidelity(filename, distribution, cwd, plot=False):
    p_n = np.loadtxt(filename)
    p_space = np.array(range(len(p_n)))
    p_stat = distribution.pmf(p_space)
    # fid = 1 - np.linalg.norm(p_n - p_stat)
    # print(f'{fid=}')
    fid = np.sum([np.sqrt(p_n * p_stat)]) ** 2
    print(f'{fid=}')
    hell_distance = 0.5 * np.sum(np.square(np.sqrt(p_n) - np.sqrt(p_stat)))
    print(f'{hell_distance=}')

    # EXPORT
    np.savetxt(cwd / 'predicted.txt', p_n)
    np.savetxt(cwd / 'theoretical.txt', p_stat)

    # PLOT
    plt.figure(1)
    plt.bar(p_space, p_n, alpha=0.4, color="blue", label="Retrieved distribution")
    plt.plot(p_stat, "ok", label="Theoretical distribution")

    plt.legend()
    plt.xlabel("n")
    plt.ylabel("p(n)")
    plt.title(f'Reconstruction fidelity: {fid*100:.4f}%')
    # plt.ylim([0, 1])

    if plot:
        plt.show()
        plt.savefig(cwd / 'plot.png')
    return fid


def main():

    # Load data from config.yaml
    yaml_data = load_yaml("config.yaml")
    photon_distribution = yaml_data.get("photon_distribution", {})
    distribution = create_distribution(photon_distribution)
    r, c, e, x = yaml_data.get("num_row"), yaml_data.get("num_col"), yaml_data.get("eta"), yaml_data.get("xtk")
    filename = f'M{r}x{c}_PDE{e*100}_XT{x*100}'
    if photon_distribution.get('type') == 'poisson':
        filename += f'_STAT_{photon_distribution.get("type")}_{photon_distribution.get("mean")}'
    else:
        stat = '-'.join([str(p) for p in photon_distribution.get("probabilities")])
        filename += f'_STAT_{photon_distribution.get("type")}_{stat}'

    cwd = pathlib.Path('.') / 'simulations' / filename
    pathlib.Path(cwd).mkdir(parents=True, exist_ok=True)

    try:
        fid = get_fidelity('p_out.txt', distribution, cwd)
        print(f"Fidelity of the reconstruction: {fid:.3f}.")
    except FileNotFoundError:
        print('File not found.')


if __name__ == "__main__":
    main()
