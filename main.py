from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scipy


def main():
    D = 64

    # INPUT PHOTON STATISTICS
    ph_stat = scipy.stats.rv_discrete(values=([20], [1]))
    # ph_stat = scipy.stats.rv_discrete(values=([1, 7, 8, 9, 10], [0.1, 0.2, 0.4, 0.2, 0.1]))
    # ph_stat = scipy.stats.poisson(mu=20)
    # ph_stat = scipy.stats.boltzmann(lambda_=0.2, N=6)

    # FIDELITY
    p_n = np.loadtxt('p.txt')
    p_space = np.array(range(D + 1))
    fid = 1 - np.linalg.norm(p_n - ph_stat.pmf(p_space))
    print(f"Fidelity of the reconstruction: {fid:.3f}.")

    # PLOT
    plt.figure()

    plt.bar(p_space, p_n, alpha=0.6, color="red", label="Analytic + EME")
    plt.plot(ph_stat.pmf(p_space), "o--k", label="Theoretical p(n)")

    plt.legend()
    plt.xlabel("n")
    plt.ylabel("p(n)")
    plt.ylim([0, 1])
    plt.show()


if __name__ == "__main__":
    main()
