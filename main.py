from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scipy

import prism
from prism import EMEsolver


def main():
    # SPAD ARRAY/MATRIX PROPERTIES
    D = 64
    # mean_dcr, dcr = get_dcr_array(D=D, dcr_min=1e-3, dcr_max=1e-2)

    # INPUT PHOTON STATISTICS
    # ph_stat = scipy.stats.rv_discrete(values=([20], [1]))
    # ph_stat = scipy.stats.rv_discrete(values=([1], [1]))
    ph_stat = scipy.stats.rv_discrete(values=([1, 7, 8, 9, 10], [0.1, 0.2, 0.4, 0.2, 0.1]))
    # ph_stat = scipy.stats.poisson(mu=40)
    # ph_stat = scipy.stats.boltzmann(lambda_=0.2, N=6)

    # GET SPAD MATRIX AND CLICKS ARRAY
    V, c = np.loadtxt('V.txt'), np.loadtxt('c.txt')

    # SOLVE THE c = Vp_n PROBLEM
    p_n = EMEsolver(V, c, alpha=1e-5, iterations=int(1e10), epsilon=1e-18)

    # FIDELITY
    p_space = np.array(range(D + 1))
    fid = 1 - np.linalg.norm(p_n - ph_stat.pmf(p_space))
    print(f"Fidelity of the reconstruction: {fid:.3f}.")
    print(f"Fidelity of the reconstruction: {prism.fidelity(p_n, ph_stat.pmf(p_space)):.3f}.")

    # PREPARE PLOT
    plt.figure()

    plt.bar(p_space, p_n, alpha=0.6, color="red", label="Analytic + EME")
    plt.plot(ph_stat.pmf(p_space), "o--k", label="Theoretical p(n)")

    plt.legend()
    plt.title(f"Reconstructed p vector")
    plt.xlabel("n")
    plt.ylabel("p(n)")
    plt.ylim([0, 1])
    plt.show()


if __name__ == "__main__":
    main()

