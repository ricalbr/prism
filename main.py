from __future__ import annotations


import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as sst
import prism as pp

def main():
    # SPAD ARRAY/MATRIX PROPERTIES
    D = 64
    eta = 0.80
    mean_dcr, dcr = pp.get_dcr_array(num_det=D, dcr_min=1e-3, dcr_max=1e-2)
    xtk = 0.05
    iterations = int(1e5)

    # INPUT PHOTON STATISTICS
    ph_stat = sst.rv_discrete(values=([20], [1]))
    # ph_stat = sst.rv_discrete(values=([2], [1]))
    # ph_stat = sst.rv_discrete(values=([1, 3], [0.5, 0.5]))
    # ph_stat = sst.poisson(mu=40)
    # ph_stat = sst.poisson(mu=2)
    # ph_stat = sst.boltzmann(lambda_=0.2, num_ph=6)
    # ph_stat = sst.rv_discrete(values=([1, 3, 5, 7, 9, 11, 13], np.ones((7,))/7))

    # GET SPAD MATRIX AND CLICKS ARRAY
    V = pp.get_spad_matrix(num_det=D, eta=eta, dcr=mean_dcr, xtk=xtk)
    c = pp.get_clicks_array(num_det=D, eta=eta, dcr=dcr, xtk=xtk, ph_stat=ph_stat, it=iterations)
    # np.savetxt('c_vec', c)
    # c = np.loadtxt('c_vec')

    # SOLVE THE c = Vp_n PROBLEM
    p_n = pp.EMEsolver(V, c, alpha=1e0, epsilon=1e-18, max_stagnation=1e6)

    # FIDELITY
    p_space = np.array(range(D + 1))
    print(f"Fidelity of the reconstruction: {pp.fidelity(p_n, ph_stat.pmf(p_space)):.3f}.")

    # PREPARE PLOT
    plt.figure()

    plt.bar(p_space, p_n, alpha=0.6, color="red", label="Analytic + EME")
    plt.plot(ph_stat.pmf(p_space), "o--k", label="Theoretical p(n)")

    plt.legend()
    plt.title(f"Reconstructed p vector")
    plt.xlabel("n")
    plt.ylabel("p(n)")
    plt.show()


if __name__ == "__main__":
    main()
