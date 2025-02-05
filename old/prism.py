from __future__ import annotations

import math
import pathlib
import random
import warnings
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
import scipy
from joblib import Parallel, delayed
from tqdm import tqdm

warnings.simplefilter(action="ignore", category=RuntimeWarning)


def fidelity(a, b):
    return np.sum(np.sqrt(np.multiply(a, b)))


def EMEsolver(
    mat: np.ndarray,
    b_vec: np.ndarray,
    alpha: float = 5e-2,
    iterations=int(1e10),
    epsilon=1e-15,
    max_stagnation=int(1e5),
):
    nMax = len(b_vec)
    EME = np.zeros(nMax)

    # Initial guess is a uniform distribution
    pn = np.array([1.0 / float(nMax)] * nMax)
    iteration = 0
    stagnation_counter = 0
    last_dist = float("inf")

    # Calculate the number of orders of magnitude to track
    start_dist = 1e0  # Assume starting `dist` is at least 1 for scaling
    if start_dist <= epsilon:
        raise ValueError("Initial `dist` must be greater than `epsilon`.")
    orders_of_magnitude = int(np.log10(start_dist / epsilon))

    with tqdm(total=orders_of_magnitude, desc="Convergence (log scale)", unit="order") as pbar:
        current_order = orders_of_magnitude

        while iteration < iterations:
            # Expectation-Maximization step
            EM = np.dot(b_vec / np.dot(mat, pn), mat) * pn
            # Additional entropy regularization
            E = alpha * (np.log(pn) - np.sum(pn * np.log(pn))) * pn
            E[np.isnan(E)] = 0.0
            EME = EM - E

            # Compute distance for convergence check
            dist = np.sqrt(np.sum((EME - pn) ** 2))

            # Update progress bar based on reduction in order of magnitude
            new_order = int(np.log10(dist / epsilon)) if dist > epsilon else 0
            if new_order < current_order:
                pbar.update(current_order - new_order)
                current_order = new_order

            # Check convergence criteria
            if dist <= epsilon:
                print("\nConverged to the desired precision.")
                break
            elif dist >= last_dist:
                stagnation_counter += 1
                if stagnation_counter >= max_stagnation:
                    print(f"\nStagnation detected: the distance has not improved for {max_stagnation} iterations.")
                    break
            else:
                stagnation_counter = 0  # Reset stagnation counter

            pn = EME
            last_dist = min(last_dist, dist)
            iteration += 1

    return EME


def binomial(n, m):
    return scipy.special.comb(n, m, exact=True)


def factorial(a):
    return scipy.special.gamma(a + 1)


def loss_matrix(pa: float, num: int) -> np.ndarray:
    mat = np.zeros([num + 1, num + 1])
    for n in range(num + 1):
        for m in range(num + 1):
            mat[n, m] = binomial(m, n) * (pa**n) * ((1 - pa) ** (m - n))
    return mat


def dark_matrix(pd: float, num: int) -> np.ndarray:
    mat = np.zeros((num + 2, num + 2))

    for n in range(num + 2):
        for m in range(n + 1):
            mat[n, m] = pd ** (n - m) * np.exp(-pd) / factorial(n - m)

    # Set the last row such that columns sum to 1, to preserve probability
    mat[-1, :] = 1 - np.sum(mat, axis=0)
    return mat[:-1, :-1]


def xtalk_matrix(xtk: float, num: int) -> np.ndarray:
    if not xtk:
        return np.eye(num + 1)

    mat = np.zeros([num + 1, num + 1])
    for n in range(num + 1):
        for m in range(num + 1):
            mat[n, m] = binomial(m, n - m) * (xtk ** (n - m)) * ((1 - xtk) ** (2 * m - n))
    return mat


@lru_cache(maxsize=None)
def P(D: int, eta: float, dcr: float, N: int, C: int, K: int) -> float:
    if N < 0 or C < 0 or C > N:
        return 0
    elif C == 0 and N == 0:
        return math.comb(D, K) * (dcr**K) * ((1 - dcr) ** (D - K))

    term1 = ((1 - eta) + eta * (K + C) / D) * P(D, eta, dcr, N - 1, C, K)
    term2 = (eta * (D - (K + C - 1)) / D) * P(D, eta, dcr, N - 1, C - 1, K)
    return term1 + term2


def analytic(D: int, eta: float, dcr: float, N: int, L: int):
    return np.sum([P(D, eta, dcr, N, i, L - i) for i in range(L + 1)])


def get_spad_matrix(
    num_det: int,
    eta: float,
    dcr: float,
    xtk: float,
    filename: str | None = None,
    plot: bool = False,
) -> np.ndarray:
    if filename and pathlib.Path(filename).is_file():
        V_th = np.loadtxt(filename)
    else:
        n_values, m_values = np.meshgrid(np.arange(num_det + 1), np.arange(num_det + 1), indexing="ij")
        y_space = np.array(
            [
                analytic(num_det, eta, dcr, m, n)
                for n, m in tqdm(
                    zip(n_values.ravel(), m_values.ravel()),
                    total=(num_det + 1) ** 2,
                    desc="Calculating V_th",
                )
            ]
        )
        V_th = y_space.reshape(num_det + 1, num_det + 1)
        if xtk:
            V_th = xtalk_matrix(xtk, num_det) @ V_th
    if filename:
        np.savetxt(filename, V_th)
    if plot:
        plt.imshow(V_th)
        plt.title("Vth")
        plt.show()
    return V_th


# def simulate_spad_matrix(
#         num_det: int,
#         eta: float,
#         dcr: float,
#         xtk: float,
#         iterations: int = int(1e5),
#         filename: str | None = None,
#         plot: bool = False,
# ):
#     if filename and pathlib.Path(filename).is_file():
#         V = np.loadtxt(filename)
#     else:
#         V = np.zeros((num_det + 1, num_det + 1), dtype=np.float64)
#         for n in tqdm(range(num_det + 1)):
#             y_space = np.zeros(num_det + 1, dtype=int)
#
#             # Precompute all simulations in a single batch
#             simulations = np.array(
#                     [
#                             np.sum(
#                                     simulate(num_det=num_det, num_ph=n, eta=eta, dcr=dcr, xtk=xtk)
#                             )
#                             for _ in range(iterations)
#                     ]
#             )
#
#             counts = np.bincount(simulations, minlength=num_det + 1)
#             y_space[: len(counts)] = counts
#             V[:, n] = y_space / iterations  # Update V matrix
#         if filename:
#             np.savetxt(filename, V)
#     if plot:
#         plt.imshow(V)
#         plt.title("Vsim")
#         plt.show()
#     return V


def get_dcr_array(D: int, dcr_min: float = 1e-3, dcr_max: float = 1e-2) -> tuple[float, np.ndarray]:
    dcr = np.logspace(np.log10(dcr_min), np.log10(dcr_max), D)
    random.shuffle(dcr)
    return np.mean(dcr), dcr


def get_clicks_array(
    num_det: int,
    eta: float,
    dcr: np.ndarray,
    xtk: float,
    ph_stat: scipy.stats.distributions,
    it: int,
    save: bool = False,
) -> np.ndarray:
    # Generate random values
    n_values = np.asarray(ph_stat.rvs(size=it))

    # Define a wrapper function for simulate
    def simulate_sum(n, num_det, eta, dcr, xtk):
        return np.sum(simulate(num_det=num_det, num_ph=n, eta=eta, dcr=dcr, xtk=xtk))

    # Run simulations in parallel
    simulation_results = Parallel(n_jobs=-1)(
        delayed(simulate_sum)(n, num_det, eta, dcr, xtk)
        for n in tqdm(n_values, desc="Running clicks simulations", total=it)
    )
    simulation_results = np.array(simulation_results)
    b = np.zeros(num_det + 1, dtype=float)
    np.add.at(b, simulation_results, 1)  # Count occurrences of each result
    b /= it
    if save:
        np.savetxt("click_vector", b)
    return b


def simulate(num_det: int, num_ph: int, eta: float, dcr: np.ndarray, xtk: float):
    """
    Simulates the click space for a detector array with dark counts, photon counts, and crosstalk.

    Parameters:
        num_det (int): Number of detectors.
        num_ph (int): Number of photons.
        eta (float): Detection efficiency for each detector.
        dcr (np.ndarray): Dark count rate for each detector.
        xtk (float): Crosstalk matrix, where xtr[m, n] is the crosstalk rate from detector n to m.

    Returns:
        np.ndarray: Updated click space after considering dark counts, photon counts, and crosstalk.
    """
    # Initialize click space
    click_space = np.zeros(num_det, dtype=int)

    # Calculate Dark Counts
    if dcr.any():
        random_values = np.random.rand(num_det)
        click_space = np.where((random_values < dcr), 1, click_space)

    # Calculate Photon Counts
    if eta:
        # select random detectors (with repetition) using random indeces
        random_index = np.random.randint(low=0, high=num_det, size=num_ph)
        miss_chances = np.random.rand(num_ph)
        click_space[random_index[miss_chances < eta]] = 1

    # Crosstalk: Victim (m) and Aggressor (n)
    if xtk:
        # Generate random values for each element in xtk
        random_matrix = np.random.rand(num_det)

        # Identify zeros with at least one neighboring 1
        padded_click = np.pad(click_space, (1, 1), "constant", constant_values=(0, 0))
        neighbors = (np.roll(padded_click, 1) == 1) | (np.roll(padded_click, -1) == 1)

        # Find indices where the value is 0 and has a neighbor marked
        mask = (click_space == 0) & neighbors[1:-1]
        click_space[mask] = (random_matrix[mask] < xtk).astype(int)
    return click_space


def main():
    num = 64
    eta = 0.9
    dcr = 0.1

    B = dark_matrix(dcr, num) @ loss_matrix(eta, num)
    plt.matshow(B)
    plt.show()


if __name__ == "__main__":
    main()
