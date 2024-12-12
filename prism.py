from __future__ import annotations

import math
import pathlib
import random
import warnings
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
import scipy
from joblib import delayed, Parallel
from tqdm import tqdm

warnings.simplefilter(action="ignore", category=RuntimeWarning)


def fidelity(a, b):
    return np.sum(np.sqrt(np.multiply(a, b)))


def EMEsolver(
    mat: np.ndarray,
    b_vec: np.ndarray,
    alpha: float = 5e-2,
    iterations: int = int(1e10),
    epsilon: float = 1e-15,
    max_stagnation: int = int(1e5),
) -> np.ndarray:
    """Estimation-Maximization algorithm with Entropy.

    Implementation of the EME algorithm from Hlousek et al. with minor tweaks for convergence stagnations.

    Parameters
    ----------
    mat: np.ndarray
        Transformation matrix.
    b_vec: np.ndarray
        Right-hand-side vector.
    alpha: float
        Regularization term for the EME algorithm. Default 0.05.
    iterations: int
        Number of iterations. Default 10^10.
    epsilon: float
        Precision for testing the convergence of the algorith. Default 10^-15.
    max_stagnation: int
        Maxmimum number of allow iterations without improvement of the convergence. Default 10^5.

    Returns
    -------
    np.ndarray
        Estimation of the unknown variable of the linear system.
    """

    nMax = len(b_vec)
    EME = np.zeros(nMax)

    # Initial guess is a uniform distribution
    pn = np.array([1.0 / float(nMax)] * nMax)
    iteration = 0
    stagnation_counter = 0
    last_dist = float("inf")

    # Calculate the number of orders of magnitude to track
    start_dist = 1.0  # Assume starting `dist` is at least 1 for scaling
    if start_dist <= epsilon:
        raise ValueError("Initial `dist` must be greater than `epsilon`.")
    orders_of_magnitude = int(np.log10(start_dist / epsilon))

    with tqdm(
        total=orders_of_magnitude, desc="Convergence (log scale)", unit="order"
    ) as pbar:
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
                    print(
                        f"\nStagnation detected: the distance has not improved for {max_stagnation} iterations."
                    )
                    break
            else:
                stagnation_counter = 0  # Reset stagnation counter

            pn = EME
            last_dist = min(last_dist, dist)
            iteration += 1

    return EME


def binomial_coefficient(n: int, m: int) -> int:
    """Binomial coefficient.

    Parameters
    ----------
    n: int
        Top number of the binomial coefficient.
    m: int
        Bottom number of the binomial coefficient.

    Returns
    -------
    int
        Binomial coefficient "n over m".

    """
    return scipy.special.comb(n, m, exact=True)


def factorial(a: float) -> float:
    """Factorial function.

    The factorial is computed using the gamma function.
    If the input number is an integer, the output coincides with the factorial of that number.

    Parameters
    ----------
    a: float
        Input number.

    Returns
    -------
    float
        Gamma(a + 1)
    """
    return scipy.special.gamma(a + 1)


def loss_matrix(eta: float, num: int) -> np.ndarray:
    """Loss matrix.

    Compute the loss matrix for a set of _num_ photons and a finite SPAD efficiency.

    Parameters
    ----------
    eta: float
        SPAD efficiency.
    num: int
        Number of photons.

    Returns
    -------
    np.ndarray
        Loss matrix.
    """
    mat = np.zeros([num + 1, num + 1])
    for n in range(num + 1):
        for m in range(num + 1):
            mat[n, m] = binomial_coefficient(m, n) * (eta**n) * ((1 - eta) ** (m - n))
    return mat


def dark_matrix(dcr: float, num: int) -> np.ndarray:
    """Dark-count matrix.

    Compute the dark-count matrix for a set of _num_ photons and a dark-count rate of k.

    Parameters
    ----------
    dcr: float
        Mean dark count rate for the SPAD array/matrix.
    num: int
        Number of photons.

    Returns
    -------
    np.ndarray
        Dark-count matrix.
    """
    mat = np.zeros((num + 2, num + 2))

    for n in range(num + 2):
        for m in range(n + 1):
            mat[n, m] = dcr ** (n - m) * np.exp(-dcr) / factorial(n - m)

    # Set the last row such that columns sum to 1, to preserve probability
    mat[-1, :] = 1 - np.sum(mat, axis=0)
    return mat[:-1, :-1]


def xtalk_matrix(xtk: float, num: int) -> np.ndarray:
    """Cross-talk matrix.

    Compute the cross-talk matrix for a set of _num_ photons and a cross-talk probability of xtlk.

    Parameters
    ----------
    xtk: float
        Cross-talk probability of the SPAD array/matrix.
    num: int
        Number of photons.

    Returns
    -------
    np.ndarray
        Cross-talk matrix.
    """
    if not xtk:
        return np.eye(num + 1)

    mat = np.zeros([num + 1, num + 1])
    for n in range(num + 1):
        for m in range(num + 1):
            mat[n, m] = (
                binomial_coefficient(m, n - m)
                * (xtk ** (n - m))
                * ((1 - xtk) ** (2 * m - n))
            )
    return mat


@lru_cache(maxsize=None)
def P(
    num_det: int, eta: float, dcr: float, num_ph: int, clicks: int, num_dark_counts: int
) -> float:
    """Photon-count probability function.

    Photon-count probability function for accounting a number C of clicks given num_ph photons, num_dark_counts dark-counts events (
    avalanches) and a finite number of detectors num_det (with efficiency eta, and dark-count rate dcr).

    Parameters
    ----------
    num_det: int
        Number of detectors.
    eta: float
        SPAD efficiency.
    dcr: float
        Mean dark count rate for the SPAD array/matrix.
    num_ph: int
        Number of photons.
    clicks: int
        Number of clicks.
    num_dark_counts: int
        Number of dark clicks events.

    Returns
    -------
    float
        The probability of having C clicks given the tuple of physical parameters (num_ph, num_dark_counts, num_det, eta, dcr).
    """

    if num_ph < 0 or clicks < 0 or clicks > num_ph:
        return 0
    elif clicks == 0 and num_ph == 0:
        return (
            math.comb(num_det, num_dark_counts)
            * (dcr**num_dark_counts)
            * ((1 - dcr) ** (num_det - num_dark_counts))
        )

    term1 = ((1 - eta) + eta * (num_dark_counts + clicks) / num_det) * P(
        num_det, eta, dcr, num_ph - 1, clicks, num_dark_counts
    )
    term2 = (eta * (num_det - (num_dark_counts + clicks - 1)) / num_det) * P(
        num_det, eta, dcr, num_ph - 1, clicks - 1, num_dark_counts
    )
    return term1 + term2


def analytic(num_det: int, eta: float, dcr: float, num_ph: int, num_dark_counts: int):
    """Analitic function to compute the matrix elements of the detection matrix.

    Parameters
    ----------
    num_det: int
        Number of detectors.
    eta: float
        SPAD efficiency.
    dcr: float
        Mean dark count rate for the SPAD array/matrix.
    num_ph: int
        Number of photons.
    num_dark_counts: int
        Number of dark clicks events.

    Returns
    -------
    float
        Matrix element (num_ph,num_dark_counts) of the detection matrix.
    """
    return np.sum(
        [
            P(
                num_det=num_det,
                eta=eta,
                dcr=dcr,
                num_ph=num_ph,
                clicks=i,
                num_dark_counts=num_dark_counts - i,
            )
            for i in range(num_dark_counts + 1)
        ]
    )


def get_spad_matrix(
    num_det: int,
    eta: float,
    dcr: float,
    xtk: float,
    filename: str | None = None,
    plot: bool = False,
) -> np.ndarray:
    """SPAD detection matrix.

    If the filename is given (and the file is present in the given direcotry), the function tries to load the matrix
    from the disk.
    Otherwise, the detection matrix is computed using the anlytical model. The analytical treatment does not model
    cross-talk, which is computed separaterly and multiplied to the overall detection matrix using matrix-to-matrix
    product.
    The user have the possibility to plot and save (export) the computed matrix to a given filename.

    Parameters
    ----------
    num_det: int
        Number of detectors.
    eta: float
        SPAD efficiency.
    dcr: float
        Mean dark count rate for the SPAD array/matrix.
    xtk: float
        Cross-talk probability.
    filename: str | None
        Filename to save the matrix in an external file. Optional, default is None.
    plot: bool
        Flag to plot the matrix. Default is False.

    Returns
    -------
    np.ndarray
        SPAD detection matrix.
    """

    if filename and pathlib.Path(filename).is_file():
        V_th = np.loadtxt(filename)
    else:
        n_values, m_values = np.meshgrid(
            np.arange(num_det + 1), np.arange(num_det + 1), indexing="ij"
        )
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


def get_dcr_array(
    num_det: int, dcr_min: float = 1e-3, dcr_max: float = 1e-2
) -> tuple[float, np.ndarray]:
    """Get dark-count rates array.

    Given a series of SPADs the dark count rate varies among the ensamble of detectors following (typically) a
    log-distribution. The function generate values of dark-count rates, one for each detector from two user-defined
    values.

    Parameters
    ----------
    num_det: int
        Number of detectors.
    dcr_min: float
        Minimum value for the dark-count rate. Default 10^-3.
    dcr_max
        Maximum value for the dark-count rate. Default 10^-2.

    Returns
    -------
    tuple[float, np.ndarray]
        Mean value of the dark-count rate, Dark-count rate array shuffled.
    """
    dcr = np.logspace(np.log10(dcr_min), np.log10(dcr_max), num_det)
    random.shuffle(dcr)
    return np.mean(dcr), dcr


def get_clicks_array(
    num_det: int,
    eta: float,
    dcr: np.ndarray,
    xtk: float,
    ph_stat: scipy.stats.rv_discrete | scipy.stats.rv_continuous,
    it: int,
) -> np.ndarray:
    """Simulates clicks.

    The function simulates clicks events for a given photon statistics and a given number of detectors with efficiency
    eta and dark-count rates (dcr).

    First, a list of random values are drawn from the photon distributions representing the number of photons preesnt
    in each pulse impinging in the system.
    For each of this values a simulation is run considering a finite number of detectors, efficiency, dark count
    rate, and cross-talk probability.

    Parameters
    ----------
    num_det: int
        Number of detectors.
    eta: float
        SPAD efficiency.
    dcr: float
        Mean dark count rate for the SPAD array/matrix.
    xtk: float
        Probability of cross-talk (between first nearest neighbours).
    ph_stat: scipy.stats.rv_discrete | scipy.stats.rv_continuous
        Photon statistics.
    it: int
        Number of iterations of the simulation.

    Returns
    -------
    np.ndarray
        Array of the clicks relative frequencies.
    """

    def simulate_sum(n, num_det, eta, dcr, xtk):
        return np.sum(simulate(num_det=num_det, num_ph=n, eta=eta, dcr=dcr, xtk=xtk))

    # Generate random values from the photon_distribution and run simulations in parallel
    n_values = np.asarray(ph_stat.rvs(size=it))
    simulation_results = Parallel(n_jobs=-1)(
        delayed(simulate_sum)(n, num_det, eta, dcr, xtk)
        for n in tqdm(n_values, desc="Running clicks simulations", total=it)
    )
    simulation_results = np.array(simulation_results)
    b = np.zeros(num_det + 1, dtype=float)
    np.add.at(b, simulation_results, 1)  # Count occurrences of each result
    return b / it


def simulate(num_det: int, num_ph: int, eta: float, dcr: np.ndarray, xtk: float):
    """Simulates the click space for a detector array with dark counts, photon counts, and crosstalk.

    Parameters
    ----------
    num_det: int
        Number of detectors.
    num_ph: int
        Number of photons.
    eta: float
        Detection efficiency for each detector.
    dcr: np.ndarray
        Dark count rates for each detector.
    xtk: float
        Crosstalk probability.

    Returns
    -------
    np.ndarray
        Updated click space after considering dark counts, photon counts, and crosstalk.
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
