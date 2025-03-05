from functools import lru_cache
from math import comb

import numpy as np
import scipy

np.set_printoptions(precision=3, suppress=True)


@lru_cache(None)
def h(n, c):
    if c == 0:
        return 1
    result = 0
    for s in range(1, c + 1):
        coeff = comb(n, s)  # Calcola il coefficiente binomiale (n su s)
        inner_sum = recursive_inner_sum(n, c, s, ())  # Usa tuple invece di liste
        result += coeff * inner_sum

    return result


@lru_cache(None)
def recursive_inner_sum(n, c, s, indices):
    # Se abbiamo raggiunto il livello più profondo della somma ricorsiva
    if len(indices) == s - 1:
        remaining = c - s - sum(indices)
        last_h = h(n, indices[-1]) if indices else 1
        return last_h * h(n, remaining) if remaining >= 0 else 0

    # Itera su tutti i valori possibili per il prossimo indice
    start = 0 if not indices else indices[-1]
    end = c - s - sum(indices) + 1
    result = 0
    for i in range(start, end):
        result += recursive_inner_sum(n, c, s, indices + (i,))  # Usa tuple invece di liste

    return result


def binomial(n, m):
    return scipy.special.comb(n, m, exact=True)


def xtalk_matrix(xtk: float, num: int) -> np.ndarray:
    if not xtk:
        return np.eye(num + 1)

    mat = np.zeros([num + 1, num + 1])
    for n in range(num + 1):
        for m in range(num + 1):
            mat[n, m] = binomial(m, n - m) * (xtk ** (n - m)) * ((1 - xtk) ** (2 * m - n))
    return mat


@lru_cache
def p_single(p: float, k: int, n: int = 2) -> float:
    if k <= 5:
        return h(n, k - 1) * p ** (k - 1) * (1 - p) ** (n * k - k + 1)
    else:
        p5 = p_single(p=p, k=5, n=n)
        s = np.sum([p_single(p=p, k=i, n=n) for i in range(1, 5)])
        return p5 * (1 - p5 / (1 - s)) ** (k - 5)


@lru_cache
def p_m(p: float, k: int, m: int, n: int = 2) -> float:
    if m < 1:
        return p_single(p, k, n)
    else:
        return np.sum([p_m(p=p, k=k - i, m=m - 1, n=n) * p_single(p, i, n) for i in range(1, k - m + 1)])


def x_gal(eps, num):
    mat = np.eye(num + 1)
    xt_prob = 1 - np.sqrt(1 - eps)

    for i in range(num + 1):
        for j in range(num):
            mat[i, j + 1] = p_m(xt_prob, i, j)
    return mat


def main():

    # n = 2
    K = 64
    eps = 0.05

    # p = 1 - np.sqrt(1 - eps)
    # for k in range(1, K):
    #     print(f"P1({k})= {p_single(p, k):.4f}")

    # print(p_m(p, 1, 5))
    # print(P(1, 5, p))

    x_art = x_gal(eps, K)
    np.savetxt("x.txt", x_art)

    # confronto articolo
    # print(f"P(1)= {(1-p)**n}")
    # print(f"P(2)= {n*p*(1-p)**(2*n-1)}")
    # print(f"P(3)= {0.5*n*(3*n-1)*p**2*(1-p)**(3*n-2)}")
    # print(f"P(4)= {1/3*(8*n**2-6*n+1)*p**3*(1-p)**(4*n-3)}")
    # print(f"P(5)= {1/4*n*(125/6*n**3 - 25*n**2+55/6*n-1)*p**4*(1-p)**(5*n-4)}")

    # XT = xtalk_matrix(eps, K)
    # print(XT)


if __name__ == "__main__":
    main()
