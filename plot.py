import matplotlib.pyplot as plt
import numpy as np


def mean(arr):
    return np.sum([i * val for i, val in enumerate(arr)])


def variance(arr):
    return np.sum([val * (i - mean(arr)) ** 2 for i, val in enumerate(arr)])


def main():
    c0 = np.loadtxt('c0.txt')
    cx = np.loadtxt('cx.txt')
    p_space = np.array(range(len(c0)))

    print(mean(c0), variance(c0))
    print(mean(cx), variance(cx))

    eps = 0.1
    n = 2
    ENF = 1 + eps + (3 * n - 3) / (2 * n) * eps**2
    borel = 1 / (1 + np.log(1 - eps))
    print(f'{ENF=}')
    print(f'{borel=}')
    print(mean(cx) / mean(c0))
    p

    # # PLOT
    # plt.figure()
    #
    # plt.plot(p_space, c0, "ro-", label="no x-talk")
    # plt.plot(p_space, cx, "bo-", label="x-talk")
    #
    # # plt.bar(p_space, c0, alpha=0.6, color="red", label="no x-talk")
    # # plt.bar(p_space, cx, alpha=0.6, color="blue", label="x-talk")
    #
    # plt.legend()
    # plt.xlabel("n")
    # plt.ylabel("freq")
    # # plt.ylim([0, 1])
    # plt.show()


if __name__ == "__main__":
    main()
