import copy
import random

import numpy as np


def find_nbors(arr):
    padded_click = np.pad(arr, (1, 1), "constant", constant_values=(0, 0))
    nbors = (np.roll(padded_click, 1) == 1).astype(int) + (np.roll(padded_click, -1) == 1).astype(int)
    return nbors


def print_events(post, mask):

    print("neighbors\t[", end=" ")
    for e in mask:
        if e == 1:
            print('v', end=" ")
        if e == 2:
            print('x', end=" ")
        if e == 0:
            print(' ', end=" ")
    print("]")

    print("config + xtalk \t[", end=" ")
    for e in post:
        if e:
            print("*", end=" ")
        else:
            print("_", end=" ")
    print("]")


def main():
    D = 40
    xtk = 0.5
    # random.seed(10)

    a = np.array([random.random() < 0.2 for _ in range(D)]).astype(int)
    mask = copy.copy(a)
    a_old = np.zeros_like(a)

    # print("click config \t[", end=" ")
    # for e in a:
    #     print(e, end=" ")
    # print("]")

    print("init \t\t[", end=" ")
    for e in a:
        if e:
            print("*", end=" ")
        else:
            print("_", end=" ")
    print("]")

    while True:
        nbors = find_nbors(a - a_old)
        mask = np.where(mask == 0, nbors[1:-1], 0)

        if not any(mask):
            break

        a_old = copy.copy(a)
        for i, num_events in enumerate(mask):
            for _ in range(num_events):
                if random.random() < xtk:
                    a[i] = 1
        print_events(a, mask)
        mask = a_old + a

        # print("mask\t\t[", end=" ")
        # for e in mask:
        #     if e:
        #         print("#", end=" ")
        #     else:
        #         print(" ", end=" ")
        # print("]")


if __name__ == "__main__":
    main()
