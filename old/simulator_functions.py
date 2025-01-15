
#%% IMPORTS

import bisect
import numpy as np
import scipy as sp
import pickle as pk
import matplotlib.pyplot as plt
from icecream import ic
from tqdm import tqdm
import math
from scipy import signal

#%% MAIN FUNCTIONS

def simulate(D, N, eta, dcr, xtr, dif, xto = -1):
    click_space = D*[0]
    cumulative_dif_space = [np.sum(dif[0:i]) for i in range(len(dif))] if type(dif)!=type(None) else None
    
    # Calculate Dark Counts.
    if np.any(dcr):
        for i, v in enumerate(click_space):
            if dcr[i] == 0:
                continue
            if np.random.rand() < dcr[i]:
                click_space[i] = 1
    
    # Calculate Photon Counts.
    for n in range(N):
        
        # Draw a detector from the diffusion space.
        if cumulative_dif_space:
            index = bisect.bisect_left(cumulative_dif_space, np.random.rand()) # Very fast search method!
            detector_id = index - 1
        else:
            detector_id = np.random.randint(0, D)
        
        # Calculate chance to miss photon for the specific detector.
        if np.random.rand() < eta[detector_id]:
            click_space[detector_id] = 1
    
    # Crosstalk. Crosstalk rate is calculated as m = victim and n = aggressor.
    if np.any(xtr):
        if xto <= 0:
            click_space_copy = np.copy(click_space)
            for n, m in np.ndindex(D, D):
                if click_space_copy[n] == 0:
                    continue
                if xtr[n, m] == 0:
                    continue
                if np.random.rand() < xtr[n, m]:
                    click_space[m] = 1
        else:
            for o in range(xto):
                click_space_copy = np.copy(click_space)
                A = sp.sparse.coo_array(xtr)
                for n, m, v in zip(A.row, A.col, A.data):
                    if click_space_copy[n] == 0:
                        continue
                    if np.random.rand() < v:
                        click_space[m] = 1
                A = A @ A
    
    return click_space

def normalize(vector):
    v = np.array(vector)
    lumpsum = np.sum(v)
    if lumpsum == 0: return 0
    else: return v / lumpsum

def P(D:int, eta:float, dcr:float, N, C, K, memo = None):
    # Probability of C actual clicks, K dark clicks for N photons.
    
    # Necessary, because Python is one of the programming languages of all time
    if memo is None:
        memo = {}
    
    if (N, C, K) in memo:
        return memo[(N, C, K)]
    
    if N < 0:
        return 0
    if C < 0:
        return 0
    if C > N:
        return 0
    elif C == 0 and N == 0:
        return math.comb(D, K) * (dcr**K) * ((1 - dcr)**(D-K))
    
    result = ((1 - eta) + eta * (K + C) / D) * P(D, eta, dcr, N - 1, C, K, memo) + \
           + (eta * (D - (K + C - 1)) / D ) * P(D, eta, dcr, N - 1, C - 1, K, memo)
    
    memo[(N, C, K)] = result
    return result

def analytic(D:int, eta:float, dcr:float, N, L):
    return np.sum([P(D, eta, dcr, N, i, L - i) for i in range(L + 1)])

def linear_array_xtr(D:int, xtr:float):
    to_return = np.zeros((D, D))
    np.fill_diagonal(to_return[:-1, 1:], xtr)
    np.fill_diagonal(to_return[1:, :-1], xtr)
    return to_return

def rectangular_array_xtr(rows:int, cols:int, xtr:float):
    D = rows*cols
    to_return = np.zeros((D, D))
    np.fill_diagonal(to_return[:-1, 1:], xtr)
    np.fill_diagonal(to_return[1:, :-1], xtr)
    np.fill_diagonal(to_return[:-cols, cols:], xtr)
    np.fill_diagonal(to_return[cols:, :-cols], xtr)
    for k in range(1, rows):
        to_return[k*cols, k*cols-1] = 0
        to_return[k*cols-1, k*cols] = 0
    return to_return

def fidelity(a, b):
    return np.sum(np.sqrt(np.multiply(a, b)))
 
def monte_carlo(it, D, N, pde, dcr, xtr, array_type = 'square', xt_order = -1):
    eta_space = pde * np.array(D*[1]) if pde is float else pde
    dcr_space = dcr * np.array(D*[1]) if dcr is float else dcr
    xtr_space = rectangular_array_xtr(int(np.sqrt(D)), int(np.sqrt(D)), xtr) if xtr is float else xtr

    dif_space = None

    y_space = (D + 1)*[0]

    for i in tqdm(range(it)):
    # for i in range(it):
        y_space[np.sum(simulate(D, N, eta_space, dcr_space, xtr_space, dif_space, xto=xt_order))] += 1
    
    return np.array(y_space) / it

def monte_carlo_fidelity_condition(fid_limit, batch, D, N, pde, dcr, xtr, array_type = 'square', xt_order = -1):
    eta_space = pde * np.array(D*[1]) if pde is float else pde
    dcr_space = dcr * np.array(D*[1]) if dcr is float else dcr
    xtr_space = rectangular_array_xtr(int(np.sqrt(D)), int(np.sqrt(D)), xtr) if xtr is float else xtr

    dif_space = None

    y_space = (D + 1)*[0]
    new_y_space = (D + 1)*[0]
    
    i = 0
    fid = 0
    while fid < fid_limit:
        for j in range(batch):
            new_y_space[np.sum(simulate(D, N, eta_space, dcr_space, xtr_space, dif_space, xto=xt_order))] += 1
        fid = fidelity(normalize(y_space), normalize(new_y_space))
        print(f"Batch #{i} \t Fidelity: {fid}")
        i += 1
        y_space = new_y_space[:] # This makes a copy! Without [:] we do pointer assignment
    
    return np.array(normalize(y_space))


# %%

