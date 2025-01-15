
#%% IMPORTS

import bisect
import numpy as np
import pickle as pk
import matplotlib.pyplot as plt
from icecream import ic
from tqdm import tqdm
import math

from simulator_functions import *

#%% MAIN FUNCTIONS

def simulate(D, N, eta, dcr, xtr, cum_dist):
    click_space = D*[0]
    
    # Calculate Dark Counts.
    if dcr:
        for i, v in enumerate(click_space):
            if dcr[i] == 0:
                continue
            if np.random.rand() < dcr[i]:
                click_space[i] = 1
    
    # Calculate Photon Counts.
    for n in range(N):
        
        # Draw a detector from the diffusion space.
        index = bisect.bisect_left(cum_dist, np.random.rand()) # Very fast search method!
        detector_id = index - 1
        
        # Calculate chance to miss photon for the specific detector.
        if np.random.rand() < eta[detector_id]:
            click_space[detector_id] = 1
    
    # Crosstalk. Crosstalk rate is calculated as m = victim and n = aggressor.
    if xtr:
        click_space_copy = np.copy(click_space)
        for n, m in np.ndindex(D, D):
            if xtr[n, m] == 0:
                continue
            if np.random.rand() < xtr[n, m] and click_space_copy[n] == 1:
                click_space[m] = 1
    
    return click_space

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
    

#%% MONTE CARLO SIMULATION w/ ANALYTIC

iterations = 100000
plot = True
ideal_eta = False
ideal_dcr = False

D = 64
N = 2

eta_space = np.array(D*[1]) if ideal_eta else np.linspace(0.6, 0.7, D)
dcr_space = np.array(D*[0]) if ideal_dcr else np.logspace(np.log(3*10e-5), np.log(2*10e-4), D)
xtr_space = np.zeros((D, D))

dif_space = normalize(D*[1])
cumulative_dif_space = [np.sum(dif_space[0:i]) for i in range(len(dif_space))]

y_space = (D + 1)*[0]
x_space = list(range(D+1))

for i in tqdm(range(iterations)):
    y_space[np.sum(simulate(D, N, eta_space, dcr_space, xtr_space, cumulative_dif_space))] += 1
y_space_analytic = [analytic(D, np.mean(eta_space), np.mean(dcr_space), N, m) for m in range(D+1)]

print(f"Monte Carlo Simulation: {np.array(y_space) / iterations}")
print(f"Analytic Method: {np.array(y_space_analytic)}")

if plot:
    plt.figure(0, (6, 3))
    
    plt.bar(x_space, [v + 1*(i == N) for i, v in enumerate((D + 1)*[0])],
            alpha = 0.6, color = 'C1', edgecolor = 'C1', linewidth = 1,
            label = 'Fock State')
    plt.bar(x_space, np.array(y_space) / iterations,
            alpha = 0.6, color = 'C0', edgecolor = 'C0', linewidth = 1,
            label = 'Monte Carlo')
    
    plt.plot(x_space, y_space_analytic, 'o', color = 'black', markersize = 3,
            label = 'Analytic')
    
    plt.xlabel("n")
    plt.ylabel("Probability")
    
    plt.xlim(-0.9, 20.9)
    plt.xticks(list(range(0, 21)))

    plt.title(f"Click distribution for {N} photons with {D} detectors")
    
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    #specify order of items in legend
    order = [1,2,0]

    #add legend to plot
    # plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order]) 
    plt.legend()
    
    filename = rf'E:\python\work\SPAD\article_graphs\distribution_{N}p_{D}d'
    filename += '_eta' if not ideal_eta else ''
    filename += '_dcr' if not ideal_dcr else ''
    filename += '.png'
    plt.savefig(filename, bbox_inches = 'tight')
    plt.show()

#%% MONTE CARLO SIMULATION w/o ANALYTIC

iterations = 10000
plot = True
ideal_eta = False
ideal_dcr = False
ideal_xtr = False

D = 64
N = 0

eta_space = np.array(D*[1]) if ideal_eta else np.linspace(0.6, 0.7, D)
dcr_space = np.array(D*[0]) if ideal_dcr else np.logspace(np.log(3*10e-5), np.log(2*10e-4), D)
xtr_space = np.zeros((D, D)) if ideal_xtr else rectangular_array_xtr(8, 8, 0.01)

dif_space = normalize(D*[1])
cumulative_dif_space = [np.sum(dif_space[0:i]) for i in range(len(dif_space))]

y_space = (D + 1)*[0]
x_space = list(range(D+1))

for i in tqdm(range(iterations)):
    y_space[np.sum(simulate(D, N, eta_space, dcr_space, xtr_space, cumulative_dif_space))] += 1

print(f"Monte Carlo Simulation: {np.array(y_space) / iterations}")

if plot:
    plt.figure(0, (6, 3))
    
    plt.bar(x_space, [v + 1*(i == N) for i, v in enumerate((D + 1)*[0])],
            alpha = 0.6, color = 'C1', edgecolor = 'C1', linewidth = 1,
            label = 'Input State')
    plt.bar(x_space, np.array(y_space) / iterations,
            alpha = 0.6, color = 'C0', edgecolor = 'C0', linewidth = 1,
            label = 'Monte Carlo')
    
    plt.xlabel("n")
    plt.ylabel("Probability")
    
    # plt.xlim(-0.9, 20.9)
    # plt.xticks(list(range(0, 21)))

    plt.title(f"Click distribution for {N} photons with {D} detectors")
    
    #get handles and labels
    handles, labels = plt.gca().get_legend_handles_labels()

    #specify order of items in legend
    order = [1,0]

    #add legend to plot
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order]) 
    
    filename = rf'E:\python\work\SPAD\article_graphs\distribution_{N}p_{D}d'
    filename += '_eta' if not ideal_eta else ''
    filename += '_dcr' if not ideal_dcr else ''
    filename += '_xtr' if not ideal_dcr else ''
    filename += '.png'
    plt.savefig(filename, bbox_inches = 'tight')
    plt.show()

# %% RECONSTRUCT

V_iterations = 100000
b_iterations = 1000000
plot = True
ideal_eta = False
ideal_dcr = False

D = 32
eta_space = np.array(D*[1]) if ideal_eta else np.linspace(0.65, 0.7, D)
dcr_space = np.array(D*[0]) if ideal_dcr else np.logspace(np.log(3*10e-3), np.log(2*10e-2), D) 
dif_space = normalize(D*[1])
cumulative_dif_space = [np.sum(dif_space[0:i]) for i in range(len(dif_space))]

V = np.zeros((D+1, D+1))
for n in tqdm(range(D + 1)):
    y_space = (D + 1)*[0]
    for t in (range(V_iterations)):
        y_space[np.sum(simulate(D, n, eta_space, dcr_space, cumulative_dif_space))] += 1
        V[:, n] = (np.array(y_space) / V_iterations).T
plt.imshow(V)
plt.show()

V_th = np.zeros((D+1, D+1))
for n in tqdm(range(D + 1)):
    y_space = [analytic(D, np.mean(eta_space), np.mean(dcr_space), n, m) for m in range(D+1)]
    V_th[:, n] = (np.array(y_space)).T
plt.imshow(V_th)
plt.show()
    
# Simulate mixed state
n_1 = 4
n_2 = 8

y_space = (D + 1)*[0]
x_space = list(range(D + 1))
for t in tqdm(range(b_iterations)):
    n = n_1 if (np.random.rand() < 0.5) else n_2
    y_space[np.sum(simulate(D, n, eta_space, dcr_space, cumulative_dif_space))] += 1
b = np.array(y_space) / b_iterations

x = ic(np.linalg.lstsq(V, b))
plt.bar(x_space, x[0], alpha = 0.6, color = 'C0', edgecolor = 'C0', linewidth = 1)

plt.title(f"Reconstructed p vector, Monte Carlo")
plt.xlabel("n")
plt.ylabel("Weight")
plt.ylim(-0.2, 0.6)

plt.show()


plt.title(f"Reconstructed p vector, Analytic")
plt.xlabel("n")
plt.ylabel("Weight")
plt.ylim(-0.2, 0.6)

x_th = ic(np.linalg.solve(V_th, b))
plt.bar(x_space, x_th, alpha = 0.6, color = 'C0', edgecolor = 'C0', linewidth = 1)
plt.show()

# %%

def monte_carlo(it, D, N, pde, dcr, xtr, array_type = 'square'):
    eta_space = pde * np.array(D*[1])
    dcr_space = dcr * np.array(D*[1])
    xtr_space = rectangular_array_xtr(int(np.sqrt(D)), int(np.sqrt(D)), xtr)

    dif_space = normalize(D*[1])
    cumulative_dif_space = [np.sum(dif_space[0:i]) for i in range(len(dif_space))]

    y_space = (D + 1)*[0]

    for i in tqdm(range(it)):
        y_space[np.sum(simulate(D, N, eta_space, dcr_space, xtr_space, cumulative_dif_space))] += 1
    
    return np.array(y_space) / it

# %%
