import itertools
import os
import random

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

path = "../data/random"
if not os.path.exists(path):
    os.makedirs(path)

seeds = [0, 1, 100, 999, 1234]
k = 3
for seed in seeds:
    np.random.seed(seed)
    ns = [1000, 2500, 5000, 10000, 25000, 50000, 75000, 100000, 125000, 150000]
    ps = [
        0.01,
        0.01,
        0.001,
        0.001,
        0.0001,
        0.0001,
        0.00005,
        0.00005,
        0.000001,
        0.000001,
    ]
    for n, p in zip(ns, ps):
        print(n, p)
        G = nx.fast_gnp_random_graph(n, p, seed=seed)
        iso = list(nx.isolates(G))
        G.remove_nodes_from(iso)
        nx.write_edgelist(G, f"{path}/er-n{n:06}_p{p}_s{seed:04}.txt", data=False)
        G = nx.barabasi_albert_graph(n, k)
        nx.write_edgelist(G, f"{path}/ba-n{n:06}_k{k}_s{seed:04}.txt", data=False)
