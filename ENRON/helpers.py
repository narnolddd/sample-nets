#Simple function which returns from a given networkx graph G a graph H, the induced subgraph from a random proportion p of H's edges.
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random

def edge_sample(G, p):
    edges = G.edges()
    sampled=[]
    for e in edges:
        r = random.random()
        if r < p:
            sampled.append(e)
    H = nx.Graph()
    H.add_edges_from(sampled)
    return H

# Returns degree distribution as a cdf
def deg_cdf(G):
    cts = Counter(sorted([d for n, d in G.degree()]))
    degs = sorted(cts.keys())
    vals = [cts[val] for val in degs]
    cdf = np.cumsum(vals)
    ccdf = cdf[-1] - cdf
    return degs, cdf

# Returns an empirical cdf for discrete integer data
def ecdf(a):
    x, counts = np.unique(a, return_counts=True)
    cusum = np.cumsum(counts)
    return x, cusum

# Sets up a single figure and axes plot
def setup_axes(x,y,xscale="linear",yscale="linear"):
    fig, ax = plt.subplots(figsize = a4_dims)
    ax.set_xlabel(x,fontsize=20)
    ax.set_ylabel(y, fontsize=20)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    return fig, ax

def gen_distribution(k,p,kmax):
    probs = np.zeros(kmax+1)
    for deg in range(k,kmax + 1):
        probs[deg] = comb(deg,k) * np.power(p,k) * np.power(1-p,deg - k)
    return probs/np.sum(probs)

def sample(dist):
    r = random.random()
    tot = 0.0
    ind = 0
    while tot < r:
        tot += dist[ind]
        ind +=1
    return ind