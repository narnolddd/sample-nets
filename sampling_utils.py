import networkx as nx
import random

from scipy import special
from numpy import exp
from numpy.lib.scimath import log
lgam = special.gammaln

# Binomial function that bypasses issues with numerics for larger n.
def binom(n:int, k:int, p):
    if p!=1.0:
        return exp(lgam(n+1) - lgam(n-k+1) - lgam(k+1) + k*log(p) + (n-k)*log(1.-p))
    else:
        if k==n:
            return 1.0
        else:
            return 0.0
        
# Returns an edge-sampled graph with retention probability p.
def edge_sample(G,p,random_seed=None):
    # optionally seed for reproducibility
    random.seed(random_seed)
    H = nx.Graph()
    for (u,v) in G.edges:
        if random.random()<p:
            H.add_edge(u,v)
    return H