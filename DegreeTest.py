import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
import math

from collections import Counter
from scipy.stats import poisson, zipf
from math import comb

a4_dims = (11.69,8.27)

def count2freq(arr):
    counter = Counter(arr)
    tot = sum(counter.values())
    for k in counter:
        counter[k]/=tot
    return counter

def edge_sample(G, p, seed=None):
    edges = G.edges()
    sampled=[]
    if seed is not None:
        random.Seed(seed)
    for e in edges:
        r = random.random()
        if r < p:
            sampled.append(e)
    H = nx.Graph()
    H.add_edges_from(sampled)
    return H

def degree_dict(G):
    degs={}
    for u in G.nodes():
        degs[u]=G.degree(u)
    return degs

def binom_prob(k,kp,p):
    if k < kp:
        return 0.0
    return comb(k,kp) * np.power(p,kp) *np.power(1-p,k-kp)

def deg_likelihood(k,kp,p):
    if k<kp:
        return 0.0
    den = 1 - np.power(1-p,k)
    num = comb(k,kp) * np.power(p,kp) *np.power(1-p,k-kp)
    return num/den

## This is the function which gets the estimated degree sequence with Bayes using the prior and LOTP to get the denominator
def get_deg_MSE_LOTP(G_deg,H,p):
    
    # H nodes and edges
    m = len(H.edges())
    n = len(H.nodes())
    
    # G nodes and edges
    M = sum(G_deg.values())/2
    N = len(G_deg)
    
    # Estimated and actual removed edges
    edges_estimated = round(m/p)
    missing_edges = edges_estimated - m
    removed_edges = M - m
    removed_nodes = N - n
    
    # Probability of being degree k in H
    k_sampled = H.degree()
    k_sampled_prop = count2freq([d for n,d in k_sampled]+[0 for _ in range(removed_nodes)])
    K_sampled = 2*m/n
    
    # Prior distribution for G
    k_prop = count2freq(G_deg.values())
    
    assert abs(sum(k_prop.values())-1.0)<0.001, print("oh no")
    assert abs(sum(k_sampled_prop.values())-1.0)<0.001, print("oh no")
    
    # Average degree estimate and real
    K_estimated = 2.0*edges_estimated/N
    K_real = 2.0*M/N
    
    # COMMENT OUT one of these depending on using true prior or poisson
    
#     maxK=0
#     while poisson.cdf(maxK,K_real)<0.9999:
#         maxK+=1
        
    maxK = max(G_deg.values())

    prob_caches={}
    Gp_deg={}
    Gp_deg_crude={}
    for u in H.nodes():
        ki = H.degree(u)
        #den = k_sampled_prop[ki]
        den=0.0
        num = 0.0
        for k in range(ki,maxK+1):
            if (k,ki) not in prob_caches:
                # comment out one of these depending on whether true prior or poisson
#                  prob_caches[(k,ki)]= binom_prob(k,ki,p) * poisson.pmf(k,K_estimated)
                prob_caches[(k,ki)]= binom_prob(k,ki,p) * k_prop[k]
            fact = prob_caches[(k,ki)]
            num+=k*fact
            den+=fact
        if den==0.0:
            Gp_deg[u]=0.0
        else:
            Gp_deg[u]=num/den
        Gp_deg_crude[u]=ki/p
    
    crude_mse = 1/n * sum([(Gp_deg_crude[x]-G_deg[x])**2 for x in Gp_deg_crude])
    post_mse = 1/n * sum([(Gp_deg[x]-G_deg[x])**2 for x in Gp_deg])
    
    return crude_mse, post_mse

# This one uses the measured degree distribution as the evidence
def get_deg_MSE_observed(G_deg,H,p):
    
    # H nodes and edges
    m = len(H.edges())
    n = len(H.nodes())
    
    # G nodes and edges
    M = sum(G_deg.values())/2
    N = len(G_deg)
    
    # Estimated and actual removed edges
    edges_estimated = round(m/p)
    missing_edges = edges_estimated - m
    removed_edges = M - m
    removed_nodes = N - n
    
    # Probability of being degree k in H
    k_sampled = H.degree()
    k_sampled_prop = count2freq([d for n,d in k_sampled]+[0 for _ in range(removed_nodes)])
    K_sampled = 2*m/n
    
    # Prior distribution for G
    k_prop = count2freq(G_deg.values())
    
    assert abs(sum(k_prop.values())-1.0)<0.001, print("oh no")
    assert abs(sum(k_sampled_prop.values())-1.0)<0.001, print("oh no")
    
    # Average degree estimate and real
    K_estimated = 2.0*edges_estimated/N
    K_real = 2.0*M/N
    
    # COMMENT OUT one of these depending on using true prior or poisson
    
#     maxK=0
#     while poisson.cdf(maxK,K_real)<0.9999:
#         maxK+=1
        
    maxK = max(G_deg.values())

    prob_caches={}
    Gp_deg={}
    Gp_deg_crude={}
    for u in H.nodes():
        ki = H.degree(u)
        den = k_sampled_prop[ki]
        num = 0.0
        for k in range(ki,maxK+1):
            if (k,ki) not in prob_caches:
                # comment out depending on whether using true prior or poisson
#                 prob_caches[(k,ki)]= binom_prob(k,ki,p) * poisson.pmf(k,K_estimated)
                prob_caches[(k,ki)]= binom_prob(k,ki,p) * k_prop[k]
            fact = prob_caches[(k,ki)]
            num+=k*fact
        Gp_deg[u]=num/den
        Gp_deg_crude[u]=ki/p
    
    crude_mse = 1/n * sum([(Gp_deg_crude[x]-G_deg[x])**2 for x in Gp_deg_crude])
    post_mse = 1/n * sum([(Gp_deg[x]-G_deg[x])**2 for x in Gp_deg])
    
    return crude_mse, post_mse

G = nx.gnm_random_graph(1000,10000)
G_deg = degree_dict(G)
p=0.1
H = edge_sample(G,p)

print(get_deg_MSE_observed(G_deg, H, p))