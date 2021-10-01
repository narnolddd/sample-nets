import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
import math

from collections import Counter
a4_dims = (11.69,8.27)
from scipy.stats import poisson, zipf

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

def count2freq(arr):
    counter = Counter(arr)
    tot = sum(counter.values())
    for k in counter:
        counter[k]/=tot
    return counter

def triangles(n1, n2, G):
    if not G.has_edge(n1,n2):
        return 0
    else:
        return len(list(nx.common_neighbors(G,n1,n2)))
    
def edge_triangle_count(G):
    tc = []
    for i in G.nodes():
        for j in G.neighbors(i):
            if i>j:
                continue
            tc.append(triangles(i,j,G))
    return tc

def edge_triangle_dict(G):
    tc={}
    for u,v in G.edges():
        tc[(min(u,v),max(u,v))]=triangles(u,v,G)
    return tc

# this one is for LOTP method
def triangle_likelihood(t,tc,p):
    return math.comb(t,tc)*math.pow((1.0-p*p),(t))

# for the empirical evidence method
def triangle_likelihood_normalised(t,tc,p):
    return math.comb(t,tc)*math.pow(p,2.0*tc)*math.pow(1.0-p*p,t-tc)

# For the empirical evidence method
def get_tri_MSE_observed(G_tri,H,p):
    m = len(H.edges())
    
    edges_estimated = round(m/p)
    missing_edges = edges_estimated - len(H.edges())
    removed_edges = len(G.edges()) - len(H.edges())
    
    tl_sampled = edge_triangle_dict(H)
    tl_sampled_prop = count2freq(list(tl_sampled.values())+[0 for _ in range(missing_edges)])
    T_sampled = sum(tl_sampled.values())/3.0
    
    # simple p3 estimate
    T_estimated = T_sampled/p**3
    T_real = sum(G_tri.values())/3
    
    T_lambda = T_estimated*3.0/(len(H.edges())+missing_edges)
    wedge_lambda = 2*wedge_count(H)/(len(H.edges())+missing_edges)
    
    # Comment out as appropriate depending on Poisson or true prior
    
    maxT=0
    while poisson.cdf(maxT,T_lambda)<0.9999:
        maxT+=1
        
#    maxT = max(G_tri.values())
    
    tl_prop = count2freq(G_tri.values())
    assert abs(sum(tl_prop.values()) - 1.0) < 0.0001, print("oh no")
    
    prob_caches={}
    Gp_tri={}
    Gp_tri_crude={}
    for u,v in H.edges():
        el = tl_sampled[(min(u,v),max(u,v))]
        den = tl_sampled_prop[el]
        num = 0.0
        for t in range(el,maxT+1):
            if (t,el) not in prob_caches:
                prob_caches[(t,el)]= triangle_likelihood_normalised(t,el,p) * poisson.pmf(t,T_lambda)
#                prob_caches[(t,el)]= triangle_likelihood_normalised(t,el,p) * tl_prop[t]
            fact = prob_caches[(t,el)]
            num+=t*fact
        Gp_tri[(min(u,v),max(u,v))]=num/den
        Gp_tri_crude[(min(u,v),max(u,v))]=el/p**3
        
    num=0.0
    el = 0
    den=tl_sampled_prop[el]
    for t in range(el,maxT+1):
        fact = triangle_likelihood_normalised(t,el,p) * poisson.pmf(t,T_lambda)
#        fact = triangle_likelihood_normalised(t,el,p) * tl_prop[t]
        num+=t*fact
    missing_triangles=removed_edges*(num/den)
    
    crude_mse = 1/m * sum([(Gp_tri_crude[x]-G_tri[x])**2 for x in Gp_tri_crude])
    post_mse = 1/m * sum([(Gp_tri[x]-G_tri[x])**2 for x in Gp_tri])
    post_T = (sum(Gp_tri.values())+missing_triangles)/3
    
    return crude_mse, post_mse, abs(T_estimated - T_real), abs(post_T - T_real)

def get_tri_MSE_LOTP(G_tri,H,p):
    m = len(H.edges())
    
    edges_estimated = round(m/p)
    missing_edges = edges_estimated - len(H.edges())
    removed_edges = len(G.edges()) - len(H.edges())
    
    tl_sampled = edge_triangle_dict(H)
    tl_sampled_prop = count2freq(list(tl_sampled.values())+[0 for _ in range(removed_edges)])
    T_sampled = sum(tl_sampled.values())/3.0
    
    # simple p3 estimate
    T_estimated = T_sampled/p**3
    T_real = sum(G_tri.values())/3
    
    T_lambda = T_estimated*3.0/(len(H.edges())+removed_edges)
    
    # Comment out appropriately
    
#     maxT=0
#     while poisson.cdf(maxT,T_lambda)<0.9999:
#         maxT+=1
    
    maxT = max(G_tri.values())
    
    tl_prop = count2freq(G_tri.values())
    
    prob_caches={}
    Gp_tri={}
    Gp_tri_crude={}
    for u,v in H.edges():
        el = tl_sampled[(min(u,v),max(u,v))]
        den = 0.0
        num = 0.0
        for t in range(el,maxT+1):
            if (t,el) not in prob_caches:
                #comment out as appropriate
                #prob_caches[(t,el)]= triangle_likelihood(t,el,p) * poisson.pmf(t,T_lambda)
                prob_caches[(t,el)]= triangle_likelihood(t,el,p) * tl_prop[t]
            fact = prob_caches[(t,el)]
            num+=t*fact
            den+=fact
        if den!=0.0:
            Gp_tri[(min(u,v),max(u,v))]=num/den
        else:
            Gp_tri[(min(u,v),max(u,v))]=0.0
        Gp_tri_crude[(min(u,v),max(u,v))]=el/p**3
        
    num=0.0
    el = 0
    den=0.0
    for t in range(el,maxT+1):
        if (t,el) not in prob_caches:
            #comment out as appropriate
            #prob_caches[(t,el)] = triangle_likelihood(t,el,p) * poisson.pmf(t,T_lambda)
            prob_caches[(t,el)]= triangle_likelihood(t,el,p) * tl_prop[t]
        fact = prob_caches[(t,el)]
        num+=t*fact
        den+=fact
    missing_triangles=removed_edges*(num/den)
    
    crude_mse = 1/m * sum([(Gp_tri_crude[x]-G_tri[x])**2 for x in Gp_tri_crude])
    post_mse = 1/m * sum([(Gp_tri[x]-G_tri[x])**2 for x in Gp_tri])
    post_T = (sum(Gp_tri.values())+missing_triangles)/3
    
    return crude_mse, post_mse, abs(T_estimated - T_real), abs(post_T - T_real)

n, m, p, seed = 1000, 25000, 0.5, 101

G = nx.gnm_random_graph(n, m, seed)
#G = nx.barabasi_albert_graph(n,10,seed)
H = edge_sample(G,p)

G_tri = edge_triangle_dict(G)
print(get_tri_MSE_LOTP(G_tri,H,p))