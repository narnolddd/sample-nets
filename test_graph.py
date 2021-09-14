#!/usr/bin/env python3
import sys
import math
import random


import matplotlib.pyplot as plt
import networkx as nx

def triangles(n1, n2, G):
    if not G.has_edge(n1,n2):
        return 0
    else:
        return len(list(nx.common_neighbors(G,n1,n2)))

n = 1000
m = 10000
seed = 20160  # seed random number generators for reproducibility
p=0.5

## Use seed for reproducibility
#G = nx.gnm_random_graph(n, m, seed=seed)
G = nx.gnm_random_graph(n, m)
num_tri = sum(nx.triangles(G).values()) / 3
print("Tri:",num_tri)
removeList=[]
tCount={}
for u,v,a in G.edges(data=True):
    if (random.random()) >= p:
        removeList.append((u,v))
    tc= triangles(u,v,G)
    if tc in tCount:
        tCount[tc]= tCount[tc]+1.0
    else:
        tCount[tc]= 1.0

maxT=0
for t in tCount:
    tCount[t]= tCount[t]/len(G.edges)
    if t > maxT:
        maxT= t


#print("Sum=",sumt,tCount)
G2= G.copy()
for e in removeList:
    G2.remove_edge(e[0],e[1])
tCount2={}
for u,v,a in G.edges(data=True):
    tc= triangles(u,v,G2)
    if tc in tCount2:
        tCount2[tc]= tCount2[tc]+1.0
    else:
        tCount2[tc]= 1.0
for t in tCount2:
    tCount2[t]= tCount2[t]/len(G.edges)
print(tCount2)
num_tri = sum(nx.triangles(G2).values()) / 3

#print("Tri maxT:",num_tri, maxT)

ElTot= 0.0
zcount=0.0
test= 0.0
tsampest= 0.0
for t in tCount:
    test+= t*tCount[t]*len(G.edges)
    try:
        tsampest+= t*tCount2[t]*len(G.edges)
    except:
        pass
print ("Test = ",test/3.0," Tsampest=",tsampest/3.0)
print(tCount)
for u,v,a in G.edges(data=True):
    tc= triangles(u,v,G2)
    toc= triangles(u,v,G)
    zcount+=math.pow((1.0-p*p*p),toc)
    den=0.0
    num=0.0
    for t in range(tc,maxT+1):
        fact= math.comb(t,tc)*math.pow(p*p,(tc))*math.pow((1.0-p*p),(t-tc))*tCount[t]
        num+= t*fact
        den=tCount2[tc]
        # if tc > 0:
            # fact= math.comb(t,tc)*math.pow(p*p,(tc))*math.pow((1.0-p*p),(t-tc))*tCount[t]
            # num+= t*fact
            # den=tCount2[tc]
        # else:
            # num+= t*((1.0-p)+p*math.pow(1.0-p*p,t))*tCount[t]
            # den=tCount2[0]

    ElTot+=num/den
print("Zcount=",zcount,"Zprop=",zcount/len(G.edges))
print("ElTot=", ElTot/3)
print("Tri /p^3=", num_tri/(p*p*p))
