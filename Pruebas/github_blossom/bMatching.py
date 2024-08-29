import networkx as nx
from gomory_hu import gomory_hu_tree, min_odd_cut
import numpy as np
from itertools import combinations
from gurobipy import *
import time
import math
import tsplib95

from blossom_separation import blossom_separation

def all_subsets(arr):
    """returns a list of all subsets of a list"""
    combs = []
    for i in xrange(0, len(arr)+1):
        listing = [list(x) for x in combinations(arr, i)]
        combs.extend(listing)
    return combs

def check_blossom(model, G, b, a):
    X = {}
    for key in a:
        X[key] = a[key].getAttr("x")
    thisG = G.copy()
    cut = blossom_separation(thisG, b, X)
    if cut == None:
        return False
    else:
        Ew = []
        bw = sum([b[i] for i in cut])
        for e in G.edges:
            if e[0] in cut and e[1] in cut:
                Ew.append(frozenset(e))
        model.addConstr(quicksum([a[e] for e in Ew]) <= (bw-1)/2 )
        print("===========CUT!=============")
        return True

def bMatchingInteger(G, b, r):
    '''
    input:  G = (V, E)
            b = integer capacity vector
    output: max cardinality matching
    '''

    IP = Model()
    IP.setParam('OutputFlag', 0)
    V = G.nodes
    E = [frozenset(e) for e in G.edges]

    a = IP.addVars(E, name = "a", vtype=GRB.INTEGER)

    obj = 0
    for e in E:
        obj += r[e] * a[e]

    IP.setObjective(obj, GRB.MAXIMIZE)

    for i in V:
        Ei = [frozenset(e) for e in G.edges(i)]
        IP.addConstr(quicksum([a[e] for e in Ei]) <= b[i])
    IP.optimize()
    action = {}
    for e in E:
        action[e] = abs(a[e].X)
    return action

def bMatching(G, b, r):
    '''
    input:  G = (V, E)
            b = integer capacity vector
    output: max cardinality matching
    '''

    LP = Model()
    LP.setParam('OutputFlag', 0)
    V = G.nodes
    E = [frozenset(e) for e in G.edges]

    a = LP.addVars(E, name = "a")

    obj = 0
    for e in E:
        obj += r[e] * a[e]
    LP.setObjective(obj, GRB.MAXIMIZE)

    for i in V:
        Ei = [frozenset(e) for e in G.edges(i)]
        LP.addConstr(quicksum([a[e] for e in Ei]) <= b[i])
    LP.optimize()
    while check_blossom(LP, G, b, a):
        LP.setParam('Method', 1)
        LP.optimize()
    action = {}
    for e in E:
        action[e] = abs(a[e].X)
    return action.values()

def bMatchingComplete(G, b):
    '''
    input:  G = (V, E)
            b = integer capacity vector
    output: max cardinality matching
    '''

    LP = Model()
    LP.setParam('OutputFlag', 0)
    V = G.nodes
    E = [frozenset(e) for e in G.edges]

    a = LP.addVars(E, name = "a")

    obj = quicksum(a)
    LP.setObjective(obj, GRB.MAXIMIZE)

    for i in V:
        Ei = [frozenset(e) for e in G.edges(i)]
        LP.addConstr(quicksum([a[e] for e in Ei]) <= b[i])

    for W in all_subsets(list(V)):
        Ew = []
        bw = sum([b[i] for i in W])
        for e in G.edges:
            if e[0] in W and e[1] in W:
                Ew.append(frozenset(e))
        if len(Ew)>0:
            LP.addConstr(quicksum([a[e] for e in Ew]) <= (bw-1)/2 )
    LP.optimize()

    print(LP.objVal)


def mycallback(model, where):
    thisG = G.copy()
    if where == GRB.Callback.MIPNODE:
        status = model.cbGet(GRB.Callback.MIPNODE_STATUS)
        if status == GRB.OPTIMAL:
            rel = model.cbGetNodeRel(model._vars)
            cut = blossom_separation(thisG, b, rel)
            if cut is None:
                return
            print('hola')
            Ew = []
            bw = sum([b[i] for i in cut])
            for e in thisG.edges:
                if e[0] in cut and e[1] in cut:
                    Ew.append(frozenset(e))
            model.cbCut(quicksum([model._vars[e] for e in Ew]) <= (bw-1)/2 )


def bMatchingUserCut(G, b,r):
    '''
    input:  G = (V, E)
            b = integer capacity vector
    output: max cardinality matching
    '''

    IP = Model()
    IP.setParam('OutputFlag', 0)
    V = G.nodes
    E = [frozenset(e) for e in G.edges]

    a = IP.addVars(E, name = "a", vtype=GRB.BINARY)

    obj = 0
    for e in E:
        obj += r[e] * a[e]
    IP.setObjective(obj, GRB.MINIMIZE)
    # IP.setParam('Method', 1)
    for i in V:
        Ei = [frozenset(e) for e in G.edges(i)]
        IP.addConstr(quicksum([a[e] for e in Ei]) == b[i])

    IP._vars = a

    IP.optimize(mycallback)

    return IP

G = nx.complete_graph(100)
G.remove_edges_from(nx.selfloop_edges(G))

b = {}
E = [frozenset(e) for e in G.edges]
for i in G.nodes():
    b[i] = 2
r = {}
for e in E:
    if (list(e)[0], list(e)[1]) in r.keys():
        r[e] = r[(list(e)[0], list(e)[1])]
    else:
        r[e] = np.random.rand()



P = bMatchingUserCut(G,b,r)



H = nx.Graph()
variables = P.getVars()
print(len(variables))
P.write("blossom.lp")
for v in variables:
    if v.X > 0:
        u = int(v.VarName.split(",")[0].split('[')[1] )
        _v = int(v.VarName.split(",")[1].split(']')[0])
        # print(u,_v,':',v.X)
        H.add_edge(u,_v)
nx.draw(H, with_labels=True)
import matplotlib.pyplot as plt
plt.show()