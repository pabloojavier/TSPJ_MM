import math
import random
from itertools import combinations
import gurobipy as gp
import tsplib95
import networkx as nx
import numpy as np
import csv
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple,Set
import os



def subtour_method_sol(subruta, vals,n):
    # obtener una lista con los arcos parte de la solucións
    arcos = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)
    noVisitados = list(range(n))
    while noVisitados: # true if list is non-empty
        ciclo = []
        vecinos = noVisitados
        while vecinos:
            actual = vecinos[0]
            ciclo.append(actual)
            noVisitados.remove(actual)
            vecinos = [j for i, j in arcos.select(actual, '*') if j in noVisitados]
        if len(subruta) > len(ciclo):
            subruta[:] = ciclo

def subtourelim_sol(model:gp.Model, donde):
    
    if donde == gp.GRB.Callback.MIPSOL:
        vals = model.cbGetSolution(model._vars)
        # for i in vals:
        #     if vals[i] > 0 and vals[i] < 1:
        #         print(i,vals[i])
        # find the shortest cycle in the selected edge list
        tour = [i for i in range(model._n+1)]
        subtour_method_sol(tour, vals,model._n)

        if len(tour) < model._n:
            # tour2 = [i for i in range(n) if i not in tour]
            model.cbLazy(gp.quicksum(model._vars[i, j]  for i, j in combinations(tour, 2) ) <= len(tour)-1)

def subtour_method_node(subruta, vals,n):
    arcos = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)
    noVisitados = list(range(n))
    while noVisitados: # true if list is non-empty
        ciclo = []
        vecinos = noVisitados
        while vecinos:
            actual = vecinos[0]
            ciclo.append(actual)
            noVisitados.remove(actual)
            vecinos = [j for i, j in arcos.select(actual, '*') if j in noVisitados]
        if len(subruta) > len(ciclo):
            subruta[:] = ciclo

def subtourelim_node(model:gp.Model, donde):
    n = model._n
    case1 = donde == gp.GRB.Callback.MIPSOL
    case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
    
    if not case1 and not case2:
        return
    
    # retrieve the LP solution
    if case1:
        valoresX = model.cbGetSolution(m._vars)
    elif case2:
        valoresX = model.cbGetNodeRel(m._vars)

    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._vars[i, j] for i in tour for j in tour2) >= 1)

def delta(G,U): 
    """dado un conjunto U de vertices y una grafo G que los contenga, nos devuelve las
aristas de G que tienen un unico extremo en U."""
    ejes = set()
    for l, nbrs in ((n, G[n]) for n in U):
        ejes.update( (l,y) if l<y else (y,l) for y in nbrs if y not in U)
    return ejes

def term(e,H):
    s = max(H[e[0]][e[1]]['weight'], 1-H[e[0]][e[1]]['weight'])
    t = min(H[e[0]][e[1]]['weight'], 1-H[e[0]][e[1]]['weight']) 
    return(s-t)

def totalcap(G,F,capacity = 'weight'):
    return sum([G[e[0]][e[1]][capacity] for e in F])

def check_blossom(handle:Set[tuple],T:List[Set[int]]):
    'Check the three conditions for handle and teeths'
    T = [set(t) for t in T]
    s = len(T)
    if handle is None or T is None:
        return False
    for j in range(s):
        if (not len(handle.intersection(T[j])) >= 1) or (not len(T[j]-handle) >= 1):
            return False
        for i in range(j):
            if not T[i].isdisjoint(T[j]):
                return False
    return True

def graph_solution(G:nx.Graph,
                   n:int,
                   valoresX:Dict[Tuple[int,int],float],
                   Fe : Set[Tuple[int]] = None,
                   H : List[int] = None,
                   ax : plt.Axes = None,
                   show_zero : bool = False):
    """
    Plot a figure with the graph G and the solution x.

    The edges from x are plotted with different styles depending on their value:
        - solid line if the edge has a value of 1
        - dotted line if the edge has a value of 0
        - dashed line if the edge has a value between 0 and 1

    Input:
        -G: Original graph, with the nodes and the positions
        -n: Number of nodes
        -valoresX: Vector of solution x
        -Fe: Subset of teeth. A set with tuple with the nodes of each tooth. This subset is colored in red.
        -ax: Axes of the plot, if we want to save the plot
        -show_zero: If True, the edges with value 0 are plotted. If False, they are not plotted

    """
    G_to_graph = nx.Graph()
    G_to_graph.add_nodes_from([i for i in range(n)])
    nx.set_node_attributes(G_to_graph,nx.get_node_attributes(G,'pos'),'pos')
    pos_g = nx.get_node_attributes(G_to_graph,'pos')
    
    for key,value in valoresX.items():
        if show_zero == True:
            G_to_graph.add_edge(key[0],key[1],weight = value)
        elif show_zero == False and value >0:
            G_to_graph.add_edge(key[0],key[1],weight = value)
    
    color_map = []
    edge_color_map = []
    for i in G_to_graph.nodes:
        if Fe is not None and i in [e  for subset in Fe for e in subset] and i in H:
            color_map.append('tab:red')
            edge_color_map.append('tab:green')
        elif Fe is not None and i in [e  for subset in Fe for e in subset] and i not in H:
            color_map.append('tab:red')
            edge_color_map.append('tab:red')
        elif Fe is not None and i not in [e  for subset in Fe for e in subset] and i not in H:
            color_map.append('tab:blue')
            edge_color_map.append('tab:blue')
        
        elif Fe is None and H is None:
            color_map.append('tab:blue')
            edge_color_map.append('tab:blue')
        
    try:
        nx.draw_networkx_nodes(G_to_graph, pos=pos_g, node_color=color_map,node_size=200,ax=ax,edgecolors=edge_color_map)
    except:
        nx.draw_networkx_nodes(G_to_graph, pos=pos_g,node_size=200,ax=ax)
    nx.draw_networkx_labels(G_to_graph, pos =pos_g, font_size=10,ax=ax)
    for e in G_to_graph.edges:
        if G_to_graph[e[0]][e[1]]['weight'] <= 0.001:
            nx.draw_networkx_edges(G_to_graph, pos=pos_g, edgelist=[e], style='dotted',ax=ax,alpha = 0.3)
        elif G_to_graph[e[0]][e[1]]['weight'] < 1 :
            nx.draw_networkx_edges(G_to_graph, pos=pos_g, edgelist=[e], style='dashed',ax=ax)
        elif G_to_graph[e[0]][e[1]]['weight'] > 0.99999:
            nx.draw_networkx_edges(G_to_graph, pos=pos_g, edgelist=[e], style='solid',ax=ax)
        else:
            raise ValueError('Error en el peso de la arista')
    
    if ax is None:
        plt.show()

def consolidated_graph(H,GomHu,G,valoresX,U_,Fe,e,save = False):
    fig,ax = plt.subplots(1,2,figsize=(12,6))
    graph_solution(H,len(G.nodes),valoresX,Fe,U_,ax[0],show_zero = False)
    graph_solution(GomHu,len(G.nodes),nx.get_edge_attributes(GomHu,'weight'),None,None,ax[1],show_zero = True)
    ax[0].set_title('Graph with vector x*',y=0)
    ax[1].set_title('Gomory Hu with weights 1-x*',y=0)
    instance_name = instancia.split('/')[-1].split('.')[0]
    plt.suptitle(instance_name+f'\n U:{U_}\n Fe:{Fe}')
    if save:
        file_path = f'cut_plots/{instance_name}'
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        plt.savefig(f'{file_path}/{instance_name}_{m._cuts_added}_cuts_{e[0]}_{e[1]}.png',dpi=300)
        plt.close()
    else:
        plt.show()

def padberg_rao(G:nx.Graph,
                valoresX:Dict[Tuple[int,int],float]
                ) -> tuple:

    c = {key:min(value,1-value) for key,value in valoresX.items()}

    H = nx.Graph()
    H_prime = nx.Graph()
    for e in G.edges:
        if valoresX[e] > 0:
            H_prime.add_edge(e[0],e[1],capacity = c[e], weight = valoresX[e])
        H.add_edge(e[0],e[1],capacity = c[e], weight = valoresX[e])

    nx.set_node_attributes(H,nx.get_node_attributes(G,'pos'),'pos')

    GomHu:nx.Graph = nx.gomory_hu_tree(H_prime, capacity = 'capacity')
    nx.set_node_attributes(GomHu,nx.get_node_attributes(G,'pos'),'pos')
    
    for e in GomHu.edges:
        Te = nx.Graph(GomHu)
        Te.remove_edge(*e)
        U,V = list(nx.connected_components(Te))
        
        if len(U)<len(V):
            U_,V_ = U,V
        else:
            V_,U_ = U,V

        delta_W = delta(H_prime,U_)
        Fe = []
        for e in delta(H_prime,U_):
            if 1 - H_prime[e[0]][e[1]]['weight'] < H_prime[e[0]][e[1]]['weight'] :
                Fe.append(e)
        Fe = set(Fe)
        # Fe = set([e for e in delta(H_prime,U_) if 1 - H[e[0]][e[1]]['weight'] < H[e[0]][e[1]]['weight']])

        if len(Fe)%2 == 0:
            Faux = sorted([(term(edge,H_prime),edge) for edge in delta_W])
            fp = set()
            fp.add(Faux[0][1])
            Fe = Fe.symmetric_difference(fp)
        
        des = totalcap(H_prime,delta_W.difference(Fe)) + len(Fe) - totalcap(H_prime,Fe)

        # aux_list = [e for edge in Fe for e in edge]
        # if len(Fe) != len(U_) or len(aux_list) != len(set(aux_list)):
        #     continue
        
        if des < 0.999:# and check_blossom(U_,Fe):
            # print([H_prime[e[0]][e[1]]['weight'] for e in Fe])
            # print(e,U_,Fe,delta(H_prime,U_),f'{des}=({totalcap(H_prime,delta_W.difference(Fe))}+{len(Fe)})-{totalcap(H_prime,Fe)}')
            # consolidated_graph(H,GomHu,G,valoresX,U_,Fe,e,save = True)
            m._cuts_added += 1
            # if m._cuts_added == 4:
            #     exit(0)
            return (U_,Fe,delta_W)


    return (None,None,None)

def delta_Fe(x,Fe):
    # def delta_W(x,u):
    #     'Dada una solución x, entrega todas las aristas incidentes en manojo U con un solo extremo en U'
    #     return set([(i,j) if i < j else (j,i) for i,j in x.keys() if x[i,j] > 0 and ((i in u and j not in u) or (i not in u and j in u))])
    'Dada una solución x, entrega todas las aristas incidentes en los dientes Fe'
    conjunto = []
    for diente in Fe:
        conjunto.append(delta_W(x,diente))  
    return [e  for subset in conjunto for e in subset]

def cutting_planes(model:gp.Model):
    
    #valoresX =     

    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._vars[i, j] for i in tour for j in tour2) >= 1)


    W,Fe,set_delta_W = padberg_rao(model._G,valoresX)
    if W != None and Fe != None:
        # print([e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]],Fe,len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2   )
        model.cbLazy(gp.quicksum(model._vars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ #Suma del manojo
                        gp.quicksum(model._vars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)




        
instancia = 'Codigos/Data/tsp/eil51.tsp'
problem = tsplib95.load(instancia)
n = len(list(problem.get_nodes()))

points = [(problem.node_coords[i][0],problem.node_coords[i][1]) for i in range(1,n+1)]
nodos = [i for i in range(n)]
arcos = [(i,j) for i in nodos for j in nodos if i !=j]
dist = {(i, j):
        round(math.sqrt(sum((points[i][k]-points[j][k])**2 for k in range(2))),0)
        for i in range(n) for j in range(i)}



env = gp.Env(empty=True)
env.setParam('OutputFlag', 0)
env.start()

m = gp.Model(env=env)

x = m.addVars(dist.keys(), obj=dist, vtype=gp.GRB.CONTINUOUS,lb = 0,ub = 1, name='e')
for i, j in x.keys():   
    x[j, i] = x[i, j]
m.addConstrs(x.sum(i, '*') == 2 for i in range(n))

G = nx.Graph()
for i in nodos:
    G.add_node(i, pos = (points[i][0],points[i][1]))
G.add_edges_from(arcos)

tour = [0]


while len(tour) < n:
    
    m.optimize()
    
    m.update()
    tour = [i for i in range(n+1)]
    valoresX = m.getAttr('X', x)
    subtour_method_sol(tour, valoresX,n)

    m.addConstr(gp.quicksum(x[i, j]  for i, j in combinations(tour, 2) ) <= len(tour)-1)
    
    print(len(tour),n,tour)
    # m.write('modelo.lp')
    # print('Callback:',callback.__name__)
print({(i,j):valoresX[i,j] for i,j in valoresX if valoresX[i,j] !=0})
print(tour)
print('Optimal cost: %g' % m.ObjVal)
print('Time: %g' % m.Runtime)


sol = nx.Graph()
sol.add_nodes_from(G)
valoresX = m.getAttr('X', x)
for i,j in valoresX:
    if valoresX[i,j] > 0:
        sol.add_edge(i,j,weight = valoresX[i,j])
nx.draw(sol,with_labels=True,pos = nx.get_node_attributes(G,'pos'))
G = nx.Graph()
for i in nodos:
    G.add_node(i, pos = (points[i][0],points[i][1]))

for i,j in arcos:
    if x[i,j].x > 0.5:
        G.add_edge(i,j)

nx.draw(G,with_labels=True,pos=nx.get_node_attributes(G,'pos'))
plt.show()
