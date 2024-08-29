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
        model._subtour_cuts += 1

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
        if Fe is None and H is None:
            color_map.append('tab:blue')
            edge_color_map.append('tab:blue')

        elif i in [e  for subset in Fe for e in subset] and i in H:
            color_map.append('tab:red')
            edge_color_map.append('tab:green')
        elif i in [e  for subset in Fe for e in subset] and i not in H:
            color_map.append('tab:red')
            edge_color_map.append('tab:red')
        elif i not in [e  for subset in Fe for e in subset] and i not in H:
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

    c = {key:min(value+valoresX[(key[1],key[0])],1-value-valoresX[(key[1],key[0])]) for key,value in valoresX.items()}

    H = nx.Graph()
    for e in G.edges:
        if valoresX[e] > 0 and e[0]<e[1]:
            H.add_edge(e[0],e[1],capacity = c[e], weight = valoresX[e]+valoresX[(e[1],e[0])])

    GomHu:nx.Graph = nx.gomory_hu_tree(H, capacity = 'capacity')
    
    for e in GomHu.edges:
        Te = nx.Graph(GomHu)
        Te.remove_edge(*e)
        U,V = list(nx.connected_components(Te))
        
        if len(U)<len(V):
            U_,V_ = U,V
        else:
            V_,U_ = U,V
        
        if len(U_)<3:
            continue

        delta_W = delta(H,U_)
        Fe = []
        nodos_agregados = []
        for e in delta(H,U_):
            if 1 - H[e[0]][e[1]]['weight'] < H[e[0]][e[1]]['weight'] :
                if e[0] in nodos_agregados or e[1] in nodos_agregados:
                    continue
                Fe.append(e)
                nodos_agregados.append(e[0])
                nodos_agregados.append(e[1])
        Fe = set(Fe)

        if len(Fe) == 0:
            continue

        if len(Fe)%2 == 0:
            Faux = sorted([(term(edge,H),edge) for edge in delta_W])
            fp = set()

            fp.add(Faux[0][1])
            Fe = Fe.symmetric_difference(fp)
        
        des = totalcap(H,delta_W.difference(Fe)) + len(Fe) - totalcap(H,Fe)

        if len(Fe) != len(U_):
            continue
        if des < 0.999:
            # consolidated_graph(H,GomHu,G,valoresX,U_,Fe,e,save = True)
            return (U_,Fe)

    return (None,None)

def subtourelim_gomhu(model:gp.Model, donde):
    n = model._n
    case1 = donde == gp.GRB.Callback.MIPSOL
    case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
    
    if not case1 and not case2:
        return
    
    if case1:
        valoresX = model.cbGetSolution(m._vars)
    elif case2:
        valoresX = model.cbGetNodeRel(m._vars)

    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._vars[i, j] for i in tour for j in tour2) >= 1)
        model._subtour_cuts += 1

    if case2:
        W,Fe = padberg_rao(model._G,valoresX)
        if W != None and Fe != None:
            model._gomhu_cuts += 1
            model.cbLazy(gp.quicksum(model._vars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ #Suma del manojo
                         gp.quicksum(model._vars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)

files = os.listdir('Codigos/Data/tsp_experiment_gomhu')
files.remove('.DS_Store')
for callback in (subtourelim_node,subtourelim_gomhu):
    for file in files:
        instancia = f'Codigos/Data/tsp_experiment_gomhu/{file}'
        problem = tsplib95.load(instancia)
        if problem.edge_weight_type != 'EUC_2D':
            continue
        n = len(list(problem.get_nodes()))

        points = [(problem.node_coords[i][0],problem.node_coords[i][1]) for i in range(1,n+1)]
        nodos = [i for i in range(n)]
        arcos = [(i,j) for i in nodos for j in nodos if i !=j]
        dist = {(i, j):
                round(math.sqrt(sum((points[i][k]-points[j][k])**2 for k in range(2))),0)
                for i in range(n) for j in range(n) if i != j}

        env = gp.Env(empty=True)
        env.setParam('LogToConsole', 0)
        env.start()

        m = gp.Model(env=env)
        x = m.addVars(dist.keys(), obj=dist, vtype=gp.GRB.BINARY, name='e')

        m.addConstrs((x.sum(i, '*') == 1 for i in range(n)), name='out')
        m.addConstrs((x.sum('*', i) == 1 for i in range(n)), name='in')

        # Optimize model
        G = nx.Graph()
        for i in nodos:
            G.add_node(i)#, pos = (points[i][0],points[i][1]))
        G.add_edges_from(arcos)


        m._vars = x
        m._n = n
        m._G = G
        m._gomhu_cuts = 0
        m._subtour_cuts = 0
        m.Params.LazyConstraints = 1
        m.Params.TimeLimit = 600
        m.Params.LogFile = f'/Users/pgutiea/Desktop/TSPJ_MM/Pruebas/experimento/logs/asymmetric_{instancia.split("/")[-1].split(".")[0]}_{callback.__name__}.log'
        m.update()
        # callback = subtourelim_gomhu
        m.optimize(callback)
        try:
            print(f'{callback.__name__} {instancia.split("/")[-1].split(".")[0] } {m.ObjVal} {m.Runtime:.3} {m.MIPGap:.2%} {m.NodeCount} {m._gomhu_cuts} {m._subtour_cuts}')
        except Exception as ex:
            print(f'{callback.__name__} {instancia.split("/")[-1].split(".")[0]}')
            continue

