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
        valoresX = model.cbGetSolution(m._xvars)
    elif case2:
        valoresX = model.cbGetNodeRel(m._xvars)

    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
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

    c = {key:min(value,1-value) for key,value in valoresX.items()}

    H = nx.Graph()
    for e in G.edges:
        if valoresX[e] > 0:
            H.add_edge(e[0],e[1],capacity = c[e], weight = valoresX[e])

    pos = nx.get_node_attributes(G,'pos')
    
    GomHu:nx.Graph = nx.gomory_hu_tree(H, capacity = 'capacity')
    nx.set_node_attributes(H, pos, 'pos')
    nx.set_node_attributes(GomHu, pos, 'pos')
    
    
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
            # consolidated_graph(H,GomHu,G,valoresX,U_,Fe,e,save = False)
            # exit(0)
            return (U_,Fe)


    return (None,None)

def grotschel_holland_algorithm(
        G:nx.Graph,
        x_values:Dict[Tuple[int,int],float],
        symmetric:bool = False) -> tuple:
    """
    Grotschel-Holland algorithm for the TSP problem.
    """
    H = nx.Graph()
    for e in G.edges:
        if x_values[e] > 0 and e[0]<e[1]:
            H.add_edge(e[0],e[1], weight = x_values[e])
    aristasN = [e for e in H.edges if totalcap(H,[e])>0]
    GN = nx.Graph()
    GN.add_edges_from(aristasN)

    conN = [a for a in list(nx.connected_components(GN)) if len(a) %2 == 1]
    
    return GN,conN

def grotschell1987(
        G:nx.Graph,
        valoresX:Dict[Tuple[int,int],float]):
    epsilon = 0.3

    H = nx.Graph()
    for e in G.edges:
        if valoresX[e] >= epsilon and valoresX[e] <= 1-epsilon and e[0]<e[1]:
            H.add_edge(e[0],e[1], weight = valoresX[e])
        else:
            H.add_node(e[0])
            H.add_node(e[1])

    conN = [a for a in list(nx.connected_components(H)) if len(a) >= 3]
    if len(conN) == 0:
        return
    print(conN)
    nx.draw(H,pos = nx.get_node_attributes(G,'pos'),with_labels=True,node_size=200)
    plt.show()
    return
    for i,component in enumerate(conN[::-1]):
        print('Componente',i+1,'de',len(conN),':',component)
        sorted_edges = sorted(
            [
                (u, v) for u, v in H.edges(component) #if H[u][v]['weight'] > 1-epsilon
            ], 
            key = lambda e: H[e[0]][e[1]]['weight'],
            reverse = True
        )
        if len(sorted_edges) != 0:

            print(sorted_edges,conN[i])
            print([H[e[0]][e[1]]['weight'] for e in sorted_edges])
            exit(0)
            
    
    pass

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
        model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
        model._subtour_cuts += 1
        # return

    if case2:
        # grotschell1987(model._G,valoresX)

        # GN,conN = grotschel_holland_algorithm(model._G,valoresX)
        # if len(conN) != 0:
        #     for a in conN:
        #         S = list(GN.edges(a))
        #         # print(S,conN)
        #         # graph_solution(model._G,n,valoresX,None,None,None,show_zero = False)
        #         # exit(0)
        #         model.cbLazy(gp.quicksum(model._vars[e] for e in S) <= (len(a)-1))
        #         model._gomhu_cuts += 1
        #     return        


        W,Fe = padberg_rao(model._G,valoresX)
        if W != None and Fe != None:
            model._gomhu_cuts += 1
            # Libro TSP
            # model.cbLazy(gp.quicksum(model._vars[e] for e in [e for e in valoresX if (valoresX[e] > 0) and (e[0] in W) and (e[1] in W) and (e[0]<e[1])] )+ #Suma del manojo
            #              gp.quicksum(model._vars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            
            #Grotschel & holland 1991, restriccion 4.1 
            #COMB INEQUALITY (PADBERG & RAO 1990A FACET IDENTIFICATION FOR THE SYMMETRIC TRAVELING SALESMAN POLYTOPE) PAGE 5
            # model.cbLazy(gp.quicksum(model._vars[e] for e in [e for e in valoresX if (valoresX[e] > 0) and (e[0] in W) and (e[1] in W) and (e[0]<e[1])] )+ 
            #              gp.quicksum(model._vars[e] for e in Fe)<=len(W)+sum([len(e)-1 for e in Fe])-(len(Fe)+1)/2) 
            
            model.cbLazy(gp.quicksum(model._vars[e] for e in [e for e in valoresX if (valoresX[e] > 0) and (e[0] in W) and (e[1] in W) and (e[0]<e[1])] )+ 
                         gp.quicksum(model._vars[e] for e in Fe)<=len(W)+(len(Fe)-1)/2) 

def blossom_method(model:gp.Model, donde):
    n = model._n
    case1 = donde == gp.GRB.Callback.MIPSOL
    case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
    
    if not case1 and not case2:
        return
    
    from common_functions import plot_generical_solution
    if case1:
        valoresX = model.cbGetSolution(model._xvars)

    elif case2:
        valoresX = model.cbGetNodeRel(model._xvars)
    
    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
        model._subtour_cuts += 1

    ############################################################################################################

    if case2:
        delta_value = 0.0000001
        aux_valoresX = {key:value for key,value in valoresX.items() if value >0 and value <1-delta_value}
        fractional_G = nx.Graph()
        fractional_G.add_nodes_from([i for i in range(n)])
        nx.set_node_attributes(fractional_G,nx.get_node_attributes(model._G,'pos'),'pos')

        complete_G = nx.Graph()
        complete_G.add_nodes_from([i for i in range(n)])
        nx.set_node_attributes(complete_G,nx.get_node_attributes(model._G,'pos'),'pos')

        for key,value in aux_valoresX.items():
            fractional_G.add_edge(key[0],key[1],weight = value)
        
        for key,value in valoresX.items():
            if value > 0:
                complete_G.add_edge(key[0],key[1],weight = value)

        connected_components = list(nx.connected_components(fractional_G))
        for component in connected_components:
            
            if len(component)<3 or len(component) % 2 == 0:
                continue
            
            T = []
            fixed_T = []
            out_V = []
            for key,value in valoresX.items():
                condition = (
                    value > 1-delta_value 
                    and (
                        (key[0] in component and key[1] not in component) 
                        or (key[0] not in component and key[1] in component)) 
                    and key[0]<key[1]
                )
                
                if condition:
                    if key[0] not in component:
                        out_V.append(key[0])
                    elif key[1] not in component:
                        out_V.append(key[1])
                    
                    T.append(key)
            
            if len(T)%2 == 0:
                continue

            repeated_nodes_out_V = [i for i in out_V if out_V.count(i) > 1]
            fixed_component = [i for i in component]+list(set(repeated_nodes_out_V))
            fixed_T = [i for i in T if i[0] not in repeated_nodes_out_V and i[1] not in repeated_nodes_out_V]

            if len(fixed_T)== 1:
                # print('fixed component violated a subtour inequality',fixed_component)
                tour2 = [i for i in range(n) if i not in fixed_component]
                model.cbLazy(gp.quicksum(model._xvars[i, j] for i in fixed_component for j in tour2) >= 1)
                model._gomhu_cuts += 1
                return
            elif len(fixed_T) >= 3:
                # print('teeth of a violated blossom inequality with Vi as the handle',fixed_T,fixed_component)
                
                # plot_generical_solution(complete_G, valoresX, title='New G')
                model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in fixed_component and e[1] in fixed_component and e[0]<e[1]] )+ 
                        gp.quicksum(model._xvars[e] for e in fixed_T)<=len(fixed_component)+sum([len(e) for e in fixed_T])-(3*len(fixed_T)+1)/2)
                model._gomhu_cuts += 1
                return
                   
    ############################################################################################################



        W,Fe = padberg_rao(model._G,valoresX)
        if W != None and Fe != None:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                        gp.quicksum(model._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            model._gomhu_cuts += 1
        return


files = os.listdir('Codigos/Data/tsp_experiment_gomhu')

files.remove('.DS_Store')
print(f'callback_name inst obj time gap nodecount blossom subtour')
# file = 'st70.tsp'
# callback = subtourelim_node
for file in files:
    for callback in (subtourelim_node,blossom_method):

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
                for i in range(n) for j in range(i)}

        env = gp.Env(empty=True)
        env.setParam('LogToConsole', 0)
        env.start()

        m = gp.Model(env=env)

        x = m.addVars(dist.keys(), obj=dist, vtype=gp.GRB.BINARY, name='e')
        for i, j in x.keys():   
            x[j, i] = x[i, j]
        m.addConstrs(x.sum(i, '*') == 2 for i in range(n))

        # Optimize model
        G = nx.Graph()
        for i in nodos:
            if len(points) > 0:
                G.add_node(i, pos=points[i])
            else:
                G.add_node(i)
        G.add_edges_from(arcos)


        m._xvars = x
        m._n = n
        m._G = G
        m._gomhu_cuts = 0
        m._subtour_cuts = 0
        m.Params.LazyConstraints = 1
        m.Params.TimeLimit = 600
        m.Params.LogFile = f'/Users/pgutiea/Desktop/TSPJ_MM/Pruebas/experimento/logs/grotschell/symmetric_{instancia.split("/")[-1].split(".")[0]}_{callback.__name__}.log'
        m.update()
        m.optimize(callback)
        try:
            print(f'{callback.__name__} {instancia.split("/")[-1].split(".")[0] } {m.ObjVal} {m.Runtime:.3} {m.MIPGap:.2%} {m.NodeCount} {m._gomhu_cuts} {m._subtour_cuts}')
        except Exception as ex:
            print(f'{callback.__name__} {instancia.split("/")[-1].split(".")[0]}')
            # continue

