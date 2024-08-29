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

def __graph_solution(G:nx.Graph,
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
    

    pos_g = nx.get_node_attributes(G, 'pos')
    if len(pos_g) == 0:
        pos_g = nx.spring_layout(G)
    

    nx.set_node_attributes(G_to_graph,pos_g,'pos')
    
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

def plot_blossom_solution(H,valoresX,U_,Fe):
    __graph_solution(H,len(H.nodes),valoresX,Fe,U_)
    plt.show()

def plot_generical_solution(G:nx.Graph,
                   valoresX:Dict[Tuple[int,int],float],
                   **kwargs):
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
    G_to_graph.add_nodes_from([i for i in range(len(G.nodes))])
    

    pos_g = nx.get_node_attributes(G, 'pos')
    if len(pos_g) == 0:
        pos_g = nx.spring_layout(G)
    

    nx.set_node_attributes(G_to_graph,pos_g,'pos')
    
    for key,value in valoresX.items():
        if value >0:
            G_to_graph.add_edge(key[0],key[1],weight = value)
    
    color_map = []
    edge_color_map = []
    for i in G_to_graph.nodes:
        color_map.append('tab:blue')
        edge_color_map.append('tab:blue')
        
        

    nx.draw_networkx_nodes(
        G_to_graph, 
        pos=pos_g, 
        node_color=color_map,
        node_size=200,
        edgecolors=edge_color_map
    )
    
    nx.draw_networkx_labels(G_to_graph, pos =pos_g, font_size=10)
    delta = 0.0000001
    flag1 = False
    flag2 = False
    for e in G_to_graph.edges:
    
        if G_to_graph[e[0]][e[1]]['weight'] < 1-delta :
            nx.draw_networkx_edges(G_to_graph, pos=pos_g, edgelist=[e], style='dashed',label='0<e<1' if flag1 == False else None)
            flag1 = True
        elif G_to_graph[e[0]][e[1]]['weight'] >= 1-delta:
            nx.draw_networkx_edges(G_to_graph, pos=pos_g, edgelist=[e], style='solid',label='e = 1' if flag2 == False else None)
            flag2 = True
        else:
            raise ValueError('Error en el peso de la arista')
    
    for key,value in kwargs.items():
        if key == 'title':
            plt.title(value)

    plt.legend(loc='upper right')
    plt.show()



