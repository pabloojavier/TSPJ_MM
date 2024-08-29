import sys
from itertools import combinations
import gurobipy as gp
from gurobipy import *
import tsplib95
import matplotlib.pyplot as plt
import networkx as nx

import numpy as np

# Definimos algunas funciones para ayudarnos a escribir las restricciones y la funcion objetivo.
def suma(S):
    s = 0
    for (a,b) in S:
        if a<b:
            s = s + x[a, b]
        else:
            s = s + x[b,a]
    return(s)

# Definimos la funcion delta , que dado un conjunto U de vertices y una grafo G que los contenga , nos devuelve las aristas de G que tienen un unico extremo en U.
def delta(G,U):
    ejes = set()
    for l, nbrs in ((n, G[n]) for n in U):
        ejes.update ((l, y) if l<y else (y,l) for y in nbrs if y not in U)
    return(ejes)

def totalcap(G,F, capacity = 'weight'):
    cap = [G[e[0]][e[1]][capacity] for e in F]
    val = sum(cap)
    return(val)

                
def Padberg_Rao(G, x):
    # Para aplicar el algoritmo , consideramos el vector de pesos
    c = np.array ([min(a, 1 - a) for a in x])
    # Construimos el nuevo grafo H sobre G con capacidad c, atributo con el que calculamos el arbol de Gomory−Hu, y con peso x.
    H = nx.Graph ()
    for i, e in enumerate(G.edges ()):
        H.add_edge(e[0],e[1],capacity = c[i], weight = x[i])

    # A continuacion, consideramos el arbol de Gomory−Hu para el grafo H
    GomHu = nx.gomory_hu_tree (H, capacity = 'capacity')
    #nx.draw(GomHu, pos=my_pos, node_size=15, width=0.5)
    #plt.show()
    # Consideramos ahora el conjunto de corte asociado a cada arista del arbol y para cada uno de ellos buscamos una desigualdad deseada.
    for e in GomHu.edges():
        # Asignamos una copia de GomHu a la cual eliminamos la  arista y calculamos las componentes conexas
        Te = nx.Graph(GomHu)
        Te.remove_edge(*e)
        U,V = list(nx.connected_components(Te))
        # Calculamos las arista del conjunto de cortes
        cutset = delta(H,U)
        #print(cutset)
        # Encontramos ahora el conjunto F con las propiedades deseadas
        Fe = set([e for e in cutset if 1 - H[e[0]][e[1]]['weight'] < H[e[0]][e[1]]['weight']])
        # Si el conjunto anterior verifica la propiedad de que la suma de los cardinales sea impar , es el optimo. En otro caso obtenemos el optimo mediante el siguiente metodo
        if ((len(U) + len(Fe)) % 2 == 0):
            def term(e):
                s = max(H[e[0]][e[1]]['weight'], 1 - H[e[0]][e[1]]['weight'])
                t = min(H[e[0]][e[1]]['weight'], 1 - H[e[0]][e[1]]['weight'])
                return(s - t)
            Faux = sorted([(term(e),e) for e in cutset])
            fp = set()
            fp.add(Faux[0][1])
            Fe = Fe.symmetric_difference(fp)
            
        # Finalmente, comprobamos si se viola la desigualdad y, en dicho caso , hemos acabado y devolvemos.
        des = totalcap(H, cutset.difference (Fe)) + len(Fe) - totalcap(H,Fe)
        if des <1:
            return((cutset, Fe))
    return()
    
    
# Leer datos
instancia = 'Codigos/Data/tsp/att48.tsp'
problem = tsplib95.load(instancia)
# Datos del TSP
n = len(list(problem.get_nodes()))
inicio = list(problem.get_nodes())[0]
ciudades = [i for i in range(n)]
arcos = dict()
# Grafo usando la biblioteca networkx 
G = nx.complete_graph(n)

info = problem.as_keyword_dict()
my_pos = {}
for i in range(1, n + 1):
    x, y = info['NODE_COORD_SECTION'][i]
    my_pos[i-1] = x, y

#nx.draw(G, pos=my_pos, node_size=15, width=0.1)
#plt.show()

inicio = list(problem.get_nodes())[0]
if inicio == 0:
    for i,j in G.edges:
        G.edges[i,j]['length'] = problem.get_weight(*(i, j))
else:
    for i,j in G.edges:
        G.edges[i,j]['length'] = problem.get_weight(i+1, j+1)

modelo = gp.Model()
#x = tupledict()
#for i, j in G.edges:
#    x[i,j] = modelo.addVar(obj = G.edges[i,j]['length'], vtype = GRB.BINARY, name = "x[%d,%d]" % (i,j))
                               
#x = modelo.addVars(G.edges, vtype=GRB.BINARY, name = "x")
x = modelo.addVars(G.edges, vtype=GRB.CONTINUOUS, lb = 0, ub = 1, name = "x")
modelo.setObjective(gp.quicksum(G.edges[e]['length'] * x[e] for e in G.edges), GRB.MINIMIZE)

# usando 2-matching
modelo.addConstrs(gp.quicksum(x[e] for e in G.edges if e in G.edges(i)) == 2 for i in G.nodes)

print(modelo.display())

modelo.optimize()      
print("Iteración 0")
print("Costo  : %g" % modelo.ObjVal)
print("Tiempo : %f" % modelo.Runtime)

# Solución relajada
ruta = [e for e in G.edges if x[e].x > 0.999]
#nx.draw(G.edge_subgraph(ruta), pos=my_pos, node_size=15, width=0.5)
#plt.show()
contar = 0
# Se debe agregar las restricciones y volver a resolver iteradamente
tiempo = 0.0  
while not nx.is_connected(G.edge_subgraph(ruta)):
    contar += 1
    
    #xx = [modelo.getVars()[i].X for i in G.edges]
    #xx = [x[e].x for e in G.edges]
    #print(xx)
    
    
    #W,F = Padberg_Rao(G, xx)
    #print("Conjunto Blossom", W, F)
    #modelo.addConstr(suma(W.difference(F)) - suma(F) >= 1 - len(F))
    #print(suma(W. difference (F)) - suma(F) >= 1 - len(F))
    
    #exit(0)
    
    
    for component in nx.connected_components(G.edge_subgraph(ruta)):
        print("Agregando restricciones para esta subruta", component)
        aristas = [(i,j) for (i,j) in G.edges if i in component and j in component]
        modelo.addConstr(gp.quicksum(x[e] for e in aristas) <= len(component) - 1)
    
    modelo.optimize()
    ruta = [e for e in G.edges if x[e].x > 0.9]
    #nx.draw(G.edge_subgraph(ruta), pos=my_pos, node_size=15, width=0.5)
    #plt.show()
    tiempo += modelo.Runtime
    print("Iteración", contar)
    print("Costo  : %g" % modelo.ObjVal)
    print("Tiempo : %f" % modelo.Runtime)

print("Tiempo Final: %f" % tiempo)
print("Costo  Final: %g" % modelo.ObjVal)

