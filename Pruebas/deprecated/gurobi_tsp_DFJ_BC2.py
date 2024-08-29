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
    cap = [G[e[0]][e[1]][ capacity] for e in F]
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
    
                              
# Callback - para usar cortes lazy de eliminación de subtour Eq DFJ
def subtour_elim_DFJ(modelo, donde):
    if donde == GRB.Callback.MIPNODE  and modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        M = nx.Graph()
        xval = modelo.cbGetNodeRel(modelo._x)
        for i, e in enumerate(modelo._G.edges):
            M.add_edge(e[0],e[1], weight = xval[e])
            #print(i, e, xval[e])
        #print(M)
        # Fijamos un parametro epsilon de tolerancia para la segunda heuristica. Tal y como Groetschel y Holland en su trabajo , tomamos 0.3.
        epsilon = 0.3
        # Obtenemos las aristas con variables asociadas en
        aristasN = [e for e in M.edges () if totalcap(M, [e]) > 0]
        aristasE = [e for e in M.edges () if totalcap(M, [e]) > epsilon]
        # Construimos nuevos grafos con estas aristas.
        GN = nx.Graph ()
        GN.add_edges_from(aristasN)
        
        GE = nx.Graph ()
        GE.add_edges_from(aristasE)

        conN = [a for a in list(nx.connected_components (GN)) if len(a) % 2 == 1]
        # Tenemos que tener en cuenta que los planos de corte que proporciona la segunda heuristica pueden no ser utiles , ya que nuestra solucion no tiene por que
        # violarlos necesariamente y seria un esfuerzo innecesario volver a reoptimizar para nada, por lo que solo nos quedamos con dichas componentes.
        # Definimos previamente una funcion que nos haga dicha comprobacion
        def comprueba (nodos):
            E = GE.subgraph(nodos).edges()
            suma = totalcap (M,E)
            card = (len(nodos) - 1) / 2
            return(suma <= card)
        
        # Generamos las componentes adecuadas
        conE = [a for a in list(nx.connected_components(GE)) if len(a) % 2 == 1 and comprueba (a)]
        # Si las lista con es no vacia , hemos encontrado planos de cortes , los añadimos a nuestro modelo y volvemos a comprobar.
        if conN != []:
            return
        elif conE != []:
            return
        # Paso 7.3: En el caso de que las heuristicas no hayan funcionado , construimos un plano de corte mediante el procedimiento de Padberg y Rao.
        else:           
            xx = [xval[e] for e in modelo._G.edges]
            W,F = Padberg_Rao(modelo._G, xx)
            print("Conjunto Blossom", W)
            ruta = [e for e in W]
            #
            #print(modelo._G.edges.edge_subgraph(ruta))
            #print(xx)
            #ruta = [e for e in modelo._G.edges if xval[e] <0.9 and xval[e] >0]
            print(ruta)
            print(len(modelo._G.edges))
            print(modelo._G.nodes())
            #nx.draw(modelo._G.edges.edge_subgraph(ruta), pos=my_pos, node_size=15, width=0.5)
            nx.draw(modelo._G.edge_subgraph(W), pos=my_pos, node_size=15, width=0.5)
            plt.show()
            exit(0)
            modelo.cbCut((suma(W.difference(F)) - suma(F) >= 1 - len(F)))
            print("Hemos usado Padberg-Rao")
            print(suma(W. difference (F)) - suma(F) >= 1 - len(F))
            

    if donde == GRB.Callback.MIPSOL:
        xval = modelo.cbGetSolution(modelo._x)
        ruta = [e for e in modelo._G.edges if xval[e] > 0.5]
        
        for componente in nx.connected_components(modelo._G.edge_subgraph(ruta)):          
            if len(componente) <= modelo._G.number_of_nodes() / 2:  
                modelo._contarCut += 1     
                # agregar cortes de elimination de subtour DFJ1    
                aristas = [(i,j) for (i,j) in modelo._G.edges if i in componente and j in componente]
                modelo.cbLazy(gp.quicksum(modelo._x[e] for e in aristas) <= len(componente) - 1 )
             
# Leer datos
instancia = 'Codigos/Data/tsp/att48.tsp'#sys.argv[1]

problem = tsplib95.load(instancia)
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

nx.draw(G, pos=my_pos, node_size=15, width=0.1)
plt.show()

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
                               
x = modelo.addVars(G.edges, vtype=GRB.BINARY, name = "x")
modelo.setObjective(gp.quicksum(G.edges[e]['length'] * x[e] for e in G.edges), GRB.MINIMIZE)

# usando 2-matching
modelo.addConstrs(gp.quicksum(x[e] for e in G.edges if e in G.edges(i)) == 2 for i in G.nodes)

print(modelo.display())

# Parámetros
modelo.Params.Threads = 1
modelo.Params.LazyConstraints = 1
modelo._x = x
modelo._G = G
modelo._contarCut = 0
# imprimir modelo
# print(modelo.display())

modelo.optimize(subtour_elim_DFJ)

print("Costo  : %g" % modelo.ObjVal)
print("Tiempo : %f" % modelo.Runtime)
print("Cortes : %g" % modelo._contarCut)
ruta = [e for e in G.edges if x[e].x > 0.9]
print(ruta)
nx.draw(G.edge_subgraph(ruta), pos=my_pos, node_size=15, width=0.5)
plt.show()
