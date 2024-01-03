import gurobipy as gp
from gurobipy import GRB
import sys
import warnings
import networkx as nx
import time
warnings.filterwarnings("ignore")

from Source.MathematicalModel import MathematicalModel as MILP

def integer_fractional_cut1(modelo:gp.Model, donde):
    initial = time.time()
    n = modelo._n
    if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
        valoresX = modelo.cbGetNodeRel(modelo._xvars)
        tour = [i for i in range(n+1)]
        MILP.subtour_method(tour, valoresX,n)
        if len(tour) < n:
            modelo._callback_count +=1 
            tour2 = [i for i in range(n) if i not in tour]
            modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)
            modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour if i != j) <= len(tour)-1)    
    modelo._callback_time += time.time()-initial

def integer_fractional_cut2(modelo:gp.Model, donde):
    initial = time.time()
    n = modelo._n
    if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
        valoresX = modelo.cbGetNodeRel(modelo._xvars)
        tour = [i for i in range(n+1)]
        MILP.subtour_method(tour, valoresX,n)
        if len(tour) < n:
            modelo._callback_count +=1 
            tour2 = [i for i in range(n) if i not in tour]
            modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)
            #modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour if i != j) <= len(tour)-1)    

        
        DG: nx.DiGraph = modelo._DG
        for i,j in DG.edges:
            DG.edges[i,j]['capacity'] = valoresX[i,j]

        for t in range(1,n):
            (cut_value, node_partition) = nx.minimum_cut( DG, _s=0, _t=t )
            # print("cut_value =",cut_value)
            if cut_value < 1 - modelo._epsilon:
                modelo._callback_count +=1 
                S = node_partition[0]  # 'left' side of the cut
                T = node_partition[1]  # 'right' side of the cut
                modelo.cbLazy( gp.quicksum( modelo._xvars[i,j] for i in S for j in T ) >= 1 )
                return

    modelo._callback_time += time.time()-initial

def integer_fractional_cut3(modelo:gp.Model, donde):
    initial = time.time()
    n = modelo._n
    if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
        valoresX = modelo.cbGetNodeRel(modelo._xvars)
        tour = [i for i in range(n+1)]
        MILP.subtour_method(tour, valoresX,n)
        if len(tour) < n:
            modelo._callback_count +=1 
            tour2 = [i for i in range(n) if i not in tour]
            modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)
            modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour if i != j) <= len(tour)-1)    

        DG: nx.DiGraph = modelo._DG
        for i,j in DG.edges:
            DG.edges[i,j]['capacity'] = valoresX[i,j]

        for t in range(1,n):
            (cut_value, node_partition) = nx.minimum_cut( DG, _s=0, _t=t )
            # print("cut_value =",cut_value)
            if cut_value < 1 - modelo._epsilon:
                modelo._callback_count +=1 
                S = node_partition[0]  # 'left' side of the cut
                T = node_partition[1]  # 'right' side of the cut
                modelo.cbLazy( gp.quicksum( modelo._xvars[i,j] for i in S for j in T ) >= 1 )
                return

    modelo._callback_time += time.time()-initial

argv = sys.argv[1:]
opts = [(argv[2*i],argv[2*i+1]) for i in range(int(len(argv)/2))]

size = "tsplib"
instance = "gr17"
output = True
subtour = "wc"
initial_sol = True
callback =  'subtourelim1'
bounds = True
new_formulation = True
time_limit = 1800
new_m = False

for i in range(len(opts)):
    if opts[i][0][1:] == "size":  size  = opts[i][1]
    elif opts[i][0][1:] == "instance": 
        try:
            instance = int(opts[i][1])  
        except ValueError:
            instance = opts[i][1]  
    elif   opts[i][0][1:] == "initial_sol": initial_sol = False if opts[i][1].lower() == "false" else True
    elif   opts[i][0][1:] == "output": output =  False if opts[i][1].lower() == "false" else True
    elif   opts[i][0][1:] == "bounds": bounds = False if opts[i][1].lower() == "false" else True
    elif   opts[i][0][1:] == "subtour": subtour = opts[i][1]
    elif   opts[i][0][1:] == "callback": 
        if opts[i][1] == "custom1":
            callback = integer_fractional_cut1
        elif opts[i][1] == "custom2":
            callback = integer_fractional_cut2
        elif opts[i][1] == "custom3":
            callback = integer_fractional_cut3
        else:
            callback = opts[i][1]
    elif   opts[i][0][1:] == "newformulation": new_formulation =  False if opts[i][1].lower() == "false" else True
    elif   opts[i][0][1:] == "newm": new_m =  False if opts[i][1].lower() == "false" else True

MM = MILP(size,
          instance,
          output,
          subtour,
          initial_sol,
          callback,
          bounds,
          new_formulation,
          time_limit,
          new_m = new_m
          )
MM.run()
MM.print_results()
