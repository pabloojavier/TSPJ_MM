import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import numpy.typing as npt
import warnings
import networkx as nx
from typing import Dict, List, Callable, Any,Tuple
from Source.Problem import Problem 
import time
import re
import math
import matplotlib.pyplot as plt
import lib.sec as sc

warnings.filterwarnings("ignore")

import numpy as np

def get_min_job(JT: np.ndarray, ja):
    """
    Returns the minimum jobtime from a relaxed job assignments.

    Parameters:
    JT (np.ndarray): The input array.
    ja (list): A list of tuples representing job assignments and his relaxed value.

    Returns:
    float: The minimum value from the modified array.
    """
    JT = JT.T
    df = JT.copy()
    assigns = [i[0] for i in ja if i[1] > 0.99]
    for i in assigns:
        df[:, i[0]] = 99999
        df[i[1], i[0]] = JT[i[1]][i[0]]

    for i in assigns:
        df[i[1]] = 99999
        df[i[1], i[0]] = JT[i[1]][i[0]]

    return df[df > 0].min()

def NNJA(route,JT):
    """
    Neares Neighbor Algorithm for Job Assignment. From the last node of the tour assigns
    the cheapest job available. Then in the next node assing the cheapest job availabe 
    and so on until the first node
    """
    job = []
    n = len(route)  
    cont = len(route)-1
    while len(job)<n:
        times = {(i,route[cont]):JT[route[cont]][i] for i in range(1,n+1)}
        new = min(times.items(),key=lambda x:x[1]) 
        job.append(new[0][1])
        cont -= 1
    return job

def sort_jobs(route,jobs):
    return [(route[i],jobs[i]) for i in range(len(route))]

def sort_arch(route):
    """
    Aux method, from a list solution to a arch solution.
    """
    archs = []
    for i in range(len(route)-1):
        archs.append((str(route[i]),str(route[i+1])))
    archs.append((str(route[-1]),str(route[0])))
    return archs

def delta(G,U): 
    ejes = set()
    for l, nbrs in ((n, G[n]) for n in U):
        ejes.update( (l,y) if l<y else (y,l) for y in nbrs if y not in U)
    return ejes

def term(e,H):
    s = max(H[e[0]][e[1]]['weight'], 1-H[e[0]][e[1]]['weight'])
    t = min(H[e[0]][e[1]]['weight'],1-H[e[0]][e[1]]['weight']) 
    return(s-t)

def totalcap(G,F,capacity = 'weight'):
    cap = [G[e[0]][e[1]][capacity] for e in F]
    val = sum(cap) 
    return(val)

def letchford_algorithm(
        G:nx.Graph,
        x_values:Dict[Tuple[int,int],float],
        symmetric:bool = False
        ) -> tuple:
    """
    Letchford algorithm for the TSP problem.

    Letchford, Adam; Reinelt, Gerhard; Theis, Dirk Oliver. (2004). A Faster Exact Separation Algorithm for Blossom Inequalities. 10.1007/978-3-540-25960-215.
    """
    if symmetric:
        c = {key:min(value+x_values[(key[1],key[0])],1-value-x_values[(key[1],key[0])]) for key,value in x_values.items()}
    else:
        c = {key:min(value,1-value) for key,value in x_values.items()}

    H = nx.Graph()
    for e in G.edges:
        if x_values[e] > 0 and e[0]<e[1]:
            if symmetric:
                H.add_edge(e[0],e[1],capacity = c[e], weight = x_values[e])
            else:
                H.add_edge(e[0],e[1],capacity = c[e], weight = x_values[e]+x_values[(e[1],e[0])])

    GomHu:nx.Graph = nx.gomory_hu_tree(H, capacity = 'capacity')
    
    for e in GomHu.edges:
        Te = nx.Graph(GomHu)
        Te.remove_edge(*e)
        U,V = list(nx.connected_components(Te))
        
        if len(U)<len(V):
            U_,V_ = U,V
        else:
            V_,U_ = U,V
        del V,V_
        
        if len(U_)<3:
            continue

        delta_W = delta(H,U_)
        Fe = []
        added_nodes = []
        for e in delta(H,U_):
            if 1 - H[e[0]][e[1]]['weight'] < H[e[0]][e[1]]['weight'] :
                if e[0] in added_nodes or e[1] in added_nodes:
                    continue
                Fe.append(e)
                added_nodes.append(e[0])
                added_nodes.append(e[1])
        del added_nodes
        
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
            return U_,Fe

    return None,None

def heuristic_separation(valoresX,n):
    delta_value = 0.0000001
    aux_valoresX = {key:value for key,value in valoresX.items() if value >0 and value <1-delta_value}
    fractional_G = nx.Graph()
    fractional_G.add_nodes_from([i for i in range(n)])

    complete_G = nx.Graph()
    complete_G.add_nodes_from([i for i in range(n)])

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
        return fixed_component,fixed_T
    return [],[]


class MathematicalModel(Problem):
    def __init__(self,size:str,
                 instance,
                 output = False,
                 subtour : str = "GG",
                 initial_solution : bool = True,
                 callback : str = "None",
                 bounds: bool = True,
                 new_formulation: bool = False,
                 time_limit : int = 7200,
                 new_m : bool = False,
                 relax : bool = False,
                 mga : bool = False,
                 exp_lb : int = 0
                 ):
        """
        Initialize mathematical model with some default values

        Default values:

            output = True\n
            subtour = GG\n
            initial_solution = True\n
        """

        super().__init__(size,instance)
        self.size = size
        self.output = output
        self.subtour = subtour.lower()
        self.initial_solution = initial_solution
        self.callback = callback.lower() if isinstance(callback,str) else "custom"
        self.bounds = bounds
        self.time_limit = time_limit
        self.new_m = new_m
        self.relax = relax
        self.new_formulation = new_formulation
        self.mga = mga
        if int(exp_lb) in [0,1,2]:
            self.exp_lb = int(exp_lb)
        else:
            raise ValueError("exp_lb must be 0, 1 or 2")

        self.callback_dict = {"cut_integer_separation":MathematicalModel.CUT_integer_separation,
                              "cut_naive_fractional_separation":MathematicalModel.CUT_naive_fractional_separation,
                              "cut_smarter_fractional_separation":MathematicalModel.CUT_smarter_fractional_separation,
                              "dfj_integer_separation":MathematicalModel.DFJ_integer_separation,
                              "dfj_naive_fractional_separation":MathematicalModel.DFJ_naive_fractional_separation,
                              "dfj_smarter_fractional_separation":MathematicalModel.DFJ_smarter_fractional_separation,
                              "subtourelim1":MathematicalModel.subtourelim1,
                              "subtourelim2":MathematicalModel.subtourelim2,
                              "subtourelim3":MathematicalModel.subtourelim3,
                              "subtourelim4":MathematicalModel.subtourelim4,
                              "subtourelim5":MathematicalModel.subtourelim5,
                              "subtourelim6":MathematicalModel.subtourelim6,  
                              "subtourelim7":MathematicalModel.subtourelim7,  
                              "exact_blossom":MathematicalModel.exact_blossom_method,
                              "heuristic_blossom":MathematicalModel.heuristic_blossom_method,
                              "both_blossom":MathematicalModel.both_blossom_method,
                              "exact_blossom2":MathematicalModel.exact_blossom_method2,
                              "heuristic_blossom2":MathematicalModel.heuristic_blossom_method2,
                              "both_blossom2":MathematicalModel.both_blossom_method2,
                              "both_blossom3":MathematicalModel.both_blossom_method3,
                              'custom':callback}
        if self.output:
            print(f'running with: {self.size} {self.instance} {self.subtour} {self.initial_solution} {self.callback} {self.bounds} {self.new_formulation} {self.time_limit} {self.new_m} {self.relax} {self.mga}')
        self.jobs = self.cities.copy()
        self.jobs_arch = [(i,k) for i in self.cities for k in self.cities]
        self.initial_solution_time = time.time()

        self.create_initial_solution()
        self.initial_solution_time = time.time() - self.initial_solution_time

        self.compute_new_M() if new_m else self.compute_M()
        self.jt_min = self.JT[self.JT>0].min()

    def compute_M(self):
        """
        This method computes the M value, it can be overwritting.
        """
        self.M = 0
        self.sum_min_row = 0
        for i in range(self.n):
            self.sum_min_row += np.min(self.TT[i][np.nonzero(self.TT[i])])
            max_t = 0
            max_tt = 0
            for j in range(self.n):
                if i!=j and self.TT[i][j]>max_t:
                    max_t = self.TT[i][j]
                if j!= 0 and self.JT[i][j]>max_tt:
                    max_tt = self.JT[i][j]
            self.M += max_t #+max_tt
        self.M = np.full((len(self.cities), len(self.cities)), self.M)

    def compute_new_M(self):
        menor_arco_depot = min(self.TT[0][i] for i in range(1,self.n))

        self.M = np.zeros((len(self.cities),len(self.cities)))
        self.sum_min_row = 0
        for j in range(1,len(self.cities)):
            self.M[0,j] = self.TT[0][j]
        for i in range(1, len(self.cities)):
            self.sum_min_row += np.min(self.TT[i][np.nonzero(self.TT[i])])
            for j in range(len(self.cities)):
                if i != j:
                    if self.exp_lb == 0:
                        self.M[i,j] = self.initial_route_fitness + self.TT[i][j] 
                    elif self.exp_lb == 1:
                        self.M[i,j] = self.initial_route_fitness + self.TT[i][j] - self.TT[self.last_node_initial_route][0]
                    elif self.exp_lb == 2:
                        self.M[i,j] = self.initial_route_fitness + self.TT[i][j] - menor_arco_depot

    def create_initial_solution(self):
        self.initial_jobs = NNJA(self.initial_route[1:],self.JT) 
        if self.mga:
            from Source.MGA import MGA
            mga = MGA(self.size,self.instance,seed = 0)
            mga.parameters['TIMELIMIT'] = 100
            mga.run()
            mga_solution = mga.get_solution()
            self.initial_route = [0]+mga_solution[0]
            self.initial_jobs = mga_solution[1]

        self.initial_fitness = self.fitness_functions([self.initial_route[1:],self.initial_jobs])[0]
        self.initial_route_fitness = self.route_fitness(self.initial_route[1:])
        self.initial_arch = sort_arch(self.initial_route)
        self.initial_job_arch = sort_jobs(self.initial_route[1:],self.initial_jobs)
        self.last_node_initial_route = self.initial_route[-1]

    def create_base_model(self):
        """
        Create the initial MILP, whithout any additional constraint.
        """
        env = gp.Env(empty=True)
        env.setParam('OutputFlag', 1 if self.output else 0)
        env.start()
        self.modelo = gp.Model(env=env)
        self.modelo._n_subtour1_constraints = 0    #Number of classic connectivity constraints added
        self.modelo._n_subtour2_constraints = 0    
        self.modelo._n_blossom_heuristic_constraints = 0 #Number of heuristic constraints added
        self.modelo._n_blossom_exact_constraints = 0     #Number of exact constraints added

        self.modelo._time_subtour1_constraints = 0           #Time for connectivity constraints
        self.modelo._time_subtour2_constraints = 0
        self.modelo._time_blossom_heuristic_constraints = 0 #Time for heuristic constraints
        self.modelo._time_blossom_exact_constraints = 0     #Time for exact constraints
        self.modelo._time_total_constraints = 0             #Total time for constraints

        self.modelo._simplex_lower_bound = 1e-9
        self.modelo._relax = self.relax
        self.modelo._stop_flag = False

        self.Cmax = self.modelo.addVar(vtype=GRB.CONTINUOUS,name="Cmax")

        var_mode = GRB.BINARY
        self.x = self.modelo.addVars(self.arch, vtype=var_mode, name='x')
        self.y = self.modelo.addVars(self.jobs_arch, vtype=var_mode, name='y')
        self.TS = self.modelo.addVars(self.cities,vtype=GRB.CONTINUOUS,name="TS")

        self.modelo.setObjective(self.Cmax, GRB.MINIMIZE)

        #Restricción adicional de reforzamiento
        self.modelo.addConstr(self.Cmax >= self.jt_min + gp.quicksum(self.x[(i,j)]*self.TT[i][j] for i in self.cities for j in self.cities[1:] if i!=j),name="Reinforcement")

        self.modelo.addConstrs((self.x.sum(i,'*') == 1 for i in self.cities) , name = 'Outgoing') # Outgoing
        self.modelo.addConstrs((self.x.sum('*', j) == 1 for j in self.cities), name = 'Incoming') # Incoming 

        for k in self.jobs[1:self.n]:
            self.modelo.addConstr(gp.quicksum(self.y[(i,k)] for i in self.cities if i != 0) == 1 , name = f'Job_{k}_out')

        for i in self.cities[1:self.n]: 
            self.modelo.addConstr(gp.quicksum(self.y[(i,k)] for k in self.jobs if k != 0) == 1 , name = f'Job_{i}_in')
        
        for i in self.cities[1:self.n]:
            self.modelo.addConstr(self.Cmax >= self.TS[i] + gp.quicksum(self.y[(i,k)]*self.JT[i][k] for k in self.jobs if k!=0 ) , name = f'Cmax_{i}')

        for i in self.cities[1:self.n]:
            self.modelo.addConstr(self.Cmax >= self.TS[i] + self.x[(i,0)]*self.TT[i][0] , name = f'Cmax_{i}_0')
        
        for i in self.cities: #16
            for j in self.cities[1:self.n]:
                if i!=j:
                    self.modelo.addConstr(self.TS[i] + self.TT[i][j] - (1-self.x[(i,j)])*self.M[i][j] <= self.TS[j] , name = f'TS_{i}_{j}')

        
        self.modelo.Params.Threads = 1
        self.modelo._callback_time = 0
        self.modelo.update()

    def create_new_formulation(self):
        """
        Create the initial new formulation MILP, whithout any additional constraint.
        """
        env = gp.Env(empty=True)
        env.setParam('OutputFlag', 1 if self.output else 0)
        env.start()
        self.modelo = gp.Model(env=env)

        self.modelo._n_subtour1_constraints = 0    #Number of classic connectivity constraints added
        self.modelo._n_subtour2_constraints = 0    
        self.modelo._n_blossom_heuristic_constraints = 0 #Number of heuristic constraints added
        self.modelo._n_blossom_exact_constraints = 0     #Number of exact constraints added

        self.modelo._time_subtour1_constraints = 0           #Time for connectivity constraints
        self.modelo._time_subtour2_constraints = 0
        self.modelo._time_blossom_heuristic_constraints = 0 #Time for heuristic constraints
        self.modelo._time_blossom_exact_constraints = 0     #Time for exact constraints
        self.modelo._time_total_constraints = 0             #Total time for constraints

        self.modelo._simplex_lower_bound = 1e-9
        self.modelo._relax = self.relax
        self.modelo._stop_flag = False

        self.Cmax = self.modelo.addVar(vtype=GRB.CONTINUOUS,name="Cmax")

        var_mode = GRB.BINARY
        self.x = self.modelo.addVars(self.arch, vtype=var_mode, name='x')
        self.y = self.modelo.addVars(self.jobs_arch, vtype=var_mode, name='y')

        self.t = self.modelo.addVars(self.arch,vtype=GRB.CONTINUOUS,name="t", lb = 0) # new variable for the time

        self.modelo.setObjective(self.Cmax, GRB.MINIMIZE)

        #Restricción adicional de reforzamiento
        self.modelo.addConstr(self.Cmax >= self.jt_min + gp.quicksum(self.x[(i,j)]*self.TT[i][j] for i in self.cities for j in self.cities[1:] if i!=j), name = "Reinforcement")
        
        for i in self.cities[1:self.n]:
            self.modelo.addConstr(self.Cmax >= gp.quicksum(self.t[(i,k)] for k in self.cities if i != k) 
                                             + gp.quicksum(self.y[(i,k)]*self.JT[i][k] for k in self.jobs if k!=0 ) , name = f'Cmax_{i}')

        for i in self.cities[1:self.n]:
            self.modelo.addConstr(self.Cmax >= gp.quicksum(self.t[(i,k)] for k in self.cities if i != k) 
                                             + self.x[(i,0)]*self.TT[0][i] , name = f'Cmax_{i}_0')

        for k in self.jobs[1:self.n]:
            self.modelo.addConstr(gp.quicksum(self.y[(i,k)] for i in self.cities if i != 0) == 1 , name = f'Job_{k}_out')

        for i in self.cities[1:self.n]: 
            self.modelo.addConstr(gp.quicksum(self.y[(i,k)] for k in self.jobs if k != 0) == 1 , name = f'Job_{i}_in')
        
        self.modelo.addConstrs((self.x.sum(i,'*') == 1 for i in self.cities) , name = 'Outgoing') # Outgoing
        self.modelo.addConstrs((self.x.sum('*', j) == 1 for j in self.cities) , name = 'Incoming') # Incoming 
        

        for k in self.cities[1:self.n]: #16
            self.modelo.addConstr(gp.quicksum(self.t[(i,k)] for i in self.cities if i != k) 
                                + gp.quicksum(self.TT[i][k]*self.x[(i,k)] for i in self.cities if i != k)
                                <= gp.quicksum(self.t[(k,l)] for l in self.cities if k != l) , name = f't_{k}')

        #parte desde 1, el primer nodo no tiene tiempo
        #for i in self.cities[1:]:
        #    LB = np.min(self.TT[i][np.nonzero(self.TT[i])])
        #    for k in self.cities:
        #        if i != k:
        #            self.modelo.addConstr(self.t[(i,k)] <= self.M[i][k]*self.x[(i,k)] , name = f't_{i}_{k}_UB')
        #            self.modelo.addConstr(self.t[(i,k)] >= 0 , name = f't_{i}_{k}_LB')
        
        self.modelo.Params.Threads = 1
        #0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent
        self.modelo.setParam("Method",1) 
        
        self.modelo.update()

    def add_subtour_constraint(self):
        """
        Adds subtour constraints, it can be GG, MTZ or DL constraint
        """
        if self.subtour == "gg":
            self.g = self.modelo.addVars(self.arch,name = "y") 
            for i in self.cities[1:self.n]: #14
                self.modelo.addConstr(gp.quicksum(self.g[(i,j)] for j in self.cities if i !=j) - gp.quicksum(self.g[(j,i)] for j in self.cities if i !=j) == 1 , name = f'GG1_{i}')

            for i in self.cities[1:self.n]: #15
                for j in self.cities:
                    if i!=j:
                        self.modelo.addConstr(self.g[(i,j)]<= self.n*self.x[(i,j)] , name=f'GG2_{i}_{j}')
            
        elif self.subtour in ("mtz","dl","dl_real"):
            self.u = self.modelo.addVars(self.cities , vtype = GRB.CONTINUOUS , name = "u")
            if self.subtour == "mtz":    
                for i,j in self.arch:
                    if i>0:
                        #self.M debería ser self.n
                        self.modelo.addConstr(self.u[i] - self.u[j] + 1 <= self.n * (1 - self.x[(i,j)]) , f"MTZ({i},{j})")
            
            elif self.subtour == "dl":
                for i in range(1,self.n):
                    self.modelo.addConstr(self.u[i] >= 1 + (self.n-3)*self.x[(i,0)] + gp.quicksum(self.x[(j,i)] for j in self.cities[1:] if j != i), name=f"DL_{i}_LB")

                for i in range(1,self.n):
                    self.modelo.addConstr(self.u[i] <= self.n-1 - (self.n-3)*self.x[(0,i)]- gp.quicksum(self.x[(i,j)] for j in self.cities[1:] if j != i), name = f"DL_{i}_UB")
            
            elif self.subtour == "dl_real":
                for i,j in self.arch:
                    if i>0:
                        self.modelo.addConstr(self.u[i] - self.u[j] + (self.n-1)*self.x[(i,j)] + (self.n-3)*self.x[(j,i)] <= self.n - 2 , name = f"DL_real_{i}_{j})")

                for i in range(1,self.n):
                    self.modelo.addConstr(self.u[i] >= 1 + (self.n-3)*self.x[(i,0)] + gp.quicksum(self.x[(j,i)] for j in self.cities[1:] if j != i), name = f"DL_real_{i}_LB")

                for i in range(1,self.n):
                    self.modelo.addConstr(self.u[i] <= self.n-1 - (self.n-3)*self.x[(0,i)]- gp.quicksum(self.x[(i,j)] for j in self.cities[1:] if j != i) , name = f"DL_real_{i}_UB")
        elif self.subtour != "wc":
            print(self.subtour)
            raise ValueError("Subtour method not implemented <wc/gg/dl/mtz/wc>")
        self.modelo.update()

    @staticmethod
    def integer_fractional_cut(modelo:gp.Model, donde):
        initial = time.time()
        n = modelo._n

        if donde == gp.GRB.Callback.SIMPLEX:
            modelo._simplex_lower_bound = modelo.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)
            if len(tour) < n:
                # modelo._callback_count +=1 
                modelo._n_subtour1_constraints += 1
                tour2 = [i for i in range(n) if i not in tour]
                modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)
                modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour if i != j) <= len(tour)-1)
            

            # DG: nx.DiGraph = modelo._DG
            # for i,j in DG.edges:
            #     DG.edges[i,j]['capacity'] = valoresX[i,j]

            # for t in range(1,n):
            #     (cut_value, node_partition) = nx.minimum_cut( DG, _s=0, _t=t )
            #     # print("cut_value =",cut_value)
            #     if cut_value < 1 - modelo._epsilon:
            #         modelo._callback_count +=1 
            #         S = node_partition[0]  # 'left' side of the cut
            #         T = node_partition[1]  # 'right' side of the cut
            #         modelo.cbLazy( gp.quicksum( modelo._xvars[i,j] for i in S for j in T ) >= 1 )
            #         return

        # modelo._callback_time += time.time()-initial
        modelo._time_subtour1_constraints += time.time()-initial
        modelo._time_total_constraints += time.time()-initial
    
    @staticmethod
    def CUT_integer_separation(m:gp.Model, where):
        initial = time.time()
        # check if LP relaxation at this branch-and-bound node has an integer solution
        
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if where == GRB.Callback.MIPSOL: 
            
            # retrieve the LP solution
            xval = m.cbGetSolution(m._xvars)
            
            # which edges are selected?
            edges_used = [ (i,j) for i,j in m._DG.edges if xval[i,j] > 0.5 ]
            
            # create support graph
            DG_soln = m._DG.edge_subgraph( edges_used )
            
            # if solution is not connected, add a (violated) CUT constraint for each subtour
            if not nx.is_strongly_connected( DG_soln ):
                for component in nx.strongly_connected_components( DG_soln ):
                    # m._callback_count +=1 
                    m._n_subtour1_constraints += 1
                    complement = [ i for i in DG_soln.nodes if i not in component ]
                    m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in complement ) >= 1 )
        # m._callback_time += time.time()-initial
        m._time_subtour1_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def CUT_naive_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case1 and not case2:
            m._time_subtour1_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        
        # retrieve the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
            
        # which edges are selected (in whole or in part)?
        DG = m._DG
        edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]

        # check if any CUT inequalities are violated (by solving some min cut problems)
        for i,j in DG.edges:
            DG.edges[i,j]['capacity'] = xval[i,j]

        s = 0
        for t in range(1,m._n):
            (cut_value, node_partition) = nx.minimum_cut( DG, _s=s, _t=t )
            # print("cut_value =",cut_value)
            if cut_value < 1 - m._epsilon:
                # m._callback_count +=1 
                m._n_subtour1_constraints += 1
                S = node_partition[0]  # 'left' side of the cut
                T = node_partition[1]  # 'right' side of the cut
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in T ) >= 1 )
                # m._callback_time += time.time()-initial
                m._time_subtour1_constraints += time.time()-initial
                m._time_total_constraints += time.time()-initial
                return
        m._time_subtour1_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def CUT_smarter_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #

        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        n = m._n
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if not case1 and not case2:
            # m._callback_time += time.time()-initial
            m._time_subtour1_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        
        # retrieve the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
        
        DG = m._DG
        # if the support graph is disconnected, then finding violated cuts is easy!
        edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]
        DG_support = DG.edge_subgraph( edges_used )
        if not nx.is_strongly_connected( DG_support ):
            for component in nx.strongly_connected_components( DG_support ):
                # m._callback_count +=1 
                m._n_subtour1_constraints += 1
                complement = [ i for i in DG.nodes if i not in component ]
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in component for j in complement ) >= 1 )
        else:
            # check if any CUT inequalities are violated (by solving some min cut problems)
            for i,j in DG.edges:
                DG.edges[i,j]['capacity'] = xval[i,j]

            s = 0
            for t in range(1,n):
                (cut_value, node_partition) = nx.minimum_cut( DG, _s=s, _t=t )
                # print("cut_value =",cut_value)
                if cut_value < 1 - m._epsilon:
                    m._callback_count +=1 
                    S = node_partition[0]  # 'left' side of the cut
                    T = node_partition[1]  # 'right' side of the cut
                    m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in T ) >= 1 )
                    # m._callback_time += time.time()-initial
                    m._time_subtour1_constraints += time.time()-initial
                    m._time_total_constraints += time.time()-initial
                    return
        m._time_subtour1_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def DFJ_integer_separation(m:gp.Model, where):
        initial = time.time()
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        
        # check if LP relaxation at this branch-and-bound node has an integer solution
        if where == GRB.Callback.MIPSOL: 
            
            # retrieve the LP solution
            xval = m.cbGetSolution(m._xvars)
            
            # which edges are selected?

            edges_used = [ (i,j) for i,j in m._DG.edges if xval[i,j] > 0.5 ]
            
            # create support graph
            DG_soln = m._DG.edge_subgraph( edges_used )
            
            # if solution is not connected, add a (violated) DFJ constraint for each subtour
            if not nx.is_strongly_connected( DG_soln ):
                for component in nx.strongly_connected_components( DG_soln ):
                    # m._callback_count +=1 
                    m._n_subtour1_constraints += 1
                    m.cbLazy( gp.quicksum( m._xvars[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )
        # m._callback_time += time.time()-initial
        m._time_subtour1_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def DFJ_naive_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case1 and not case2:
            # m._callback_time += time.time()-initial
            m._time_subtour1_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        
        # retrieve the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
            
        # which edges are selected (in whole or in part)?
        DG = m._DG
        edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]

        # check if any CUT inequalities are violated (by solving some min cut problems)
        for i,j in DG.edges:
            DG.edges[i,j]['capacity'] = xval[i,j]

        s = 0
        for t in range(1,m._n):
            (cut_value, node_partition) = nx.minimum_cut( DG, _s=s, _t=t )
            # print("cut_value =",cut_value)
            if cut_value < 1 - m._epsilon:
                # m._callback_count +=1 
                m._n_subtour1_constraints += 1
                S = node_partition[0]  # 'left' side of the cut
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
        # m._callback_time += time.time()-initial
        m._time_subtour1_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def DFJ_smarter_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case1 and not case2:
            # m._callback_time += time.time()-initial
            m._time_subtour1_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial

            return
        
        # retrieve the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
        
        # if the support graph is disconnected, then finding violated cuts is easy!
        DG = m._DG
        edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]
        DG_support = DG.edge_subgraph( edges_used )
        if not nx.is_strongly_connected( DG_support ):
            for component in nx.strongly_connected_components( DG_support ):
                # m._callback_count +=1 
                m._n_subtour1_constraints += 1
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )
            m._time_subtour1_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
        else:
            # check if any CUT inequalities are violated (by solving some min cut problems)
            for i,j in DG.edges:
                DG.edges[i,j]['capacity'] = xval[i,j]

            s = 0
            for t in range(1,m._n):
                (cut_value, node_partition) = nx.minimum_cut( DG, _s=s, _t=t )
                # print("cut_value =",cut_value)
                if cut_value < 1 - m._epsilon:
                    S = node_partition[0]  # 'left' side of the cut
                    # m._callback_count +=1 
                    m._n_subtour1_constraints += 1
                    m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
                    # m._callback_time += time.time()-initial
                    m._time_subtour1_constraints += time.time()-initial
                    m._time_total_constraints += time.time()-initial
                    return
        
    @staticmethod
    def subtour_method(subruta, vals,n):
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
    
    @staticmethod
    def subtourelim1(modelo:gp.Model, donde):
        initial = time.time()
        n = modelo._n
        if donde == gp.GRB.Callback.SIMPLEX:
            modelo._simplex_lower_bound = modelo.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)

            if len(tour) < n:
                # modelo._callback_count +=1
                modelo._n_subtour1_constraints += 1
                tour2 = [i for i in range(n) if i not in tour]
                modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)
            modelo._time_subtour1_constraints += time.time()-initial
        # modelo._callback_time += time.time()-initial
        modelo._time_total_constraints += time.time()-initial

    @staticmethod
    def subtourelim2(modelo:gp.Model, donde):        
        initial = time.time()
        n = modelo._n
        if donde == gp.GRB.Callback.SIMPLEX:
            modelo._simplex_lower_bound = modelo.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        
        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)
            if len(tour) < n:
                # modelo._callback_count +=1 
                modelo._n_subtour1_constraints += 1
                tour2 = [i for i in range(n) if i not in tour]
                modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)

            valoresY = modelo.cbGetNodeRel(modelo._yvars)
            solucion = [(arco,solucion) for arco,solucion in valoresY.items() if solucion >0 and solucion <1]
            if len(solucion) >0:
                solucion = [(arco,solucion) for arco,solucion in valoresY.items() if solucion >0 ]
                new_lb = get_min_job(modelo._JT,solucion)
                if new_lb > modelo._jt_min:
                    modelo._jt_min = new_lb
                    #print(modelo._jt_min)
                    modelo.cbLazy(modelo._Cmax >= modelo._jt_min + gp.quicksum(modelo._xvars[i,j]*modelo._TT[i][j] for i in modelo._cities for j in modelo._cities[1:] if i!=j))
                    # modelo._callback_count +=1 
                    modelo._n_subtour1_constraints += 1


        # modelo._callback_time += time.time()-initial  
        modelo._time_subtour1_constraints += time.time()-initial
        modelo._time_total_constraints += time.time()-initial

    @staticmethod
    def subtourelim3(m:gp.Model, where):
        initial = time.time()
        case2 = (where == GRB.Callback.MIPNODE) and (m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL)
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        # obtains the LP solution
        if case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            return
            
        n = m._n
        tour = sc.SEC(xval, 0.00001, n) # using SECpy by C power
        m._time_subtour2_constraints += time.time()-initial
        if len(tour) > 0:
            tour2 = [i for i in range(n) if i not in tour]    
            m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
            # m._callback_count +=1
            # m._callback_time += time.time()-initial
            m._n_subtour2_constraints += 1
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def subtourelim4(m:gp.Model, where):
        initial = time.time()
        case1 = where == GRB.Callback.MIPSOL
        case2 = (where == GRB.Callback.MIPNODE) and (m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL)
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        if not case1 and not case2:
            m._time_subtour1_constraints += time.time()-initial
            m._time_subtour2_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        # obtains the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
        
        n = m._n
        tour = [i for i in range(n+1)]
        MathematicalModel.subtour_method(tour, xval, n)
        m._time_subtour1_constraints += time.time()-initial
        if len(tour) < n:
         # adds cut DFJ
            tour2 = [i for i in range(n) if i not in tour]
            m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
            # m._callback_count +=1
            m._n_subtour1_constraints += 1
        else:
            tour = sc.SEC(xval, 0.00001, n)
            m._time_subtour2_constraints += time.time()-initial
            if len(tour) == 0:
                m._time_total_constraints += time.time()-initial
                return

            tour2 = [i for i in range(n) if i not in tour]    
            m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
            # m._callback_count +=1
            m._n_subtour2_constraints += 1
        # m._callback_time += time.time()-initial
        m._time_total_constraints += time.time()-initial
 
    @staticmethod
    def subtourelim5(m:gp.Model, where):
        case1 = where == GRB.Callback.MIPSOL
        case2 = (where == GRB.Callback.MIPNODE) and (m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL)
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        initial = time.time()
        # obtains the LP solution
        if case1:
            xval = m.cbGetSolution(m._xvars)
        elif case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            m._time_subtour1_constraints += time.time()-initial
            m._time_subtour2_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        flag = False
        if case1:
            flag = True
            n = m._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
                m._n_subtour1_constraints += 1
                flag = False
            m._time_subtour1_constraints += time.time()-initial
        
        if case2 or flag:
            n = m._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
                m._n_subtour2_constraints += 1
            m._time_subtour2_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def subtourelim6(m:gp.Model, where):
        #case1 = where == GRB.Callback.MIPSOL
        case2 = (where == GRB.Callback.MIPNODE) and (m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL)
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        initial = time.time()
        # obtains the LP solution
        if case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            m._time_subtour1_constraints += time.time()-initial
            m._time_subtour2_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        flag = True
        if case2:
            #flag = True
            n = m._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
                m._n_subtour1_constraints += 1
                flag = False
            m._time_subtour1_constraints += time.time()-initial
        
        if flag:
            n = m._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
                m._n_subtour2_constraints += 1
            m._time_subtour2_constraints += time.time()-initial
        m._time_total_constraints += time.time()-initial

    @staticmethod
    def subtourelim7(m:gp.Model, where):
        #case1 = where == GRB.Callback.MIPSOL
        case2 = (where == GRB.Callback.MIPNODE) and (m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL)
        if where == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        initial = time.time()
        # obtains the LP solution
        if case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            m._time_subtour1_constraints += time.time()-initial
            m._time_subtour2_constraints += time.time()-initial
            m._time_total_constraints += time.time()-initial
            return
        flag = True     
        if case2:
            n = m._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
                m._n_subtour2_constraints += 1
                flag = False
            m._time_subtour2_constraints += time.time()-initial
        if flag:
            n = m._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
                m._n_subtour1_constraints += 1
                
            m._time_subtour1_constraints += time.time()-initial
            
        m._time_total_constraints += time.time()-initial
        
    @staticmethod
    def both_blossom_method(model:gp.Model, donde):
        initial_subtour = time.time()
        n = model._n
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        
        if donde == gp.GRB.Callback.SIMPLEX:
            model._simplex_lower_bound = model.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case2:
            # En caso de que no cumpla ninguna condición, el tiempo se agrega a todos los contadores
            model._time_subtour1_constraints += time.time()-initial_subtour
            model._time_blossom_heuristic_constraints += time.time()-initial_subtour
            model._time_blossom_exact_constraints += time.time()-initial_subtour
            model._time_total_constraints += time.time()-initial_subtour
            return
        
        xval = model.cbGetNodeRel(model._xvars)
        
        tour = [i for i in range(n+1)]
        MathematicalModel.subtour_method(tour, xval,n)
        # Se aumenta el contador de separación de subtour
        model._time_subtour1_constraints += time.time()-initial_subtour

        if len(tour) < n:
            tour2 = [i for i in range(n) if i not in tour]
            model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
            model._n_subtour1_constraints += 1

        # Blossom inequalities
        initial_heuristic_blossom = time.time()
        W,T = heuristic_separation(xval,n)
        # Se aumenta el contador de la separación heuristica
        model._time_blossom_heuristic_constraints += time.time()-initial_heuristic_blossom
        
        if len(T) >= 3:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in T)<=len(W)+sum([len(e) for e in T])-(3*len(T)+1)/2)
            model._n_blossom_heuristic_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        initial_exact_blossom = time.time()
        W,Fe = letchford_algorithm(model._G,xval)
        # Se aumenta el contador de separación exacta
        model._time_blossom_exact_constraints += time.time()-initial_exact_blossom
        
        if W != None and Fe != None:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            model._n_blossom_exact_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        #En caso que no se agreguen restricciones, se aumenta el contador global
        model._time_total_constraints += time.time() - initial_subtour
    
    @staticmethod
    def both_blossom_method2(model:gp.Model, donde):
        if model._stop_flag == True:
            model.terminate()
            return
        
        if model._relax == True and donde == gp.GRB.Callback.MESSAGE:
            gurobi_msj = model.cbGet(gp.GRB.Callback.MSG_STRING)
            if 'Root relaxation: objective' in gurobi_msj:
                objective_index = gurobi_msj.index('objective')
                seconds_index = gurobi_msj.index('seconds')
                objective_number = float(gurobi_msj[objective_index + len('objective'):seconds_index].strip().split(',')[0])
                seconds_number = float(gurobi_msj[seconds_index - 1:objective_index:-1].strip()[::-1].split(',')[-1])
                model._objective_root_relaxation = objective_number
                model._time_root_relaxation = seconds_number
                model._where_it_stopped = 'msj_r_rl'
                model._stop_flag = True
                model.terminate()
                return
            
            mensaje1 = ['Expl','Unexpl','Obj','Depth','IntInf','Incumbent','BestBd','Gap','It/Node','Time']
            condition = all([msj in gurobi_msj for msj in mensaje1])
            if condition :
                model._stop_flag = True
                model._where_it_stopped = 'msj_iter'
        
        initial_subtour = time.time()
        n = model._n
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        if case2:
            xval = model.cbGetNodeRel(model._xvars)
        else:
            model._time_subtour1_constraints += time.time()-initial_subtour
            model._time_subtour2_constraints += time.time()-initial_subtour
            model._time_total_constraints += time.time()-initial_subtour
            return
        flag = True     
        if case2:
            n = model._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                model.cbLazy(gp.quicksum(model._xvars[i,j] for i in tour for j in tour2) >= 1)
                model._n_subtour2_constraints += 1
                flag = False
            model._time_subtour2_constraints += time.time()-initial_subtour
        if flag:
            n = model._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
                model._n_subtour1_constraints += 1
                
            model._time_subtour1_constraints += time.time()-initial_subtour

        # Blossom inequalities
        initial_heuristic_blossom = time.time()
        W,T = heuristic_separation(xval,n)
        # Se aumenta el contador de la separación heuristica
        model._time_blossom_heuristic_constraints += time.time()-initial_heuristic_blossom
        
        if len(T) >= 3:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in T)<=len(W)+sum([len(e) for e in T])-(3*len(T)+1)/2)
            model._n_blossom_heuristic_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        initial_exact_blossom = time.time()
        W,Fe = letchford_algorithm(model._G,xval)
        # Se aumenta el contador de separación exacta
        model._time_blossom_exact_constraints += time.time()-initial_exact_blossom
        
        if W != None and Fe != None:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            model._n_blossom_exact_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        #En caso que no se agreguen restricciones, se aumenta el contador global
        model._time_total_constraints += time.time() - initial_subtour

    @staticmethod
    def both_blossom_method3(model:gp.Model, donde):
        initial_subtour = time.time()
        n = model._n
        #case1 = donde == gp.GRB.Callback.MIPSOL
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        
        if donde == gp.GRB.Callback.SIMPLEX:
            model._simplex_lower_bound = model.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if case2:
            xval = model.cbGetNodeRel(model._xvars)
        else:
            model._time_subtour1_constraints += time.time()-initial_subtour
            model._time_subtour2_constraints += time.time()-initial_subtour
            model._time_total_constraints += time.time()-initial_subtour
            return
        flag = True
        if case2:
            n = model._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
                model._n_subtour1_constraints += 1
                flag = False
            model._time_subtour1_constraints += time.time()-initial_subtour
        
        if flag:
            n = model._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                model.cbLazy(gp.quicksum(model._xvars[i,j] for i in tour for j in tour2) >= 1)
                model._n_subtour2_constraints += 1       
            model._time_subtour2_constraints += time.time()-initial_subtour

        # Blossom inequalities
        initial_heuristic_blossom = time.time()
        W,T = heuristic_separation(xval,n)
        # Se aumenta el contador de la separación heuristica
        model._time_blossom_heuristic_constraints += time.time()-initial_heuristic_blossom
        
        if len(T) >= 3:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in T)<=len(W)+sum([len(e) for e in T])-(3*len(T)+1)/2)
            model._n_blossom_heuristic_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        initial_exact_blossom = time.time()
        W,Fe = letchford_algorithm(model._G,xval)
        # Se aumenta el contador de separación exacta
        model._time_blossom_exact_constraints += time.time()-initial_exact_blossom
        
        if W != None and Fe != None:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                    gp.quicksum(model._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            model._n_blossom_exact_constraints += 1
            model._time_total_constraints += time.time() - initial_subtour
            return
        
        #En caso que no se agreguen restricciones, se aumenta el contador global
        model._time_total_constraints += time.time() - initial_subtour

    @staticmethod
    def exact_blossom_method(model:gp.Model, donde):
        initial_subtour = time.time()
        n = model._n
        #case1 = donde == gp.GRB.Callback.MIPSOL
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        
        if donde == gp.GRB.Callback.SIMPLEX:
            model._simplex_lower_bound = model.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case2:
            model._time_subtour1_constraints += time.time()-initial_subtour
            model._time_blossom_exact_constraints += time.time()-initial_subtour
            model._time_total_constraints += time.time()-initial_subtour
            return

        valoresX = model.cbGetNodeRel(model._xvars)
        
        tour = [i for i in range(n+1)]
        MathematicalModel.subtour_method(tour, valoresX,n)
        model._time_subtour1_constraints += time.time()-initial_subtour

        if len(tour) < n:
            tour2 = [i for i in range(n) if i not in tour]
            model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
            model._n_subtour1_constraints += 1

        #if case2:                    
        initial_blossom = time.time()
        W,Fe = letchford_algorithm(model._G,valoresX)
        model._time_blossom_exact_constraints += time.time()-initial_blossom
        if W != None and Fe != None:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                gp.quicksum(model._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            model._n_blossom_exact_constraints += 1
        
        model._time_total_constraints += time.time() - initial_subtour
            
    @staticmethod
    def exact_blossom_method2(m:gp.Model, donde):
        initial_subtour = time.time()
        n = m._n
        #case1 = donde == gp.GRB.Callback.MIPSOL
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( m.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        
        if donde == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            m._time_subtour1_constraints += time.time()-initial_subtour
            m._time_subtour2_constraints += time.time()-initial_subtour
            m._time_total_constraints += time.time()-initial_subtour
            return
        flag = True     
        if case2:
            n = m._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
                m._n_subtour2_constraints += 1
                flag = False
            m._time_subtour2_constraints += time.time()-initial_subtour
        if flag:
            n = m._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
                m._n_subtour1_constraints += 1
                
            m._time_subtour1_constraints += time.time()-initial_subtour
                               
        initial_blossom = time.time()
        W,Fe = letchford_algorithm(m._G,xval)
        m._time_blossom_exact_constraints += time.time()-initial_blossom
        if W != None and Fe != None:
            m.cbLazy(gp.quicksum(m._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                gp.quicksum(m._xvars[e] for e in Fe)<=len(W)+sum([len(e) for e in Fe])-(3*len(Fe)+1)/2)
            m._n_blossom_exact_constraints += 1
        
        m._time_total_constraints += time.time() - initial_subtour
        
    @staticmethod
    def heuristic_blossom_method(model:gp.Model, donde):
        initial_subtour = time.time()
        n = model._n
        #case1 = donde == gp.GRB.Callback.MIPSOL
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        
        if donde == gp.GRB.Callback.SIMPLEX:
            model._simplex_lower_bound = model.cbGet(gp.GRB.Callback.SPX_OBJVAL)

        if not case2:
            model._time_subtour1_constraints += time.time()-initial_subtour
            model._time_blossom_heuristic_constraints += time.time()-initial_subtour
            model._time_total_constraints += time.time()-initial_subtour
            return
        
        valoresX = model.cbGetNodeRel(model._xvars)
        
        tour = [i for i in range(n+1)]
        MathematicalModel.subtour_method(tour, valoresX,n)
        model._time_subtour1_constraints += time.time()-initial_subtour

        if len(tour) < n:
            tour2 = [i for i in range(n) if i not in tour]
            model.cbLazy(gp.quicksum(model._xvars[i, j] for i in tour for j in tour2) >= 1)
            model._n_subtour1_constraints += 1

        
        initial_blossom = time.time()
        W,T = heuristic_separation(valoresX,n)
        model._time_blossom_heuristic_constraints += time.time()-initial_blossom
        if len(T) >= 3:
            model.cbLazy(gp.quicksum(model._xvars[e] for e in [e for e in valoresX if valoresX[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                gp.quicksum(model._xvars[e] for e in T)<=len(W)+sum([len(e) for e in T])-(3*len(T)+1)/2)    
            model._n_blossom_heuristic_constraints += 1
            
        model._time_total_constraints += time.time() - initial_subtour
       #     return

    @staticmethod
    def heuristic_blossom_method2(m:gp.Model, donde):
        initial_subtour = time.time()
        n = m._n
        #case1 = donde == gp.GRB.Callback.MIPSOL
        case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( m.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
        if donde == gp.GRB.Callback.SIMPLEX:
            m._simplex_lower_bound = m.cbGet(gp.GRB.Callback.SPX_OBJVAL)
        if case2:
            xval = m.cbGetNodeRel(m._xvars)
        else:
            m._time_subtour1_constraints += time.time()-initial_subtour
            m._time_subtour2_constraints += time.time()-initial_subtour
            m._time_total_constraints += time.time()-initial_subtour
            return
        flag = True     
        if case2:
            n = m._n
            tour = sc.SEC(xval, 0.00001, n)
            if len(tour) > 0:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i,j] for i in tour for j in tour2) >= 1)
                m._n_subtour2_constraints += 1
                flag = False
            m._time_subtour2_constraints += time.time()-initial_subtour
        if flag:
            n = m._n
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, xval, n)
            
            if len(tour) < n:
                tour2 = [i for i in range(n) if i not in tour]
                m.cbLazy(gp.quicksum(m._xvars[i, j] for i in tour for j in tour2) >= 1)
                m._n_subtour1_constraints += 1
                
            m._time_subtour1_constraints += time.time()-initial_subtour
      
        initial_blossom = time.time()
        W,T = heuristic_separation(xval,n)
        m._time_blossom_heuristic_constraints += time.time()-initial_blossom
        if len(T) >= 3:
            m.cbLazy(gp.quicksum(m._xvars[e] for e in [e for e in xval if xval[e] > 0 and e[0] in W and e[1] in W and e[0]<e[1]] )+ 
                gp.quicksum(m._xvars[e] for e in T)<=len(T)+sum([len(e) for e in T])-(3*len(T)+1)/2)    
            m._n_blossom_heuristic_constraints += 1
            
        m._time_total_constraints += time.time() - initial_subtour

    @staticmethod
    def get_lower_bound_callback(model:gp.Model, donde):
        if model._stop_flag == True:
            model.terminate()
            return
        
        if model._relax == True and donde == gp.GRB.Callback.MESSAGE:
            gurobi_msj = model.cbGet(gp.GRB.Callback.MSG_STRING)
            if 'Root relaxation: objective' in gurobi_msj:
                objective_index = gurobi_msj.index('objective')
                seconds_index = gurobi_msj.index('seconds')
                objective_number = float(gurobi_msj[objective_index + len('objective'):seconds_index].strip().split(',')[0])
                seconds_number = float(gurobi_msj[seconds_index - 1:objective_index:-1].strip()[::-1].split(',')[-1])
                model._objective_root_relaxation = objective_number
                model._time_root_relaxation = seconds_number
                model._where_it_stopped = 'msj_r_rl'
                model._stop_flag = True
                model.terminate()
                return
            
            mensaje1 = ['Expl','Unexpl','Obj','Depth','IntInf','Incumbent','BestBd','Gap','It/Node','Time']
            condition = all([msj in gurobi_msj for msj in mensaje1])
            if condition :
                model._stop_flag = True
                model._where_it_stopped = 'msj_iter'
   
    def add_new_constraint(self):
        """
        Add new constraints, builts in this work.
        """
        
        menor_arco_depot = min(self.TT[0][i] for i in range(1,self.n))
        
        if hasattr(self,'TS'):
            for i in self.cities:
                if self.exp_lb == 0:
                    self.modelo.addConstr(self.TS[i] <= self.initial_route_fitness , name = f'bounds1_TS_{i}')
                elif self.exp_lb == 1:
                    self.modelo.addConstr(self.TS[i] <= self.initial_route_fitness - self.TT[self.last_node_initial_route][0] , name = f'bounds1_TS_{i}')
                elif self.exp_lb == 2:
                    self.modelo.addConstr(self.TS[i] <= self.initial_route_fitness - menor_arco_depot , name = f'bounds1_TS_{i}')

            for i in range(1,self.n):
                self.modelo.addConstr(self.TS[i] >= menor_arco_depot , name = f'bounds2_TS_{i}') 
        
        elif hasattr(self,'t'):
            for i in self.cities:
                LB = np.min(self.TT[i][np.nonzero(self.TT[i])])
                for j in self.cities:
                    if i!=j:
                        if self.exp_lb == 0:

                            self.modelo.addConstr(self.t[(i,j)] <= self.initial_route_fitness * self.x[(i,j)], name = f'bounds1_t_{i}_{j}')
                        elif self.exp_lb == 1:
                            self.modelo.addConstr(self.t[(i,j)] <= (self.initial_route_fitness-self.TT[self.last_node_initial_route][0]) * self.x[(i,j)], name = f'bounds1_t_{i}_{j}')
                        elif self.exp_lb == 2:
                            self.modelo.addConstr(self.t[(i,j)] <= (self.initial_route_fitness - menor_arco_depot ) * self.x[(i,j)], name = f'bounds1_t_{i}_{j}')
                        
                    if i!=j and i>0:
                        self.modelo.addConstr(self.t[(i,j)] >= LB*self.x[(i,j)] , name = f'bouns3_t_{i}_{j}')
     
        # Cota superior de solucion inicial (lkh+NNJ)
        self.modelo.addConstr(self.Cmax<=self.initial_fitness)
        # Cota inferior de solucion inicial
        self.modelo.addConstr(self.Cmax>= self.sum_min_row + self.jt_min)
        
        for i in self.cities:
            for j in self.cities:
                if i!=j:
                    self.modelo.addConstr(self.x[(i,j)] + self.x[(j,i)] <= 1 , name = f'bouns2_t_{i}_{j}')
        self.modelo.update()

    def add_initial_solution(self):
        """
        Add initial solution to MILP from LKH and NNJA
        """

        self.modelo.NumStart = 1
        self.modelo.update()
        for s in range(self.modelo.NumStart):
            self.modelo.Params.StartNumber = s
        
            for var in self.modelo.getVars(): 
                if var.VarName[0] == "x":
                    arch_name = tuple(var.varName.split("[")[1][:-1].split(","))
                    if arch_name in self.initial_arch:
                        var.Start = 1

                elif var.VarName[0] == "y":
                    arch_name = tuple(var.varName.split("[")[1][:-1].split(","))
                    if arch_name in self.initial_job_arch:
                        var.Start = 1
        self.modelo.setParam("Cutoff",self.initial_fitness) 
        self.modelo.update()

    def optimize(self):
        self.modelo.update()
        if self.callback in ("none",None) :
            self.modelo.update()
            self.modelo.optimize(MathematicalModel.get_lower_bound_callback)
        elif self.callback in self.callback_dict.keys():
            self.modelo.Params.LazyConstraints = 1
            self.modelo._xvars = self.x
            self.modelo._yvars = self.y
            self.modelo._JT = self.JT
            self.modelo._TT = self.TT
            self.modelo._Cmax = self.Cmax
            self.modelo._jt_min = self.jt_min
            self.modelo._cities = self.cities
            self.modelo._n = self.n
            self.modelo._coords = self.coords
            
            if "blossom" in self.callback :
                self.modelo._G = nx.Graph(nx.complete_graph(self.n))
            elif 'separation' in self.callback:
                self.modelo._epsilon = 0.00001
                self.modelo._DG = nx.DiGraph(nx.complete_graph(self.n))
            else:
                self.modelo._DG = None
            
            try:
                self.modelo.update()
                self.modelo.optimize(self.callback_dict[self.callback])
            except Exception as e:
                raise Exception("Callback error") from e
        else:
            raise Exception(f"Callback {self.callback} not found")

    def run(self):
        self.create_new_formulation() if self.new_formulation else self.create_base_model()
        self.add_subtour_constraint()
        if self.bounds:
            self.add_new_constraint()
        
        if self.initial_solution:
            self.add_initial_solution()

        self.modelo.Params.TimeLimit = max(self.time_limit,self.time_limit - self.initial_solution_time)
        self.optimize()
        self.modelo.update()

    def get_integer_results(self):
        if self.modelo.Status == GRB.OPTIMAL or self.modelo.SolCount > 0:
            objective = round(self.modelo.ObjVal, 2)
            lower = round(self.modelo.ObjBound, 4)
            gap = round(self.modelo.MIPGap * 100, 2)
        else:
            objective = round(getattr(self.modelo, 'ObjVal', self.initial_fitness), 2)
            lower = round(getattr(self.modelo, 'ObjBound', self.modelo.ObjBoundC), 4)
            gap = round((objective - lower) / lower * 100, 4)

        return objective, lower, gap

    def get_lb_results(self):
        if not hasattr(self.modelo,'_where_it_stopped'):
            lower = round(self.modelo._simplex_lower_bound, 4)
            self.modelo._where_it_stopped = 'NO_LB'

        if self.modelo._where_it_stopped != 'msj_iter':
            lower = round(self.modelo._objective_root_relaxation, 4)
        else:
            try:
                lower = round(getattr(self.modelo, 'ObjBound', self.modelo.ObjBoundC), 2)
            except AttributeError:
                lower = round(self.modelo._objective_root_relaxation, 2)
        objective = float('inf')
        gap = float('inf')
        return objective, lower, gap

    def print_results(self):
        objective, lower, gap = self.get_lb_results() if self.relax else self.get_integer_results()

        dict_status = {
            1: 'LOADED', 2: 'OPTIMAL', 3: 'INFEASIBLE', 4: 'INF_OR_UNBD',
            5: 'UNBOUNDED', 6: 'CUTOFF', 7: 'ITERATION_LIMIT', 8: 'NODE_LIMIT',
            9: 'TIME_LIMIT', 10: 'SOLUTION_LIMIT', 11: 'INTERRUPTED', 12: 'NUMERIC',
            13: 'SUBOPTIMAL', 14: 'INPROGRESS', 15: 'USER_OBJ_LIMIT'
        }
        
        time = round(self.modelo.Runtime + self.initial_solution_time, 4)
        time_subtour1 = round(self.modelo._time_subtour1_constraints, 4)
        time_subtour2 = round(self.modelo._time_subtour2_constraints, 4)
        time_heuristic_blossom = round(self.modelo._time_blossom_heuristic_constraints, 4)
        time_exact_blossom = round(self.modelo._time_blossom_exact_constraints, 4)
        time_total_callback = round(self.modelo._time_total_constraints, 4)

        status = dict_status[self.modelo.Status]
        if status == 'OPTIMAL':
            gap = 0
            lower = objective

        if self.relax == True and self.modelo.Status == GRB.INTERRUPTED:
            status = self.modelo._where_it_stopped.upper()
        
        lower = round(lower, 2)
        print("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format(
            self.size, self.instance, objective, lower, gap, time, status,
            int(self.modelo.NodeCount), time_subtour1, time_subtour2, time_heuristic_blossom,
            time_exact_blossom, time_total_callback, self.modelo._n_subtour1_constraints,
            self.modelo._n_subtour2_constraints, self.modelo._n_blossom_heuristic_constraints,
            self.modelo._n_blossom_exact_constraints
        ))

    def old_print_results(self):
        dict_status = {1: 'LOADED', 2: 'OPTIMAL', 3: 'INFEASIBLE', 4: 'INF_OR_UNBD', 5: 'UNBOUNDED', 6: 'CUTOFF', 7: 'ITERATION_LIMIT', 8: 'NODE_LIMIT', 9: 'TIME_LIMIT', 10: 'SOLUTION_LIMIT', 11: 'INTERRUPTED', 12: 'NUMERIC', 13: 'SUBOPTIMAL', 14: 'INPROGRESS', 15: 'USER_OBJ_LIMIT'}
        try:
            lower = self.modelo.ObjBound
        except AttributeError:
            lower = self.modelo._simplex_lower_bound
        if self.modelo.Status == GRB.OPTIMAL or self.modelo.SolCount > 0:
            objective = self.modelo.getObjective().getValue()
            gap = round((objective-lower)/lower*100,4)
            lower = round(lower,2)
            objective = round(objective,2)

        else:
            lower = round(lower,2)
            objective = round(self.initial_fitness,2)
            gap = round((objective-lower)/lower*100,4)

        time = round(self.modelo.Runtime + self.initial_solution_time,4) 
        time_subtour1 = round(self.modelo._time_subtour1_constraints,4)
        time_subtour2 = round(self.modelo._time_subtour2_constraints,4)
        time_heuristic_blossom = round(self.modelo._time_blossom_heuristic_constraints,4)
        time_exact_blossom = round(self.modelo._time_blossom_exact_constraints,4)
        time_total_callback = round(self.modelo._time_total_constraints,4)
        if dict_status[self.modelo.Status] == 'OPTIMAL':
            gap = 0

        lower = round(lower,2)
        print("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format(
            self.size,
            self.instance,
            objective,
            lower,
            gap,
            time,
            dict_status[self.modelo.Status],
            int(self.modelo.NodeCount),
            time_subtour1,
            time_subtour2,
            time_heuristic_blossom,
            time_exact_blossom,
            time_total_callback,
            self.modelo._n_subtour1_constraints,
            self.modelo._n_subtour2_constraints,
            self.modelo._n_blossom_heuristic_constraints,
            self.modelo._n_blossom_exact_constraints,
        )
        )
        #print("{}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}".format(self.size,self.instance,objective,lower,gap,time,dict_status[self.modelo.Status],self.modelo.NodeCount,self.modelo._callback_count,self.modelo._callback_time))

    def get_solution(self):
        tsp_sol =  [key for key,value in self.modelo.getAttr('x', self.x).items() if value>0]
        actual = 0
        tsp_solution = []
        while len(tsp_solution) < self.n:
            tsp_solution.append(actual)
            for i in tsp_sol:
                if i[0]==actual and i[1] not in tsp_solution:
                    _next = i[1]
                    actual = _next
                    break
        
        job_sol = {key[0]:key[1] for key,value in self.modelo.getAttr('x', self.y).items() if value>0}
        job_sol = [0]+[job_sol[i] for i in tsp_solution[1:]]

        # t_values = {key:value for key,value in self.modelo.getAttr('x', self.t).items()}
        # x_values = {key:value for key,value in self.modelo.getAttr('x', self.x).items()}
        # print(tsp_sol)
        # for key in tsp_sol:
        #     print(key,t_values[key],x_values[key])

        # exit(0)

        return tsp_solution,job_sol
    
    def print_solution(self):
        print(self.get_solution())
