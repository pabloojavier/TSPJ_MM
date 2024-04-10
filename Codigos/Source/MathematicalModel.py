import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import numpy.typing as npt
import warnings
import networkx as nx
from typing import Dict, List, Callable, Any
from Source.Problem import Problem 
import time
import re
import math
import matplotlib.pyplot as plt
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

def padberg_rao(valoresX):
    weights = {(u[0],u[1]):v  for u,v in valoresX.items()}
    capacity = {(u[0],u[1]):min(v,1-v) for u,v in valoresX.items()}

    # fig,ax = plt.subplots(1,4,figsize = (15,6))
    
    G = nx.Graph()
    for e,v in weights.items():
        if v > 0:
            G.add_edge(e[0],e[1],weight = v)
    
    H = nx.Graph()
    for e,v in capacity.items():
        if v > 0:
            H.add_edge(e[0],e[1],capacity = v,weight = weights[e])
        else:
            H.add_nodes_from([e[0],e[1]])

    # pos = nx.spring_layout(G)
    
    # g_edge_label = {k: round(v, 3) for k, v in nx.get_edge_attributes(G, 'weight').items() if v > 0}
    # nx.draw(G,ax = ax[0],with_labels=True,pos = pos)
    # nx.draw_networkx_edge_labels(G, pos, edge_labels=g_edge_label,ax=ax[0],font_size=6)
    # ax[0].set_title('Original Graph\nwith weights as x*')

    # h_edge_label = {k: round(v, 3) for k, v in nx.get_edge_attributes(H, 'capacity').items() if v > 0 }
    # nx.draw(H,ax = ax[1],with_labels=True,pos = pos)
    # nx.draw_networkx_edge_labels(H, pos, edge_labels=h_edge_label,ax=ax[1],font_size=6)
    # ax[1].set_title('H Graph\nwith capacities as min(x*,1-x*)')
    
    gomHu = nx.gomory_hu_tree(H, capacity = 'capacity')
    # nx.draw(gomHu,ax = ax[2],with_labels = True,pos = pos)
    # ax[2].set_title('Gomory-Hu Tree\nfrom H Graph')

    for e in gomHu.edges():
        # ax[3].cla()
        Te = nx.Graph(gomHu)
        Te.remove_edge(e[0],e[1])
        U,V = list(nx.connected_components(Te))
        # nx.draw(Te,ax = ax[3],with_labels = True,pos = pos)
        # nx.draw_networkx_edges(Te, ax=ax[3], pos=pos, edgelist=[e], edge_color='r',style='dashed')
        cutset = delta(H,U)

        Fe = set([e for e in cutset if 1-H[e[0]][e[1]]['weight'] < H[e[0]][e[1]]['weight']])
        if (len(U) + len(Fe))%2 == 0 and len(cutset)>0 :
            Faux = sorted([(term(edge,H),edge) for edge in cutset])
            fp = set()
            fp.add(Faux[0][1])
            Fe = Fe.symmetric_difference(fp)

        fig = None
        des = totalcap(H,cutset.difference(Fe)) + len(Fe) - totalcap(H,Fe)
        if des<1 and len(cutset)>0 and len(Fe)>0:
            # print('constraint:',cutset.difference(Fe),Fe,len(Fe))
            # ax[3].set_title(f'Edge {e} removed from Gomory-Hu Tree\nU:{U}\n V:{V}\n Cutset:{cutset}-{Fe}')
            return (cutset,Fe,fig)
            
    return (None,None,fig)

class MathematicalModel(Problem):
    def __init__(self,size:str,
                 instance,
                 output = False,
                 subtour : str = "GG",
                 initial_solution : bool = True,
                 callback : str = "None",
                 bounds: bool = True,
                 new_formulation: bool = False,
                 time_limit : int = 1800,
                 new_m : bool = False,
                 relax : bool = False
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
        self.callback_dict = {"cut_integer_separation":MathematicalModel.CUT_integer_separation,
                              "cut_naive_fractional_separation":MathematicalModel.CUT_naive_fractional_separation,
                              "cut_smarter_fractional_separation":MathematicalModel.CUT_smarter_fractional_separation,
                              "dfj_integer_separation":MathematicalModel.DFJ_integer_separation,
                              "dfj_naive_fractional_separation":MathematicalModel.DFJ_naive_fractional_separation,
                              "dfj_smarter_fractional_separation":MathematicalModel.DFJ_smarter_fractional_separation,
                              "subtourelim1":MathematicalModel.subtourelim1,
                              "subtourelim_gomhu":MathematicalModel.subtourelim_gomhu,
                              "subtourelim2":MathematicalModel.subtourelim2,
                              #'integer_fractional_cut':MathematicalModel.integer_fractional_cut,
                              'custom':callback}
        if self.output:
            print(f'running with: {self.size} {self.instance} {self.subtour} {self.initial_solution} {self.callback} {self.bounds} {self.new_formulation} {self.time_limit} {self.new_m}')
        self.compute_new_M() if new_m else self.compute_M()
        self.jobs = self.cities.copy()
        self.jobs_arch = [(i,k) for i in self.cities for k in self.cities]
        self.heuristic_jobs = NNJA(self.lkh_route[1:],self.JT)
        self.initial_fitness = self.fitness_functions([self.lkh_route[1:],self.heuristic_jobs])[0]
        self.initial_solution_time = 0
        
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
        self.M = np.zeros((len(self.cities),len(self.cities)))
        route_fitness = self.route_fitness(self.lkh_route[1:])
        self.sum_min_row = 0
        # for j in range(1,len(self.cities)):
        #     self.M[0,j] = self.TT[0][j]
        for i in range(len(self.cities)):
            self.sum_min_row += np.min(self.TT[i][np.nonzero(self.TT[i])])
            for j in range(len(self.cities)):
                if i != j:
                    self.M[i,j] = route_fitness + self.TT[i][j]

    def create_base_model(self):
        """
        Create the initial MILP, whithout any additional constraint.
        """
        env = gp.Env(empty=True)
        env.setParam('OutputFlag', 1 if self.output else 0)
        env.start()
        self.modelo = gp.Model(env=env)
        self.modelo._callback_count = 0
        self.modelo._callback_time = 0

        self.Cmax = self.modelo.addVar(vtype=GRB.CONTINUOUS,name="Cmax")

        var_mode = GRB.CONTINUOUS if self.relax else GRB.BINARY
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
        self.modelo._callback_count = 0
        self.modelo._callback_time = 0

        self.Cmax = self.modelo.addVar(vtype=GRB.CONTINUOUS,name="Cmax")

        var_mode = GRB.CONTINUOUS if self.relax else GRB.BINARY
        self.x = self.modelo.addVars(self.arch, vtype=var_mode, name='x')
        self.y = self.modelo.addVars(self.jobs_arch, vtype=var_mode, name='y')

        self.t = self.modelo.addVars(self.arch,vtype=GRB.CONTINUOUS,name="t") # new variable for the time

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
                                + gp.quicksum(self.TT[k][i]*self.x[(i,k)] for i in self.cities if i != k)
                                <= gp.quicksum(self.t[(k,l)] for l in self.cities if k != l) , name = f't_{k}')

        #parte desde 1, el primer nodo no tiene tiempo
        for i in self.cities[1:]:
            LB = np.min(self.TT[i][np.nonzero(self.TT[i])])
            for k in self.cities:
                if i != k:
                    self.modelo.addConstr(self.t[(i,k)] <= self.M[i][k]*self.x[(i,k)] , name = f't_{i}_{k}_UB')
                    self.modelo.addConstr(self.t[(i,k)] >= 0 , name = f't_{i}_{k}_LB')
        
        self.modelo.Params.Threads = 1
        self.modelo._callback_time = 0
        #0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent
        self.modelo.setParam("Method",2) 
        self.modelo.setParam("Cutoff",self.initial_fitness) 
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

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)
            if len(tour) < n:
                modelo._callback_count +=1 
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

        modelo._callback_time += time.time()-initial
    
    @staticmethod
    def CUT_integer_separation(m:gp.Model, where):
        initial = time.time()
        # check if LP relaxation at this branch-and-bound node has an integer solution
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
                    m._callback_count +=1 
                    complement = [ i for i in DG_soln.nodes if i not in component ]
                    m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in complement ) >= 1 )
        m._callback_time += time.time()-initial

    @staticmethod
    def CUT_naive_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if not case1 and not case2:
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
                m._callback_count +=1 
                S = node_partition[0]  # 'left' side of the cut
                T = node_partition[1]  # 'right' side of the cut
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in T ) >= 1 )
                m._callback_time += time.time()-initial
                return

    @staticmethod
    def CUT_smarter_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        n = m._n
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if not case1 and not case2:
            m._callback_time += time.time()-initial
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
                m._callback_count +=1 
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
                    m._callback_time += time.time()-initial
                    return

    @staticmethod
    def DFJ_integer_separation(m:gp.Model, where):
        initial = time.time()
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
                    m._callback_count +=1 
                    m.cbLazy( gp.quicksum( m._xvars[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )
        m._callback_time += time.time()-initial

    @staticmethod
    def DFJ_naive_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if not case1 and not case2:
            m._callback_time += time.time()-initial
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
                m._callback_count +=1 
                S = node_partition[0]  # 'left' side of the cut
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
        m._callback_time += time.time()-initial

    @staticmethod
    def DFJ_smarter_fractional_separation(m:gp.Model, where):
        initial = time.time()
        # We must *separately* handle these two cases of interest:
        #    1. We encounter an integer point that might replace our incumbent
        #    2. We encounter a fractional point *that is LP optimal*
        #
        case1 = where == GRB.Callback.MIPSOL
        case2 = ( where == GRB.Callback.MIPNODE ) and ( m.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL )
        
        if not case1 and not case2:
            m._callback_time += time.time()-initial
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
                m._callback_count +=1 
                m.cbLazy( gp.quicksum( m._xvars[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )
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
                    m._callback_count +=1 
                    m.cbLazy( gp.quicksum( m._xvars[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
                    m._callback_time += time.time()-initial
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

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)
            if len(tour) < n:
                modelo._callback_count +=1 
                tour2 = [i for i in range(n) if i not in tour]
                modelo.cbLazy(gp.quicksum(modelo._xvars[i, j] for i in tour for j in tour2) >= 1)

        modelo._callback_time += time.time()-initial

    @staticmethod
    def subtourelim_gomhu(modelo:gp.Model, donde):
        initial = time.time()
        n = modelo._n

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)

            cutset,Fe,fig = padberg_rao(valoresX)
            if cutset != None and Fe != None:
                modelo._callback_count +=1
                modelo.cbLazy(gp.quicksum(modelo._xvars[e] for e in cutset.difference(Fe))-
                              gp.quicksum(modelo._xvars[e] for e in Fe) 
                              >= 1-len(Fe))

        modelo._callback_time += time.time()-initial

    @staticmethod
    def subtourelim2(modelo:gp.Model, donde):        
        initial = time.time()
        n = modelo._n

        if donde == GRB.Callback.MIPNODE  and ( modelo.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL):
            valoresX = modelo.cbGetNodeRel(modelo._xvars)
            tour = [i for i in range(n+1)]
            MathematicalModel.subtour_method(tour, valoresX,n)
            if len(tour) < n:
                modelo._callback_count +=1 
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
                    modelo._callback_count +=1 

        modelo._callback_time += time.time()-initial    
    
    def add_new_constraint(self):
        """
        Add new constraints, builts in this work.
        """
        
        route_fitness = self.route_fitness(self.lkh_route[1:])
        menor_arco_depot = min(self.TT[0][i] for i in range(1,self.n))

        if hasattr(self,'TS'):
            for i in self.cities:
                self.modelo.addConstr(self.TS[i] <= route_fitness , name = f'bounds1_TS_{i}')

            for i in range(1,self.n):
                self.modelo.addConstr(self.TS[i] >= menor_arco_depot , name = f'bounds2_TS_{i}') 
        
        elif hasattr(self,'t'):
            for i in self.cities:
                LB = np.min(self.TT[i][np.nonzero(self.TT[i])])
                for j in self.cities:
                    if i!=j:
                        self.modelo.addConstr(self.t[(i,j)] <= route_fitness, name = f'bounds1_t_{i}_{j}')            
                        self.modelo.addConstr(self.x[(i,j)] + self.x[(j,i)] <= 1 , name = f'bouns2_t_{i}_{j}')
                    if i!=j and i*j>0:
                        self.modelo.addConstr(self.t[(i,j)] >= LB*self.x[(i,j)] , name = f'bouns3_t_{i}_{j}')
                    
                
        # Cota superior de solucion inicial (lkh+NNJ)
        self.modelo.addConstr(self.Cmax<=self.initial_fitness)
        # Cota inferior de solucion inicial
        self.modelo.addConstr(self.Cmax>= self.sum_min_row + self.jt_min)
    
        self.modelo.update()

    def add_initial_solution(self):
        """
        Add initial solution to MILP from LKH and NNJA
        """
        self.initial_solution_time = time.time()
        self.initial_arch = sort_arch(self.lkh_route)
        if self.heuristic_jobs is None:
            self.heuristic_jobs = NNJA(self.lkh_route[1:],self.JT) 
        self.initial_job_arch = sort_jobs(self.lkh_route[1:],self.heuristic_jobs)
        
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
        self.initial_solution_time = time.time()-self.initial_solution_time
        self.modelo.update()

    def optimize(self):
        self.modelo.update()
        if self.callback in ("none",None) :
            self.modelo.optimize()
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
            
            self.modelo._epsilon = 0.00001
            self.modelo._DG = nx.DiGraph(nx.complete_graph(self.n))
            
            try:
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

        self.modelo.Params.TimeLimit = self.time_limit - self.initial_solution_time
        self.optimize()
        self.modelo.update()

        # dict_values = {}
        # for v in self.modelo.getVars():
        #     if v.VarName[0] in ("y"):
                
        #         if v.X > 0.5:
        #             match = re.search(r'\[(\d+),(\d+)\]', v.VarName)
        #             number1 = int(match.group(1))
        #             number2 = int(match.group(2))
        #             #if number1 == 0 and number2 == 15:
        #             dict_values[v.VarName] = math.ceil(v.X)
        
        # for i in sorted(dict_values, key=dict_values.get):
        #     print(i,dict_values[i])

    def print_results(self):
        dict_status = {1: 'LOADED', 2: 'OPTIMAL', 3: 'INFEASIBLE', 4: 'INF_OR_UNBD', 5: 'UNBOUNDED', 6: 'CUTOFF', 7: 'ITERATION_LIMIT', 8: 'NODE_LIMIT', 9: 'TIME_LIMIT', 10: 'SOLUTION_LIMIT', 11: 'INTERRUPTED', 12: 'NUMERIC', 13: 'SUBOPTIMAL', 14: 'INPROGRESS', 15: 'USER_OBJ_LIMIT'}
        try:
            lower = self.modelo.ObjBoundC
        except AttributeError:
            lower = self.modelo.ObjBound
        objective = float("inf")
        gap =  float("inf")
        if self.modelo.Status == GRB.OPTIMAL or self.modelo.SolCount > 0:
            objective = self.modelo.getObjective().getValue()
            gap = round((objective-lower)/lower*100,4)
            lower = round(lower,2)
            objective = round(objective,2)

        else:
            #print(self.instance, ':Optimization ended with status %d' % self.modelo.Status)
            if self.modelo.SolCount > 0:
                objective = self.modelo.getObjective().getValue()
                lower = self.modelo.getObjective().getValue()
                gap = round((objective-lower)/lower*100,4)
                lower = round(lower,2)
                objective = round(objective,2)
        
        time = round(self.modelo.Runtime + self.initial_solution_time,4) 
        lower = round(lower,2)
        print("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<15}{:<10}{:<10}{:<10}".format(self.size,self.instance,objective,lower,gap,time,dict_status[self.modelo.Status],self.modelo.NodeCount,self.modelo._callback_count,self.modelo._callback_time))
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
        return tsp_solution,job_sol
    
    def print_solution(self):
        print(self.get_solution())