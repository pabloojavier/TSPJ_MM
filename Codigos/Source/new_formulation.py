import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np
import os

if os.getcwd().split("/")[-1] != "Codigos":
    path = "Codigos/"
else:
    path = ""

def parameters(size,instance):
    
    batch = (int(instance)-1)//25+1 if size in ["small","medium","large"] else ""
    if size == "tsplib":
        location = path+"Data/Tsplib_problems/"
        TT = pd.read_csv(location+"TT_"+instance+".csv",index_col= None, header = None).fillna(0).to_numpy()
        JT = pd.read_csv(location+"JT_"+instance+".csv",index_col= None, header = None).fillna(0).to_numpy()

    elif size in ("small","medium","large"):
        location = path+"Data/"+str(size.capitalize())+"_problems/Batch_0"+str(batch)+"/TSPJ_"+str(instance)+size.capitalize()[0]
        TT = pd.read_csv(location+"_cost_table_by_coordinates.csv",index_col= None, header = None).fillna(0).to_numpy()
        JT = pd.read_csv(location+"_tasktime_table.csv"           ,index_col= None, header = None).fillna(0).to_numpy()

    else:
        TT = pd.read_csv(f"{path}Data/test/1_TT_paper.csv",index_col= None, header = None).fillna(0).to_numpy()
        JT = pd.read_csv(f"{path}Data/test/1_JT_paper.csv",index_col= None, header = None).fillna(0).to_numpy()

    n = len(TT)
    cities = [i for i in range(n)]
    arch = [(i,j) for i in cities for j in cities if i !=j]
    return n,cities,arch,TT,JT

instance = 1
size = "small"

n,cities,arch,TT,JT = parameters(size,str(instance))
jobs = cities.copy()
jobs_arch = [(i,k) for i in cities for k in cities]
jt_min = JT[JT>0].min()


M = 0
for i in range(n):
    max_t = 0
    max_tt = 0
    for j in range(n):
        if i!=j and TT[i][j]>max_t:
            max_t = TT[i][j]
        
        if j!= 0 and JT[i][j]>max_tt:
            max_tt = JT[i][j]
    M += max_t #+max_tt


env = gp.Env(empty=True)
env.setParam('OutputFlag', 1)
env.start()
modelo = gp.Model(env=env)

Cmax = modelo.addVar(vtype=GRB.CONTINUOUS,name="Cmax")

x = modelo.addVars(arch, vtype=GRB.BINARY, name='x')
z = modelo.addVars(jobs_arch, vtype=GRB.BINARY, name='z')

t = modelo.addVars(arch,vtype=GRB.CONTINUOUS,name="t") # new variable for the time

modelo.setObjective(Cmax, GRB.MINIMIZE)

#Restricción adicional de reforzamiento
modelo.addConstr(Cmax >= jt_min + gp.quicksum(x[(i,j)]*TT[i][j] for i in cities for j in cities[1:] if i!=j), name = "Reinforcement")

for i in cities[1:n]:
    modelo.addConstr(Cmax >= gp.quicksum(t[(i,k)] for k in cities if i != k) 
                                        + gp.quicksum(z[(i,k)]*JT[i][k] for k in jobs if k!=0 ) , name = f'Cmax_{i}')

for i in cities[1:n]:
    modelo.addConstr(Cmax >= gp.quicksum(t[(i,k)] for k in cities if i != k) 
                                        + x[(i,0)]*TT[0][i] , name = f'Cmax_{i}_0')

for k in jobs[1:n]:
    modelo.addConstr(gp.quicksum(z[(i,k)] for i in cities if i != 0) == 1 , name = f'Job_{k}_out')

for i in cities[1:n]: 
    modelo.addConstr(gp.quicksum(z[(i,k)] for k in jobs if k != 0) == 1 , name = f'Job_{i}_in')

modelo.addConstrs((x.sum(i,'*') == 1 for i in cities) , name = 'Outgoing') # Outgoing
modelo.addConstrs((x.sum('*', j) == 1 for j in cities) , name = 'Incoming') # Incoming 

for k in cities[1:n]: #16
    modelo.addConstr(gp.quicksum(t[(i,k)] for i in cities if i != k) 
                        + gp.quicksum(TT[k][i]*x[(i,k)] for i in cities if i != k)
                        <= gp.quicksum(t[(k,l)] for l in cities if k != l) , name = f't_{k}')

#parte desde 1, el primer nodo no tiene tiempo
for i in cities[1:]:
    LB = np.min(TT[i][np.nonzero(TT[i])])
    for k in cities:
        if i != k:
            modelo.addConstr(t[(i,k)] <= M*x[(i,k)] , name = f't_{i}_{k}_UB')
            modelo.addConstr(t[(i,k)] >= LB*x[(i,k)] , name = f't_{i}_{k}_LB')    

modelo.Params.Threads = 1
modelo.optimize()
print(size,instance,modelo.getObjective().getValue())