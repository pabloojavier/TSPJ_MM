import numpy as np
import random
import gurobipy as gp
from gurobipy import GRB
import os
SEED = 202
np.random.seed(SEED)
random.seed(SEED)

def export_data_to_csv(data, filename):
    np.savetxt(filename, data, delimiter=',', fmt='%10.0f')

def import_data_from_csv(filename)->np.ndarray:
    return np.genfromtxt(filename, delimiter=',')

def node_table_by_coordinates(node_number, lower_number, upper_number, csv_fname):
    data = lower_number+upper_number*np.random.rand((node_number*2),2)
    data = data.astype(int)
    a = np.ascontiguousarray(data)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    data = unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    data = data[:node_number]

    export_data_to_csv(data, csv_fname)
    
    data = open(csv_fname, 'r').read()
    os.remove(csv_fname)
    archivo = open(csv_fname, 'w')
    archivo.write(data.replace('nan','').replace(' ',''))
    archivo.close()

def cost_table_by_coordinates(csv_fname, csv_fname_out):
    data = import_data_from_csv(csv_fname).astype(int)

    x = np.asarray([x[0] for x in data],dtype=int)
    y = np.asarray([x[1] for x in data],dtype=int)

    dx = np.sqrt(
        (x[...,np.newaxis] - x[np.newaxis,...])**2+
        (y[...,np.newaxis] - y[np.newaxis,...])**2
    ).astype(float)
    export_data_to_csv(dx, csv_fname_out)
    data = open(csv_fname_out, 'r').read()
    os.remove(csv_fname_out)
    archivo = open(csv_fname_out, 'w')
    archivo.write(data.replace('nan','').replace(' ',''))
    archivo.close()

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

def subtourelim(model:gp.Model, donde):
    n = model._n
    case1 = donde == gp.GRB.Callback.MIPSOL
    case2 = ( donde == gp.GRB.Callback.MIPNODE ) and ( model.cbGet(gp.GRB.Callback.MIPNODE_STATUS) == gp.GRB.OPTIMAL )
    
    if not case1 and not case2:
        return
    
    # retrieve the LP solution
    if case1:
        valoresX = model.cbGetSolution(model._vars)
    elif case2:
        valoresX = model.cbGetNodeRel(model._vars)

    tour = [i for i in range(n+1)]
    subtour_method_node(tour, valoresX,n)

    if len(tour) < n:
        tour2 = [i for i in range(n) if i not in tour]
        model.cbLazy(gp.quicksum(model._vars[i, j] for i in tour for j in tour2) >= 1)

def TSP(filename):
    # Import cost table from CSV file
    cost_table = import_data_from_csv(filename)
    n = len(cost_table[0])
    # Create a dictionary from the cost table
    cost_dict = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                cost_dict[(i, j)] = cost_table[i, j]
    # Create model

    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()
    m = gp.Model(env=env)

    x = m.addVars(cost_dict.keys(), obj=cost_dict, vtype=gp.GRB.BINARY, name='e')
    for i, j in x.keys():   
        x[j, i] = x[i, j]
    m.addConstrs(x.sum(i, '*') == 2 for i in range(n))
    m._vars = x
    m._n = n
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)
    tour_length = m.objVal
    for i in x:
        if x[i].x > 0.5:
            print(i)
    # 0-1-2-3-4-5-0  
    m.close()
    env.close()

    return tour_length

def NN(cost_table,csv_fname_out):
    TT = import_data_from_csv(cost_table)
    n = len(TT)
    sequence = np.zeros((n+1,),dtype=int)
    result = np.zeros(shape = (n,n+2),dtype=int)

    for w in range(n):
        sequence = np.zeros((n+1,),dtype=int)
        sequence[0] = w
        sequence[n] = w
        L = len(TT)
        k = 1

        TT2 = np.array(TT)
        TT2[:,sequence[0]] = np.nan

        for j in range(1,L):
            sequence[k] = np.nanargmin(TT2[sequence[k-1]])
            TT2[:,sequence[k]] = np.nan
            k += 1
        
        result[w,:-1] = sequence
        cost = 0
        len_sequence = len(sequence)
        for s in range(0,len_sequence-1):
            cost = np.nansum([cost,TT[sequence[s],sequence[s+1]]])
        result[w,n+1:n+2]  = cost
    
    if result[np.nanargmin([result[:,-1]]),:][0] == 0 :
        NN_best_route = result[np.nanargmin([result[:,-1]]),:-1]
        export_data_to_csv(NN_best_route,csv_fname_out)
    
    else:
        for i in range(n+1):
            if result[np.nanargmin([result[:,-1]]),:][i] == 0:
                NN_best_route = np.append(result[np.nanargmin([result[:,-1]]),i:-1],result[np.nanargmin([result[:,-1]]),:i+1])
                for j in range(n-1):
                    if NN_best_route[j] == NN_best_route[j+1]:
                        
                        NN_best_route = np.delete(NN_best_route,j)
                export_data_to_csv(NN_best_route,csv_fname_out)
                break
    return (result[np.nanargmin([result[:,-1]]),-1:])

def tasktime_table_maker(node_number, lower_number, upper_number, csv_fname):
    if node_number == 0 or lower_number == 0 or upper_number == 0:
        return

    data = np.zeros((node_number, node_number))
    for i in range(node_number):
        for j in range(node_number):
            if i == 0:
                continue
            elif j ==0:
                data[i][j] = np.nan
            else:
                data[i][j] = int(np.random.randint(lower_number, upper_number))

    export_data_to_csv(data, csv_fname)
    data = open(csv_fname, 'r').read()
    os.remove(csv_fname)
    archivo = open(csv_fname, 'w')
    archivo.write(data.replace('nan','').replace(' ',''))
    archivo.close()

def create_histogram(csv_fname):
    import matplotlib.pyplot as plt

    TT = import_data_from_csv(csv_fname)
    TT = np.triu(TT)
    TT = TT.flatten()
    TT = TT[~np.isnan(TT)]
    TT = TT.astype(int)
    plt.hist(TT[TT>0], bins=100)
    plt.show()

start = 100
end = 390
num_intervals = 10

interval_size = (end - start) / num_intervals

intervals = []
for i in range(num_intervals):
    interval_start = start + i * interval_size
    interval_end = interval_start + interval_size
    if i !=0:
        intervals.append((interval_start+1, interval_end))
    else:
        intervals.append((interval_start, interval_end))

path = 'Transitional/'

if not os.path.exists(path):
    os.makedirs(path)
else:
    folder_names = os.listdir(path)
    folder_names = [folder for folder in folder_names if '.txt' not in folder]
    for folder in folder_names:
        folder_path = os.path.join(path, folder)
        for file in os.listdir(folder_path):
            os.remove(os.path.join(folder_path, file))
        os.rmdir(folder_path)
    if os.path.exists(path+'tour_length.txt'):
        os.remove(path+'tour_length.txt')
os.makedirs(path+'best_routes')

tour_length_file = open(path+'tour_length.txt', 'w')

group = 'T' #Preffix of the new instance group
lower_dist = 50
upper_dist = 200

COT = "/Users/pgutiea/Desktop/TSPJ_MM/Codigos/Data/test/4_TT_paper.csv"
print(TSP(COT))
exit(0)


for instance in range(100):

    batch = (int(instance))//25+1 
    if not os.path.exists(path+'Batch_0'+str(batch)):
        os.makedirs(path+'Batch_0'+str(batch))

    NBC = f'{path}Batch_0{batch}/'+'TSPJ_'+str(instance+1)+str(group)+'_nodes_table_by_coordinates.csv'
    COT = f'{path}Batch_0{batch}/'+'TSPJ_'+str(instance+1)+str(group)+'_cost_table_by_coordinates.csv'  
    TAT = f'{path}Batch_0{batch}/'+'TSPJ_'+str(instance+1)+str(group)+'_tasktime_table.csv'
    NNB = f'{path}best_routes/'+'TSPJ_'+str(instance+1)+str(group)+'_NN_best_route.csv'


    L_node_number = int(intervals[instance//10][0])
    U_node_number = int(intervals[instance//10][1])

    node_number = random.randint(L_node_number, U_node_number)
    node_table_by_coordinates(node_number, lower_dist, upper_dist,NBC)
    cost_table_by_coordinates(NBC , COT)
    # tour_length = TSP(COT)
    tour_length = NN(COT,NNB)[0]
    tour_length_file.write(f'{instance+1}: {L_node_number}-{U_node_number}->{node_number}: {tour_length} \n')
    tasktime_table_maker(node_number , tour_length *.5 ,tour_length*.8 , TAT)
    print(f'{instance+1}: {L_node_number}-{U_node_number}->{node_number}: {tour_length}')
tour_length_file.close()

