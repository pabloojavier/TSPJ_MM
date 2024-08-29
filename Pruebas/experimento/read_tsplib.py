import tsplib95
from tsplib95.distances import euclidean, manhattan, geographical, pseudo_euclidean, xray, pseudo_euclidean,TYPES
import os

def read_problem(file):
    map_functions ={ 
        'ATT': pseudo_euclidean,
        'EUC_2D': euclidean,
        'EXPLICIT': lambda x,y: 0,
        'GEO': geographical,

    }

    problem = tsplib95.load(file)
    
    # Get the number of nodes
    n = len(list(problem.get_nodes()))
    nodes = list(problem.get_nodes())
    if len(problem.node_coords) != 0:
        coords = [problem.node_coords[i] for i in range(1,n+1)]
    else:
        coords = None
    edges = list(problem.get_edges())

    if problem.edge_weight_type in TYPES.keys():
        weights = {(edge[0], edge[1]): TYPES[problem.edge_weight_type](problem.node_coords[edge[0]], problem.node_coords[edge[1]]) for edge in edges}
    else:
        weights = {(edge[0], edge[1]): problem.get_weight(*edge) for edge in edges}

    # elif problem.edge_weight_type == 'EUC_2D':
    #     weights = {(edge[0], edge[1]): euclidean(problem.node_coords[edge[0]], problem.node_coords[edge[1]]) for edge in edges}

    # elif problem.edge_weight_type == 'EXPLICIT':
    #     weights = {e: problem.get_weight(*e) for e in edges}

    # elif problem.edge_weight_type == 'GEO':
    #     weights = {(edge[0], edge[1]): geographical(problem.node_coords[edge[0]], problem.node_coords[edge[1]]) for edge in edges}

    # else:
    #     raise ValueError('Edge weight type not supported')


    nodes = [node - 1 for node in nodes]
    edges = [(edge[0] - 1, edge[1] - 1) for edge in edges]
    weights = {(edge[0] - 1, edge[1] - 1): weights[edge] for edge in weights}

    return n, nodes, coords, edges, weights

def test():
    path = '/Users/pgutiea/Desktop/TSPJ_MM/Codigos/Data/tsp_experiment_gomhu'

    files = os.listdir(path)
    files.remove('.DS_Store')

    for file in files:
        try:
            n, nodes, coords, edges, weights = read_problem(path + '/' + file)
        except:
            print(file)
            raise
        # print(file,n,weights)

test()
for i in TYPES.keys():
    print(i, TYPES[i])