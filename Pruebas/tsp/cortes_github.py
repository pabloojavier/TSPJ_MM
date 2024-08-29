from gurobipy import GRB
import gurobipy as gp
import networkx as nx


def CUT_integer_callback(m, where):
    
    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where == GRB.Callback.MIPSOL: 
        
        # retrieve the LP solution
        xval = m.cbGetSolution(m._x)
        
        # which edges are selected?
        edges_used = [ (i,j) for i,j in m._DG.edges if xval[i,j] > 0.5 ]
        
        # create support graph
        DG_soln = m._DG.edge_subgraph( edges_used )
        
        # if solution is not connected, add a (violated) CUT constraint for each subtour
        if not nx.is_strongly_connected( DG_soln ):
            for component in nx.strongly_connected_components( DG_soln ):
                complement = [ i for i in DG_soln.nodes if i not in component ]
                m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in complement ) >= 1 )

def CUT_naive_fractional_callback(m, where):
    
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
        xval = m.cbGetSolution(m._x)
    elif case2:
        xval = m.cbGetNodeRel(m._x)
        
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
            S = node_partition[0]  # 'left' side of the cut
            T = node_partition[1]  # 'right' side of the cut
            m.cbLazy( gp.quicksum( m._x[i,j] for i in S for j in T ) >= 1 )
            return
        
def CUT_smarter_fractional_callback(m, where):
    
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
        xval = m.cbGetSolution(m._x)
    elif case2:
        xval = m.cbGetNodeRel(m._x)
    
    DG = m._DG
    # if the support graph is disconnected, then finding violated cuts is easy!
    edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]
    DG_support = DG.edge_subgraph( edges_used )
    if not nx.is_strongly_connected( DG_support ):
        for component in nx.strongly_connected_components( DG_support ):
            complement = [ i for i in DG.nodes if i not in component ]
            m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in complement ) >= 1 )
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
                T = node_partition[1]  # 'right' side of the cut
                m.cbLazy( gp.quicksum( m._x[i,j] for i in S for j in T ) >= 1 )
                return
            
def DFJ_integer_callback(m, where):
    
    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where == GRB.Callback.MIPSOL: 
        DG = m._DG
        # retrieve the LP solution
        xval = m.cbGetSolution(m._x)
        
        # which edges are selected?
        edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > 0.5 ]
        
        # create support graph
        DG_soln = m._DG.edge_subgraph( edges_used )
        
        # if solution is not connected, add a (violated) DFJ constraint for each subtour
        if not nx.is_strongly_connected( DG_soln ):
            for component in nx.strongly_connected_components( DG_soln ):
                m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )

def DFJ_naive_fractional_callback(m, where):
    
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
        xval = m.cbGetSolution(m._x)
    elif case2:
        xval = m.cbGetNodeRel(m._x)
        
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
            S = node_partition[0]  # 'left' side of the cut
            m.cbLazy( gp.quicksum( m._x[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
            return
        
def DFJ_smarter_fractional_callback(m, where):
    
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
        xval = m.cbGetSolution(m._x)
    elif case2:
        xval = m.cbGetNodeRel(m._x)
    
    DG = m._DG
    # if the support graph is disconnected, then finding violated cuts is easy!
    edges_used = [ (i,j) for i,j in DG.edges if xval[i,j] > m._epsilon ]
    DG_support = DG.edge_subgraph( edges_used )
    if not nx.is_strongly_connected( DG_support ):
        for component in nx.strongly_connected_components( DG_support ):
            m.cbLazy( gp.quicksum( m._x[i,j] for i in component for j in component if i!=j ) <= len(component) - 1 )
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
                m.cbLazy( gp.quicksum( m._x[i,j] for i in S for j in S if i!=j ) <= len(S) - 1 )
                return
            
all_functions = [
    CUT_integer_callback, 
    CUT_naive_fractional_callback,
    CUT_smarter_fractional_callback, 
    DFJ_integer_callback, 
    DFJ_naive_fractional_callback,
    DFJ_smarter_fractional_callback
]