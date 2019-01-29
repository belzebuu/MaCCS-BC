from gurobipy import *
import common.graph as gh


# model functions
def cut(subset, edges):
    cut_set = list()
    for c in subset:
        for e in edges.select('*', c):
            if e[0] not in subset:
                cut_set.append(e[0])
        for e in edges.select(c, '*'):
            if e[1] not in subset:
                cut_set.append(e[1])
    return cut_set

# callback functions
def tolerant_int(x, int_feas_tol):
    '''
    :param x:               a list of values in the solution
    :param int_feas_tol:    feasibility tolerance, default:
    :return:                a list of binary values
    '''
    for i in range(len(x)):
        if x[i] >= (1-int_feas_tol):
            x[i] = 1.0
        elif x[i] <= int_feas_tol:
            x[i] = 0.0
        else:
            print("WRONG")
            print(i)
    return x


def make_int_check_connectivity(x_values, int_feas_tol, interactions, node_list):
    '''
    :param x_values:        list of values of x (counting from 0)
    :param int_feas_tol:    feasibility tolerance for getting a binary solution
    :param interactions:    edges between genes, given by data instance
    :param node_list:       translation from indices to nodes in data.nodes (node names) - for interctions
    :return:                list of 0-1 values, x_in_sol - indices of nodes in solution,
                            and True if solution connected, False otherwise
    '''
    in_sol = list()
    x_in_sol = list()
    edges = tuplelist([])
    for i in range(len(x_values)):
        if x_values[i] <= int_feas_tol:
            x_values[i] = 0
        elif x_values[i] >= (1-int_feas_tol):
            x_values[i] = 1
            node_name = node_list[i]
            if in_sol:
                for u in in_sol:
                    add_edges = tuplelist([])
                    # in interactions use name of node
                    if interactions.select('*', node_name):
                        if (u, node_name) in interactions.select('*', node_name):
                            add_edges.append((u, node_name))
                            # print "<%d, %d >" % (u, i+1)
                    if interactions.select(node_name, '*'):
                        if (node_name, u) in interactions.select(node_name, '*'):
                            add_edges.append((node_name, u))
                            # print "<%d, %d >" % (u, i+1)
                    if add_edges:
                        edges += add_edges
            # SOLUTION - indices, from 0
            in_sol.append(node_name)
            x_in_sol.append(i)
        else:
            print("WRONG VALUE OF x AT POSITION %d" % i)
    return x_values, x_in_sol, gh.bfs(in_sol, edges)


def tolerant_int_solution_keys(x, int_feas_tol, node_list):
    '''
    :param x:               a list of values in the solution
    :param int_feas_tol:    feasibility tolerance, default:
    :param node_list:       list of node keys
    :return:                a list of binary values
    '''
    solution_keys = set()
    for i in range(len(x)):
        if x[i] >= (1-int_feas_tol):
            solution_keys.add(node_list[i])
            x[i] = 1.0
        elif x[i] >= int_feas_tol:
            print("WRONG")
            print(i)
    return solution_keys


def check_connectivity(x_values, interactions):
    '''
    :param x_values:        list of values of x (counting from 0)
    :param interactions:    edges between genes, given by data instance
    :return:                is the solution connected or not
    '''
    in_sol = list()
    edges = tuplelist([])
    for i in range(0, len(x_values)):
        if x_values[i] == 1:
            if in_sol:
                for u in in_sol:
                    add_edges = tuplelist([])
                    if interactions.select('*', i+1):
                        if (u, i+1) in interactions.select('*', i+1):
                            add_edges.append((u, i+1))
                            # print "<%d, %d >" % (u, i+1)
                    if interactions.select(i+1, '*'):
                        if (i+1, u) in interactions.select(i+1, '*'):
                            add_edges.append((u, i+1))
                            # print "<%d, %d >" % (u, i+1)
                    if add_edges:
                        edges += add_edges
            in_sol.append(i+1)
    return gh.bfs(in_sol, edges)
