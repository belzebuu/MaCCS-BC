from gurobipy import *


def connectivity_bfs_adj_list(v, edges, k):
    """
    :param v:       first of the k nodes
    :param edges:   dictionary of edges of a approximation induced by the set of k nodes
    :param k:       number of nodes that should be connected
    :return:        is the solution connected or not
    """
    black = set()
    old_gray = set()
    old_gray.add(v)
    while len(old_gray) > 0:
        gray = set()
        for i in old_gray:
            for w in edges[i]:
                gray.add(w)
        black |= old_gray
        gray = gray - black
        old_gray = gray
    if len(black) == k:
        return True
    return False


def paths_bfs_adj_list(v, edges, k):
    """
    :param v:       first of the k nodes
    :param edges:   dictionary of edges of a approximation induced by the set of k nodes
    :param k:       max depth of bfs tree
    :return:        dictionary of shortest paths from v
    """
    black = set()
    old_gray = set()
    old_gray.add(v)
    paths = dict()
    paths[v] = list()
    paths[v].append(v)
    degree = 0
    while degree < k:
        gray = set()
        for i in old_gray:
            if i in edges:
                for w in edges[i]:
                    gray.add(w)
                    if w not in paths:
                        paths[w] = paths[i] + [w]
        degree += 1
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return paths


def best_paths_bfs_adj_list(v, edges, k, node_mut):
    """
    :param v:       first of the k nodes
    :param edges:   dictionary of edges of a approximation induced by the set of k nodes
    :param k:       max depth of bfs tree
    :return:        dictionary of shortest paths from v
    """
    black = set()
    old_gray = set()
    old_gray.add(v)
    paths = dict()
    paths[v] = [v]
    paths_coverage = dict()
    paths_coverage[v] = len(coverage({v}, node_mut))
    degree = 0
    while degree < k:
        gray = set()
        for i in old_gray:
            if i in edges:
                for w in edges[i]:
                    gray.add(w)
                    new_path = paths[i] + [w]
                    if w not in paths:
                        paths[w] = new_path
                        paths_coverage[w] = len(coverage(set(new_path), node_mut))
                    else:
                        current_coverage = paths_coverage[w]
                        new_path_coverage = coverage(set(new_path), node_mut)
                        if len(new_path_coverage) > current_coverage:
                            paths[w] = new_path
                            paths_coverage[w] = len(coverage(new_path, node_mut))
        degree += 1
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return paths


def paths_bfs_adj_list_from_set(black, old_gray, paths, edges, k):
    """
    :param black:       seen, not to be looked at again
    :param old_gray:    seen
    :param paths:       already found paths
    :param edges:       dictionary of edges of a approximation induced by the set of k nodes
    :param k:           max depth of bfs tree
    :return:            dictionary of shortest paths from v
    """
    #print old_gray
    #print black
    old_gray -= black
    gray = old_gray
    for v in old_gray:
        if v not in paths:
            paths[v] = list()
            paths[v].append(v)
    degree = 0
    while degree < k:
        gray = set()
        for i in old_gray:
            if i in edges:
                for w in edges[i]:
                    gray.add(w)
                    if w not in paths:
                        paths[w] = paths[i] + [w]
        degree += 1
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return paths, black, gray


def coverage(a_set, node_mut):
    cov = set()
    if a_set:
        for v in a_set:
            if v in node_mut:
                cov |= node_mut[v]
    return cov


def paths_bfs_adj_list_ratio(v, edges, k, node_mut):
    """
    :param v:       first of the k nodes
    :param edges:   dictionary of edges of a approximation induced by the set of k nodes
    :param k:       max depth of bfs tree
    :return:        dictionary of shortest paths from v, shortest is the one maximizing ratio pv(u)/lv(u)
    """
    black = set()
    old_gray = set()
    old_gray.add(v)
    paths = dict()
    ratios = dict()
    paths[v] = list()
    paths[v].append(v)
    ratios[v] = len(coverage(set([v]), node_mut))
    degree = 0
    while degree < k:
        gray = set()
        for i in old_gray:
            for w in edges[i]:
                gray.add(w)
                if w not in paths:
                    p_v_u = paths[i] + [w]
                    sc_v_u_elements = set(p_v_u)  # = gh.SetCov(set(p_v_u), node_mut)
                    sc_v_u_coverage = coverage(sc_v_u_elements, node_mut)
                    paths[w] = p_v_u
                    ratios[w] = 1.0 * len(sc_v_u_coverage) / len(sc_v_u_elements)
                else:
                    p_v_u = paths[i] + [w]
                    #sc_v_u = gh.SetCov(set(p_v_u), node_mut)
                    sc_v_u_elements = set(p_v_u)
                    sc_v_u_coverage = coverage(sc_v_u_elements, node_mut)
                    ratio = 1.0 * len(sc_v_u_coverage) / len(sc_v_u_elements)
                    if ratio > ratios[w]:
                        # update to better path
                        paths[w] = paths[i] + [w]
                        ratios[w] = ratio
        degree += 1
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return paths


# ## finds shortest paths from v up to depth k using bfs and tuplelists
def k_bfs_shortest_paths(v, k, interactions):
    """
    :param v:               node to find shortest paths from
    :param k:               max depth of bfs tree
    :param interactions:    tuplelist of edges between nodes in Vg
    :return:                dictionary with paths (lists of nodes)
    """
    black = set()
    old_gray = set()
    old_gray.add(v)
    paths = dict()
    paths[v] = list()
    paths[v].append(v)
    degree = 0
    while degree < k:
        gray = set()
        for i in old_gray:
            for (u, w) in interactions.select(i, '*'):
                gray.add(w)
                if w not in paths:
                    paths[w] = paths[i] + [w]
                    #print len(paths[w])
            for (w, u) in interactions.select('*', i):
                gray.add(w)
                if w not in paths:
                    paths[w] = paths[i] + [w]
                    #print len(paths[w])
        degree += 1
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return paths

# TODO update to use adjacency list and not tuplelists
# ## checking solution found in FORBIDDEN_SOL for connectivity using bfs and tuplelists
def bfs(nodes, edges):
    """
    bfs check if all nodes connected
    :param nodes:   set or list - of k node indices (counting from 1)
    :param edges:   tuplelist - edges of a approximation induced by the set of k nodes
    :return:        is the solution connected or not
    """
    black = set()
    old_gray = set()
    old_gray.add(nodes[0])
    while len(old_gray) > 0:
        gray = set()
        for i in old_gray:
            for (u, w) in edges.select(i, '*'):
                gray.add(w)
            for (w, u) in edges.select('*', i):
                gray.add(w)
        black |= old_gray
        gray = gray - black
        old_gray = gray
    if len(black) == len(nodes):
        return True
    return False

# TODO update to use adjacency list and not tuplelists
def make_int_check_connectivity(x_values, int_feas_tol, interactions):
    """
    :param x_values:        list of values of x (counting from 0)
    :param int_feas_tol:    feasibility tolerance for getting a binary solution
    :param interactions:    edges between genes, given by data instance
    :return:                list of 0-1 values, in_sol - indices of nodes in solution,
                            and True if solution connected, False otherwise
    """
    in_sol = list()
    edges = tuplelist([])
    for i in range(0, len(x_values)):
        if x_values[i] <= int_feas_tol:
            x_values[i] = 0
        elif x_values[i] >= (1-int_feas_tol):
            x_values[i] = 1
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
        else:
            print("WRONG VALUE OF x AT POSITION %d" % (i+1))
    return x_values, in_sol, bfs(in_sol, edges)


def check_connectivity(x_values, interactions):
    """
    :param x_values:        list of values of x (counting from 0)
    :param interactions:    edges between genes, given by data instance
    :return:                is the solution connected or not
    """
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
    return bfs(in_sol, edges)
