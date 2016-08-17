from itertools import combinations


class NodeSep:
    def __init__(self, i, j, separator):
        self.node_i = i
        self.node_j = j
        self.separator = separator

    def __str__(self):
        return "[" + str(self.node_i) + ", " + str(self.node_j) + " ]\n" + str(self.separator)


def unvisited_node(component_idx):
    for k in component_idx.keys():
        if component_idx[k] == 0:
            return k
    return None


def bfs_k(k, x_sol_keys, node_neighbors):
    black = set([k])
    seen = set()
    # accessible from k
    gray = x_sol_keys & node_neighbors[k]
    # find all the neighbors
    while gray:
        for x in gray:
            seen |= node_neighbors[x] & x_sol_keys
        black |= gray
        gray = seen - black
    return black


def solution_components(x_sol_keys, node_neighbors):
    '''
    :param x_sol_keys:          indices of nodes in solution
    :param node_neighbors:      dictionary with neighbors set per each node
    :return:                    dictionary with sets of connected nodes
    '''
    components = dict()
    component_idx = {node: 0 for node in x_sol_keys}
    while True:
        k = unvisited_node(component_idx)
        if k:
            if k not in node_neighbors:
                print "%d has no neighbors" % k
                component_idx[k] = k
                continue
            accessible_k = bfs_k(k, x_sol_keys, node_neighbors)
            components[k] = accessible_k
            for n in accessible_k:
                component_idx[n] = k
        else:
            break
    return components


def component_neighbors(component, node_neighbors):
    '''
    :param component:           set of connected nodes
    :param node_neighbors:      dictionary with neighbors set per each node
    :return:                    set of neighboring nodes to the component
    '''
    neighbors = set()
    for i in component:
        neighbors |= node_neighbors[i]
    neighbors -= component
    return neighbors


def node_separators(component_i, component_j, node_neighbors):
    '''
    :param component_i:          set of connected nodes
    :param component_j:          set of connected nodes, disjoint with component_i
    :param node_neighbors:       dictionary with neighbors set per each node
    :return:
    '''
    A_i = component_neighbors(component_i, node_neighbors)
    separator = set()

    black = set()
    old_gray = component_j
    while len(old_gray) > 0:
        gray = set()
        for i in old_gray:
            # nodes of A_Ci never put in BFS queue
            to_add = node_neighbors[i] - A_i
            separator |= node_neighbors[i] & A_i
            gray |= to_add
        black |= old_gray
        gray = gray - black
        old_gray = gray
    return separator


def violated_separators(x_sol_keys, node_neighbors, add_all):
    components = solution_components(x_sol_keys, node_neighbors)
    if len(components) == 1:
        return None
    n_current_violations = 0
    violations = list()
    components_combinations = combinations(sorted(components.keys()), 2)
    lower_bnd = 0
    for cc in components_combinations:
        if cc[0] < cc[1]:
            c_i = components[cc[0]]
            c_j = components[cc[1]]
        else:
            c_i = components[cc[1]]
            c_j = components[cc[0]]
        lower_bnd += len(c_i) * len(c_j)

        separator_Ci_Cj = node_separators(c_i, c_j, node_neighbors)
        # print separator_Ci_Cj
        c_i_sorted = sorted(c_i)
        c_j_sorted = sorted(c_j)
        for i in c_i_sorted:
            for j in c_j_sorted:
                violations.append(NodeSep(i, j, separator_Ci_Cj))
                n_current_violations += 1
        if add_all:
            separator_Cj_Ci = node_separators(c_j, c_i, node_neighbors)
            if separator_Cj_Ci != separator_Ci_Cj:
                for j in c_j_sorted:
                    for i in c_i_sorted:
                        violations.append(NodeSep(j, i, separator_Cj_Ci))
                        n_current_violations += 1

    # assert lower_bnd <= len(violations)
    # assert len(violations) <= 2*lower_bnd
    # print "viols lower bound: %d <= number viols: %d <= upper bound: %d" %(lower_bnd, len(violations), lower_bnd*2)

    return violations
