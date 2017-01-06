import sys
from gurobipy import tuplelist
import random
from collections import OrderedDict
import networkx as nx
import copy
import data as d
from time import time 

class Preprocess:
    def __init__(self):
        return

    # done
    def preprocess_with_weights(self, data, k):
        data_components = list()
        small_components = set()
        g = self.make_graph(data)
        components = nx.connected_components(g)
        for c in components:

            if len(c) < k:
                small_components |= {c}
                continue

            print "connected component"
            nodes_to_remove = set()
            nodes_degree_one = set()
            nodes_covering_less = set()
            nodes_noncovering = set()

            for i in c:
                if i in small_components:
                    continue

                if i in data.node_mutations:
                    for u in c:
                        if u in data.node_mutations and i != u:
                            if data.node_neighbors[i] <= data.node_neighbors[u]\
                                    and data.node_mutations[i] <= data.node_mutations[u]:
                                # 4: node u strictly "better"
                                if data.weights[i] <= data.weights[u]:
                                    nodes_covering_less.add(i)
                                    break
                else:
                    # 1: nodes of degree 1 not covering any patients
                    if len(data.node_neighbors[i]) == 1:
                        nodes_degree_one.add(i)
                    else:
                        for u in c:
                            if data.node_neighbors[i] <= data.node_neighbors[u]:
                                if u in data.node_mutations:
                                    # 3: more useful (covering) node that also connects nodes
                                    nodes_noncovering.add(i)
                                    break
                                elif data.node_neighbors[i] != data.node_neighbors[u] or u < i:
                                    # 2: another, smaller index node that has same neighbors, and covers nothing
                                    nodes_noncovering.add(i)
                                    break

            print "before, in component : %d " % len(c)
            print "degree 1, no cover: %d " % len(nodes_degree_one)
            print "exists better node: %d " % len(nodes_covering_less)
            print "no cover, same nbs: %d " % len(nodes_noncovering)
            print "small components:   %d " % len(small_components)
            nodes_to_remove |= small_components
            nodes_to_remove |= nodes_degree_one
            nodes_to_remove |= nodes_covering_less
            nodes_to_remove |= nodes_noncovering
            print "all nodes to remove: %d " % len(nodes_to_remove)
            if len(c) - len(nodes_to_remove) >= k:
                # # nodes_to_keep = [x for x in c if x not in nodes_to_remove]
                # print "nodes after prep:   %d " % len(nodes_to_keep)
                data_components.append(data.data_from_nodelist(c, nodes_to_remove, weights=data.weights))
        return data_components

    # # # PREPROCESS
    # done
    def make_graph(self, data):
        graph = nx.Graph()
        graph.add_nodes_from(data.nodes)
        for [u, w] in data.interactions:
            graph.add_edge(u, w)
        return graph

    # done
    def preprocess(self, data, k):
        if data.weights:
            return data.preprocess_with_weights(k)
        data_components = list()
        small_components = set()
        g = self.make_graph(data)
        components = nx.connected_components(g)

	start_time= time()
        for c in components:

            if len(c) < k:
                small_components |= set(c)
                continue

            nodes_to_remove = set()
            nodes_degree_one = set()
            nodes_covering_less = set()
            nodes_noncovering = set()

	    # TEST
	    # data_components.append(self.data_from_nodelist(data, c, {'8152'}))
	    # return data_components
	    # TEST
            for i in c:
                if i in data.node_mutations:
                    for u in c:
                        if u in data.node_mutations and i != u:
                            if data.node_neighbors[i] <= data.node_neighbors[u]\
                                    and data.node_mutations[i] <= data.node_mutations[u]:
                                # 4: node u strictly "better"
                                if data.node_neighbors[i] != data.node_neighbors[u] \
                                        or data.node_mutations[i] != data.node_mutations[u]:
                                    nodes_covering_less.add(i)
                                    break
                                # 5: u at least as good and with smaller index
                                elif u < i:
                                    nodes_covering_less.add(i)
                                    break
                else:
                    # 1: nodes of degree 1 not covering any patients
                    if len(data.node_neighbors[i]) == 1:
                        nodes_degree_one.add(i)
                    else:
                        for u in c:
                            if data.node_neighbors[i] <= data.node_neighbors[u]:
                                if u in data.node_mutations:
                                    # 3: more useful (covering) node that also connects nodes
                                    nodes_noncovering.add(i)
                                    break
                                elif data.node_neighbors[i] != data.node_neighbors[u] or u < i:
                                    # 2: another, smaller index node that has same neighbors, and covers nothing
                                    nodes_noncovering.add(i)
                                    break

            print "before, in component : %d " % len(c)
            print "degree 1, no cover: %d " % len(nodes_degree_one)
            print "exists better node: %d " % len(nodes_covering_less)
            print "no cover, same nbs: %d " % len(nodes_noncovering)
            print "small components:   %d " % len(small_components)
            nodes_to_remove |= small_components
            nodes_to_remove |= nodes_degree_one
            nodes_to_remove |= nodes_covering_less
            nodes_to_remove |= nodes_noncovering
            print "all nodes to remove: %d " % len(nodes_to_remove)
            print nodes_covering_less 
            if len(c) - len(nodes_to_remove) >= k:
                # nodes_to_keep = [x for x in c if x not in nodes_to_remove]
                # print "nodes after prep:   %d " % len(nodes_to_keep)
                data_components.append(self.data_from_nodelist(data, c, nodes_to_remove))
        print 'Execution finished in %.3f seconds' % (time() - start_time)	
        return data_components

    # done
    def data_from_nodelist(self, data, all_nodes, nodes_to_remove=set(), weights=None):
        if not nodes_to_remove:
            return self
        new_data = d.Data()
        print new_data
        # sys.exit(0)
        component_patient_set = set()
        for i in all_nodes:
            if i in nodes_to_remove:
                continue
            
            new_data.nodes.append(i)

            # update interactions
            if i in data.node_neighbors:
                for u in data.node_neighbors[i]:
                    if u not in nodes_to_remove:
                        # save interaction only once
                        if i < u:
                            new_data.interactions.append((i, u))
                        # update auxiliary adjacency list data structure
                        if i in new_data.node_neighbors:
                            new_data.node_neighbors[i].add(u)
                        else:
                            new_data.node_neighbors[i] = {u}
                if weights:
                    new_data.weights[i] = weights[i]

            # update mutations
            if i in data.node_mutations:
                new_data.node_mutations[i] = data.node_mutations[i]
                for p in new_data.node_mutations[i]:
                    new_data.mutations.append((p, i))
                    component_patient_set.add(p)

        # print new_data.nodes
        # for i in new_data.node_mutations.keys():
        #     print i,i in new_data.nodes
        new_data.patients = list(component_patient_set)
        return new_data

    # # TODO ?
    # def preprocess_genes(self, k):
    #     preprocessed = dict()
    #
    #     connected_components = list()
    #     data_components = list()
    #     small_components = set()
    #     g = self.make_graph(data)
    #     components = nx.connected_components(g)
    #     for c in components:
    #         if len(c) < k:
    #             small_components |= set(c)
    #         else:
    #             connected_components.append(c)
    #     for c in connected_components:
    #         nodes_to_remove = set()
    #         nodes_degree_one = set()
    #         nodes_covering_less = set()
    #         nodes_noncovering = set()
    #
    #         for i in c:
    #             if i in small_components:
    #                 continue
    #             if i in data.node_mut:
    #                 for u in c:
    #                     if u in data.node_mut and i != u:
    #                         if data.node_neighbors[i] <= data.node_neighbors[u]\
    #                                 and data.node_mut[i] <= data.node_mut[u]:
    #                             # 4: node u strictly "better"
    #                             if data.node_neighbors[i] != data.node_neighbors[u] \
    #                                     or data.node_mut[i] != data.node_mut[u]:
    #                                 nodes_covering_less.add(i)
    #                                 # if u in preprocessed:
    #                                 #     preprocessed[u].add(i)
    #                                 # else:
    #                                 #     preprocessed[i] = {i}
    #                                 # break
    #                             # 5: u at least as good and with smaller index
    #                             elif u < i:
    #                                 nodes_covering_less.add(i)
    #                                 if u in preprocessed:
    #                                     preprocessed[u].add(i)
    #                                 else:
    #                                     preprocessed[u] = {i}
    #                                 break
    #             else:
    #                 # 1: nodes of degree 1 not covering any patients
    #                 if len(data.node_neighbors[i]) == 1:
    #                     nodes_degree_one.add(i)
    #                 else:
    #                     for u in c:
    #                         if data.node_neighbors[i] <= data.node_neighbors[u]:
    #                             if u in data.node_mut:
    #                                 # 3: more useful (covering) node that also connects nodes
    #                                 nodes_noncovering.add(i)
    #                                 break
    #                             elif data.node_neighbors[i] != data.node_neighbors[u] or u < i:
    #                                 # 2: another, smaller index node that has same neighbors, and covers nothing
    #                                 nodes_noncovering.add(i)
    #                                 break
    #
    #         print "before, in component : %d " % len(c)
    #         print "degree 1, no cover: %d " % len(nodes_degree_one)
    #         print "exists better node: %d " % len(nodes_covering_less)
    #         print "no cover, same nbs: %d " % len(nodes_noncovering)
    #         print "small components:   %d " % len(small_components)
    #     return preprocessed
