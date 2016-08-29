import sys
from gurobipy import tuplelist
import random
from collections import OrderedDict
import networkx as nx
import copy


class Data:
    def __init__(self, nodes=None, patients=None, mutations=None, interactions=None,
                 node_mutations=None, node_neighbors=None, gene_nodes=None):
        """
        :param nodes: ([str])                   gene nodes
        :param patients: ([str])                patient ids
        :param mutations: ([(str, str)])        node-patient links
        :param interactions: ([(str, str)])     node-node links
        :param node_mutations: ({str: {str,}})  node to patient adjacency list
        :param node_neighbors: ({str: {str,}})  node adjacency list
        :return
        """
        self.nodes = list() if nodes is None else nodes
        self.patients = list() if patients is None else patients
        self.interactions = tuplelist([]) if interactions is None else interactions
        self.mutations = tuplelist([]) if mutations is None else mutations

        # auxiliary data structures
        self.node_mutations = OrderedDict() if node_mutations is None else node_mutations
        self.node_neighbors = OrderedDict() if node_neighbors is None else node_neighbors
        self.gene_nodes = OrderedDict() if gene_nodes is None else gene_nodes
        self.weights = dict()

    def __str__(self):
        return "Instance: %d genes, %d patients, %d mutations, %d interactions"\
               % (len(self.nodes), len(self.patients), len(self.mutations), len(self.interactions))

    # # # READ IN DATA
    @classmethod
    def from_file(cls, interaction_file, mutation_file, gene_file, weight_file=None):
        print
        data = cls()
        used_nodes = set()

        with open(gene_file) as gf:
            for line in gf:
                n_str, g = line.split()
                n = int(n_str)
                if g in data.gene_nodes:
                    data.gene_nodes[g].append(n)
                else:
                    data.gene_nodes[g] = [n]

        with open(interaction_file) as inf:
            for line in inf:
                node_a, node_b = line.split()
                if node_a in data.gene_nodes.values() and node_b in data.gene_nodes.values():
                    used_nodes.add(node_a)
                    used_nodes.add(node_b)
                    # add interaction
                    data.append_interaction(node_a, node_b)
        data.nodes = list(used_nodes)

        with open(mutation_file) as mf:
            for line in mf:
                fields = line.split()
                # # ctype = fields[0]
                patient = fields[0]
                for g in fields[1:]:
                    for n in data.gene_nodes[g]:
                        if n in used_nodes:
                            data.patients.append(patient)
                            data.append_mutation(patient, n)

        if weight_file:
            with open(weight_file) as wf:
                for line in wf:
                    gene, weight = line.split()
                    if gene in data.gene_nodes and float(weight) != 0.0:
                        for node in data.gene_nodes[gene]:
                            if node in used_nodes:
                                data.weights[node] = float(weight)
        print data
        return data

    # # # HELPER FUNCTIONS
    def append_interaction(self, node_a, node_b):
        self.interactions.append((node_a, node_b))
        if node_a in self.node_neighbors:
            self.node_neighbors[node_a].add(node_b)
        else:
            self.node_neighbors[node_a] = {node_b}
        if node_b in self.node_neighbors:
            self.node_neighbors[node_b].add(node_a)
        else:
            self.node_neighbors[node_b] = {node_a}

    def append_mutation(self, patient, g_node):
        self.mutations.append((patient, g_node))
        if g_node in self.node_mutations:
            self.node_mutations[g_node].add(patient)
        else:
            self.node_mutations[g_node] = {patient}

    # # # SHUFFLE
    def permute_mutations(self, seed):
        print "shuffling keys"
        original = copy.deepcopy(self)

        index_dictionary = {n: n for n in self.nodes}
        keys = copy.deepcopy(self.nodes)
        random.seed(seed)
        random.shuffle(keys)
        shuffled_dictionary = dict(zip(keys, index_dictionary.values()))

        self.mutations = tuplelist()
        self.interactions = tuplelist()
        self.node_mutations = OrderedDict()
        self.node_neighbors = OrderedDict()

        for n in original.nodes:
            if n in original.node_mutations:
                for p in original.node_mutations[n]:
                    self.append_mutation(p, shuffled_dictionary[n])
            for u in original.node_neighbors[n]:
                self.append_interaction(shuffled_dictionary[n], shuffled_dictionary[u])

    def make_graph(self):
        G = nx.Graph()
        G.add_nodes_from(self.nodes)
        for [u, w] in self.interactions:
            G.add_edge(u, w)
        return G
