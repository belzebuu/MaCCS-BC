import sys
from gurobipy import tuplelist
import random
from collections import OrderedDict
import networkx as nx
import copy
import pprint


class Data:
    def __init__(self, nodes=None, patients=None, mutations=None, interactions=None,
                 node_mutations=None, node_neighbors=None):
        """
        :param nodes: ([str])                   gene nodes
        :param patients: ([str])                patient ids
        :param mutations: ([(str, str)])        node-patient links
        :param interactions: ([(str, str)])     node-node links
        :param node_mutations: ({str: {str,}})  node to patient adjacency list
        :param node_neighbors: ({str: {str,}})  node adjacency list
        :return
        """
        self.gene_file = None
        self.interaction_file = None
        self.mutation_file = None
        self.weight_file = None
        
        self.nodes = list() if nodes is None else nodes
        self.patients = list() if patients is None else patients
        self.interactions = tuplelist([]) if interactions is None else interactions
        self.mutations = tuplelist([]) if mutations is None else mutations

        # auxiliary data structures
        self.node_mutations = OrderedDict() if node_mutations is None else node_mutations
        self.node_neighbors = OrderedDict() if node_neighbors is None else node_neighbors
        self.weights = dict()

    def __str__(self):
        m = "Read instance:\n"
        if self.gene_file:
            m += "Gene file: \'%s\'," % (self.gene_file)
        m += "Interaction file: \'%s\', Mutation file: \'%s\'\n" % (self.interaction_file, self.mutation_file)
        if self.weight_file:
            m += "Weight file: \'%s\'," % (self.weight_file)
        m += "Genes: %d, Interactions: %d, Patients: %d, Mutations: %d"\
               % (len(self.nodes), len(self.interactions), len(self.patients), len(self.mutations))
        return m

    # READ IN DATA
    def from_two_files(self, interaction_file, mutation_file, weight_file=None):
        self.interaction_file = interaction_file
        self.mutation_file = mutation_file

        used_genes = set()

        with open(interaction_file) as inf:
            for line in inf:
                node_a, node_b = line.split()
                used_genes.add(node_a)
                used_genes.add(node_b)
                # add interaction
                self.append_interaction(node_a, node_b)
        self.nodes = list(used_genes)

        with open(mutation_file) as mf:
            for line in mf:
                fields = line.split()
                patient = fields[0]
                if patient not in self.patients:
                    self.patients.append(patient)
                else:
                    print "Patient "+patient+" repeated."
                #patient_added = False
                for g_node in fields[1:]:
                    if g_node in used_genes:
                        # add patient
                        #if not patient_added:
                        #patient_added = True
                        self.append_mutation(patient, g_node)


        if weight_file:
            with open(weight_file) as wf:
                for line in wf:
                    g_node, weight = line.split()
                    if not float(weight) == 0.0:
                        self.weights[g_node] = float(weight)
        print self
        sys.exit(0)

    # READ IN DATA
    def from_three_files(self, gene_file, interaction_file, mutation_file, weight_file=None):
        """
        We know interactions among proteins. Proteins contain genes
        hence the interactions between genes can be deduced
        """
        self.gene_file = gene_file
        self.interaction_file = interaction_file
        self.mutation_file = mutation_file
        
        used_genes = set()
        used_nodes = set()

        gene_node = {}
        node_gene = {}

        with open(gene_file) as gf:
            for line in gf:
                n_str, g = line.split()
                # if gene_weights and g not in gene_weights:
                #     # print g
                #     continue
                n = n_str
                used_genes.add(g)
                used_nodes.add(n)
                if g in gene_node:
                    gene_node[g].append(n)
                else:
                    gene_node[g] = [n]
                    
                if n in node_gene:
                    node_gene[n].append(g)
                else:
                    node_gene[n] = [g]

        self.nodes = list(used_genes)

        pprint.pprint(node_gene)

        with open(interaction_file) as inf:
            for line in inf:
                a_str, b_str, c = line.split()
                a = a_str
                b = b_str
                if a in node_gene and b in node_gene:
                    # self.append_interaction(node_gene[a][0], node_gene[b][0])
                    for x in node_gene[a]:
                        for y in node_gene[b]:
                            self.append_interaction(x, y)



        # patients = list()
        # mutations = tuplelist([])
        # node_mut = dict()
        # gene_n_mut = dict()  - update
        with open(mutation_file) as mf:
            for line in mf:
                fields = line.split()
                #cancer_type = fields[0]
                patient = fields[0]
                if patient not in self.patients:
                    self.patients.append(patient)
                for gene in fields[1:]:
                    if gene in used_genes:
                        self.append_mutation(patient, gene)
                        

        if weight_file:
            max_weight = 0.0
            with open(weight_file) as wf:
                for line in wf:
                    g_node, weight = line.split()
                    if not float(weight) == 0.0:
                        self.weights[g_node] = float(weight)
                        if float(weight) > max_weight:
                            max_weight = float(weight)

            if max_weight > 1:
                print "Normalize weights from [0,%f] to [0,1]" % max_weight
                for gene in self.weights:
                    self.weights[gene] /= max_weight

        # TODO Remove genes (both from mutations and interaction graph) without a weight or with weight 0
        #if gene in gene_weights:
        #    if gene_weights[gene] == 0.0:
        #        continue
        #else:
        #    continue

                      
        print self
        #        print self.nodes
        #print self.patients
        sys.exit(0)



    # # # HELPER FUNCTIONS
    # done
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

    # done
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
