import sys
from gurobipy import tuplelist
import random
from collections import OrderedDict
import networkx as nx
import copy
import pprint


class Data:
    def __init__(self):
        #nodes=None, patients=None, mutations=None, interactions=None, node_mutations=None, node_neighbors=None):
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
        
        self.nodes = list() # if nodes is None else nodes
        self.patients = list() # if patients is None else patients
        #self.interactions = tuplelist([]) if interactions is None else interactions
        #self.mutations = tuplelist([]) if mutations is None else mutations

        # interactions are included twice, once in each direction
        self.interactions = list() # if interactions is None else interactions
        # mutations are included once, in the direction patient to gene
        self.mutations = list() # if mutations is None else mutations


        # auxiliary data structures
        self.node_neighbors = OrderedDict() # if node_neighbors is None else node_neighbors
        self.coverage = OrderedDict() # if node_mutations is None else node_mutations
        self.mutated_genes = OrderedDict() # if node_mutations is None else node_mutations

        self.weights = dict()

        self.components = list()



    def __str__(self):
        m = "Read instance:\n"
        if self.gene_file:
            m += "Gene file: \'%s\'," % (self.gene_file)
        m += "Interaction file: \'%s\', Mutation file: \'%s\'\n" % (self.interaction_file, self.mutation_file)
        if self.weight_file:
            m += "Weight file: \'%s\'," % (self.weight_file)
        m += "Genes: %d, Interactions: %d, Patients: %d,"\
               % (len(self.nodes), len(self.interactions), len(self.patients))
        m += " Mutations: %d, Mutated (given) genes: %d," % (len(self.mutations), len(self.coverage))
        if len(self.weights)>1:
            m += " Weights: %d" % ( len(self.weights))
        else:
            m+=" Unweighted"
        return m

    # READ IN DATA
    def from_two_files(self, interaction_file, mutation_file, weight_file=None):
        self.interaction_file = interaction_file
        self.mutation_file = mutation_file

        given_genes = set()

        with open(interaction_file) as inf:
            for line in inf:
                node_a, node_b = line.split()
                given_genes.add(node_a)
                given_genes.add(node_b)
                # add interaction
                self.append_interaction(node_a, node_b)
        self.nodes = list(given_genes)

        seen_genes = set()

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
                    seen_genes.add(g_node)
                    if g_node in given_genes:
                        # add patient
                        #if not patient_added:
                        #patient_added = True
                        self.append_mutation(patient, g_node)
            if len(self.patients) > len(self.mutated_genes.keys()):
                print "Discarded %d patients because they contain only mutations of non given genes" % (len(self.patients) - len(self.mutated_genes.keys()))
                del self.patients
                self.patients = self.mutated_genes.keys()
            if len(seen_genes) != len(given_genes):
                print "Genes in mutations %d and given %d" % (len(seen_genes), len(given_genes))

        if weight_file:
            with open(weight_file) as wf:
                for line in wf:
                    g_node, weight = line.split()
                    if not float(weight) == 0.0:
                        self.weights[g_node] = float(weight)
        print self


    # READ IN DATA
    def from_three_files(self, gene_file, interaction_file, mutation_file, weight_file=None):
        """
        We know interactions among proteins. Proteins contain genes
        hence the interactions between genes can be deduced
        """
        self.gene_file = gene_file
        self.interaction_file = interaction_file
        self.mutation_file = mutation_file
        
        given_genes = set()
        #given_labels = set()

        gene_node = {}
        label_gene = {}

        with open(gene_file) as gf:
            for line in gf:
                n_str, g = line.split()
                # if gene_weights and g not in gene_weights:
                #     # print g
                #     continue
                n = n_str
                given_genes.add(g)
                #given_labels.add(n)
                if g in gene_node:
                    gene_node[g].append(n)
                else:
                    gene_node[g] = [n]
                    
                if n in label_gene:
                    label_gene[n].append(g)
                else:
                    label_gene[n] = [g]

        self.nodes = list(given_genes)


        if False: # check relation:
            # Each protein has one gene only
            for a in label_gene:
                if len(label_gene[a])>1:
                    print(a,label_gene[a])

            # A gene can be in more than one protein
            for b in gene_node:
                if len(gene_node[b])>1:
                    print(b,gene_node[b])

        with open(interaction_file) as inf:
            for line in inf:
                a_str, b_str, c = line.split()
                a = a_str
                b = b_str
                if a in label_gene and b in label_gene:
                    # self.append_interaction(label_gene[a][0], label_gene[b][0])
                    for x in label_gene[a]:
                        for y in label_gene[b]:
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
                else:
                    print "Patient "+patient+" repeated"
                if len(fields)==1:
                    print patient
                for gene in fields[1:]:
                    if gene in given_genes:
                        self.append_mutation(patient, gene)
            if len(self.patients) > len(self.mutated_genes.keys()):
                print "Discarded %d patients because mutations of non given genes" % (len(self.patients) - len(self.mutated_genes.keys()))                
                del self.patients
                self.patients = self.mutated_genes.keys()

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



    # HELPER FUNCTIONS
    def append_interaction(self, node_a, node_b):
        #if (node_a, node_b) not in self.interactions:
        #self.interactions.append((node_a, node_b))
        #if (node_b, node_a) not in self.interactions:
        #self.interactions.append((node_b, node_a))

        if node_a in self.node_neighbors:
            if node_b not in self.node_neighbors[node_a]:
                self.interactions.append((node_a, node_b))
                self.node_neighbors[node_a].add(node_b)
        else:
            self.interactions.append((node_a, node_b))
            self.node_neighbors[node_a] = {node_b}
        if node_b in self.node_neighbors:
            if node_a not in self.node_neighbors[node_b]:
                self.interactions.append((node_b, node_a))
                self.node_neighbors[node_b].add(node_a)
        else:
            self.interactions.append((node_b, node_a))
            self.node_neighbors[node_b] = {node_a}
        #pprint.pprint(self.interactions)

    
    def append_mutation(self, patient, g_node):
        #if (patient, g_node) not in self.mutations:
        if g_node in self.coverage:
            if patient not in self.coverage[g_node]:
                self.mutations.append((patient, g_node))
                self.coverage[g_node].add(patient)
        else:
            self.mutations.append((patient, g_node))
            self.coverage[g_node] = {patient}

        if patient in self.mutated_genes:
            self.mutated_genes[patient].add(g_node)
        else:
            self.mutated_genes[patient] = {g_node}
            

    # SHUFFLE
    def permute_mutations(self, seed):
        print "shuffling keys"
        original = copy.deepcopy(self)

        index_dictionary = {n: n for n in self.nodes}
        keys = copy.deepcopy(self.nodes)
        random.seed(seed)
        random.shuffle(keys)
        shuffled_dictionary = dict(zip(keys, index_dictionary.values()))

        self.mutations = list()
        self.interactions = list()
        self.coverage = OrderedDict()
        self.mutated_genes = OrderedDict()
        self.node_neighbors = OrderedDict()

        for n in original.nodes:
            if n in original.coverage:
                for p in original.coverage[n]:
                    self.append_mutation(p, shuffled_dictionary[n])
            for u in original.node_neighbors[n]:
                self.append_interaction(shuffled_dictionary[n], shuffled_dictionary[u])


    def make_graph(self):
        G = nx.Graph()
        G.add_nodes_from(self.nodes)
        for [u, w] in self.interactions:
            G.add_edge(u, w)
        return G
