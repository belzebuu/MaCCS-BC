#!/usr/bin/python3

import sys
from collections import defaultdict
import pprint


class Problem:
    def __init__(self, gene_file, interaction_file, mutation_file):
        self.gene_node = self.read_genes(gene_file)
        self.edge_set = self.read_interactions(interaction_file)
        self.adj_node = self.read_mutations(mutation_file)


    def read_genes(self, inputfile):
        gene_node = {} #defaultdict(int)
        with open(inputfile) as gf:
            for line in gf:
                n, g = line.split()
                if g not in gene_node:
                    gene_node[g] = n
                else:
                    raise Exception("A gene represented more than once.")
        return gene_node


    def read_interactions(self, inputfile):
        node_set = {self.gene_node[g] for g in self.gene_node}
        print(node_set)
        edge_set = set()
        with open(sys.argv[2]) as inf:
            for line in inf:
                v1, v2, w = line.split()
                if v1 in node_set and v2 in node_set:
                    edge_set.add((v1,v2)) # genes)
        return edge_set


    def read_mutations(self, inputfile):
        patient_genes = {}
        patient_node = {}
        counter = reduce(lambda x,y: max(self.gene_node[x],self.gene_node[y]), self.gene_node)
        with open(inputfile) as inf:
            for line in inf:
                elements = line.split()
                counter = counter + 1
                if elements[0] not in patient_node:
                    patient_node[ elements[0] ] = counter
                else:
                    raise Exception("a repeated patient")
                try:
                    patient_genes[ counter ] = [ self.gene_node[g] for g in elements[1:] ]
                except IndexError as ie:
                    raise IndexError("a gene in mutations not listed in genes" + ie)
                except Exception as ie:
                    raise Exception("a gene in mutations not listed in genes" + ie)
        return patient_genes




def main(argv):
    print(argv)
    if len(argv) != 3:
        print("\nUsage: transform3toSTP.py gene_file interaction_file mutation_file\n")
        sys.exit(0)
    p = Problem(argv[0], argv[1], argv[2])

    pp = pprint.PrettyPrinter(indent=4)
    pp.pprint(p.gene_node)
    pp.pprint(p.edge_set)
    pp.pprint(p.patients)

    sys.exit(0)
    f= open("./interactions.txt","w")
    for e in edge_set:
        f.write(str(e[0])+" "+str(e[1])+"\n")
    f.close()


    f = open("./mutations.txt","w")
    for a in adj_node:
        print(adj_node[a])
        f.write(a+" "+" ".join([str(x) for x in adj_node[a] if x in gene_node])+"\n")
        #map(lambda x: gene_node[x] if x in gene_node, adj[1:len(adj)])))
    f.close()





if __name__ == "__main__":
        main(sys.argv[1:])
