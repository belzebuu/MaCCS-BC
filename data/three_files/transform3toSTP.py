#!/usr/bin/python3

import sys
from collections import defaultdict
from functools import reduce
import pprint


class Problem:
    def __init__(self, gene_file, interaction_file, mutation_file):
        self.gene2nodes, self.node2gene = self.read_mapping(gene_file)
        self.gene_edge_set, self.node_edge_set = self.read_interactions(interaction_file)
        self.mutations = self.read_mutations(mutation_file)
        self.patient_genes, self.node_adj = self.shrink_mutations_to_nodes()

    def read_mapping(self, inputfile):
        """
        return mappings:
        a gene can be mapped to more than one node
        a node can be mapped to only one gene.
        """
        gene_nodes = defaultdict(list)
        with open(inputfile) as gf:
            for line in gf:
                n, g = line.split()
                gene_nodes[g].append(int(n))
                
        node_gene = {v:g  for g in gene_nodes for v in gene_nodes[g] }

        node_set = set(reduce(lambda a,b: a+b, [gene_nodes[g] for g in gene_nodes]))
        assert(len(node_gene)==len(node_set))

        print("Genes read: {}, Nodes read: {}".format(len(gene_nodes),len(node_gene)))
        return gene_nodes, node_gene


    def read_interactions(self, inputfile):
        """
        return edges read for which a mapping exists.
        Repetitions are removed
        """
        node_edge_set = set()
        gene_edge_set = set()
        with open(sys.argv[2]) as inf:
            for line in inf:
                v1, v2, w = map(int,line.split())
                try:
                    if  v1 not in self.node2gene or v2 not in self.node2gene:
                        raise Exception("Node not in mapping")
                    if (v1,v2) not in node_edge_set:
                        node_edge_set.add((v1,v2))
                    else:
                        raise Exception("Edge already in node_edge_set: "+str(v1)+" "+str(v2))
                    if (self.node2gene[v1], self.node2gene[v2]) not in gene_edge_set:
                        gene_edge_set.add((self.node2gene[v1], self.node2gene[v2]))
                    else:
                        raise Exception("Edge already included: "+str(v1)+" "+str(v2))
                except Exception as ie:
                    continue
                    #print(str(ie))

        print("Gene inter. read: {}; Node edges read: {}".format(len(gene_edge_set),len(node_edge_set)))
        return gene_edge_set, node_edge_set


    def read_mutations(self, inputfile):
        """
        return dict of patients with mutated genes
        can contain genes not in mapping
        """
        patient_genes = {}
        genes = set()
        with open(inputfile) as inf:
            for line in inf:
                elements = line.split()
                genes.update(set(elements[1:]))
                if elements[0] not in patient_genes:
                    patient_genes[ elements[0] ] = [ g for g in elements[1:] ]
                else:
                    raise Exception("A repeated patient! Not handled")
        print("Patients read: {}, Genes read: {}".format(len(patient_genes),len(genes)))
        return patient_genes



    def shrink_mutations_to_nodes(self):
        """
        return mutations patients-genes in terms of nodes.
        """
        node_adj = {}
        patient_genes = {}

        counter = reduce(lambda x,y: max(x,y), [self.gene2nodes[g] for g in self.gene2nodes])[0]
        for m in self.mutations:
            genes = set(self.gene2nodes) & set(self.mutations[m]) # intersection     
            if len(genes)==0:
                continue
            counter = counter + 1
            node_adj[ counter ] = reduce(lambda x,y: x+y, [ self.gene2nodes[g] for g in genes ]) ##TODO: all or only one?
            patient_genes[ m ] = list(genes)

        print("After removing genes without mapping:")
        print("Patients read: {}, In nodes terms: {}".format(len(patient_genes),len(node_adj)))
        return patient_genes, node_adj

   
   




def writer_STP_format(prob):

    with open("./cmcp.stp","w") as f:
        f.write("33D32945 STP File, STP Format Version 1.0\n")
        f.write("SECTION Comment\n")
        f.write("Name    \"Connected Maximum Coverage Problem\"\n")
        f.write("Creator \"M. Chiarandini\"\n")
        f.write("Remark  \"Example of Connected Maximum Coverage Problem\"\n")   
        f.write("END\n\n")

        n_edges = len(prob.node_edge_set) + sum([len(prob.node_adj[m]) for m in prob.node_adj])

        f.write("SECTION Graph\n")
        f.write("Nodes {}\n".format(len(prob.node2gene)+len(prob.node_adj)))
        f.write("Edges {}\n".format(n_edges))
        count=0
        for e in prob.node_edge_set:
            f.write("E {} {}\n".format(str(e[0]),str(e[1])))               
            count+=1
        for t in prob.node_adj:
            for a in prob.node_adj[t]:
                f.write("E {} {}\n".format(str(t),str(a)))
                count+=1
        f.write("END\n\n")

        assert(count==n_edges)

        f.write("SECTION Terminals\n")
        f.write("Terminals {}\n".format(len(prob.node_adj)))
        for t in prob.node_adj:
            f.write("T {}\n".format(str(t)))
        f.write("END\n\n")


        f.write("SECTION CoverCardinality\n")
        f.write("Card {}\n".format(20))
        f.write("END\n\n")





def writer_two_files(prob):
    f= open("./interactions.txt","w")
    for e in prob.gene_edge_set:
        f.write(str(e[0])+" "+str(e[1])+"\n")
    f.close()

    f = open("./mutations.txt","w")
    for a in prob.patient_genes:
        f.write(a+" "+" ".join([str(x) for x in prob.patient_genes[a] if x in prob.gene2nodes])+"\n")
        #map(lambda x: gene2nodes[x] if x in gene2nodes, adj[1:len(adj)])))
    f.close()





def main(argv):
    print(argv)
    if len(argv) != 3:
        print("\nUsage: transform3toSTP.py gene_file interaction_file mutation_file\n")
        sys.exit(0)
    p = Problem(argv[0], argv[1], argv[2])

    #pp = pprint.PrettyPrinter(indent=4)
    #pp.pprint(p.gene2nodes)
    #pp.pprint(p.node_edge_set)
    #pp.pprint(p.patients)

    
    writer_STP_format(p)
    writer_two_files(p)



if __name__ == "__main__":
        main(sys.argv[1:])
