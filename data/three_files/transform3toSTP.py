#!/usr/bin/python
import pprint
import sys



class Problem:
        def __init__(self, gene_file, interaction_file, mutation_file):
            self.__dict__.update(
                nodes = self.read_genes(gene_file)
                edge_set = self.read_interactions(interaction_file)
                adj_node = self.read_mutations(mutation_file)
                )
            
        def read_genes(self):
		with open(inputfile) as gf:
		    for line in gf:
		        n_str, g = line.split()
		        # if gene_weights and g not in gene_weights:
		        #     # print g
		        #     continue
		        n = n_str
		        nodes.append(n)
		        genes.append(g)
		        if g in gene_node:
		            gene_node[g].append(n)
		        else:
		            gene_node[g] = [n]
		
		        if n in node_gene:
		            node_gene[n].append(g)
		        else:
		            node_gene[n] = [g]
		
		#gene_mut[g] = set()
		
        def read_interactions(self, inputfile):
		# interactions = tuplelist([])
		# node_neighbors = dict()
		edges = []
		edges_id=[]
		edge_set = set()
		with open(sys.argv[2]) as inf:
		    for line in inf:
		        aa, bb, c = line.split()
		        a = aa
		        b = bb
		        if a in node_gene and b in node_gene:
		            edges.append((node_gene[a][0],node_gene[b][0])) # genes
		            edge_set.add((node_gene[a][0],node_gene[b][0])) # genes)
		        if a in node_gene and b in node_gene:
		            edges_id.append((a,b)) # proteins
		
		
		#pp = pprint.PrettyPrinter(indent=4)
		##pp.pprint(gene_node)
		#pp.pprint(edges)
		#pp.pprint(node_gene)
		return edge_set


        def read_mutations(self, inputfile):
            adj_node={}
            with open(inputfile) as inf:
                for line in inf:
                    adj = line.split()
                    adj_node[adj[0]]=adj[1:len(adj)]
            return adj_node







def __main__():
    
    p = Problem(sys.argv[1], sys.argv[2], sys.argv[3])

    print(p)
    sys.exit(0)
    f= open("./interactions.txt","w")
    for e in edge_set:
        f.write(str(e[0])+" "+str(e[1])+"\n")
    f.close()
        
        
    f = open("./mutations.txt","w")
    for a in adj_node:
        print adj_node[a]
        f.write(a+" "+" ".join([str(x) for x in adj_node[a] if x in gene_node])+"\n")
        #map(lambda x: gene_node[x] if x in gene_node, adj[1:len(adj)])))
    f.close()
