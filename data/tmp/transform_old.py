#!/usr/bin/python
import pprint
import sys

nodes = []
genes = []
gene_node = {}
node_gene = {}

with open(sys.argv[1]) as gf:
    for line in gf:
        n_str, g = line.split()
        # if gene_weights and g not in gene_weights:
        #     # print g
        #     continue
        n = int(n_str)
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


# interactions = tuplelist([])
# node_neighbors = dict()
edges = []
edges_id=[]
with open(sys.argv[2]) as inf:
    for line in inf:
        aa, bb, c = line.split()
        a = int(aa)
        b = int(bb)
        if a in node_gene and b in node_gene:
            edges.append((node_gene[a][0],node_gene[b][0]))
        if a in node_gene and b in node_gene:
            edges_id.append((a,b))


#pp = pprint.PrettyPrinter(indent=4)
##pp.pprint(gene_node)
#pp.pprint(edges)
#pp.pprint(node_gene)

#sys.exit(0)

f= open("./interaction.txt","w")
for e in edges_id:
    f.write(str(e[0])+" "+str(e[1])+"\n")
f.close()

edges = []
edges_id=[]
with open(sys.argv[2]) as inf:
    for line in inf:
        aa, bb, c = line.split()
        a = int(aa)
        b = int(bb)
        if a in node_gene and b in node_gene:
            edges.append((node_gene[a][0],node_gene[b][0]))
        if a in node_gene and b in node_gene:
            edges_id.append((a,b))




adj_node={}
with open(sys.argv[3]) as inf:
    for line in inf:
        adj = line.split()
        adj_node[adj[0]]=adj[1:len(adj)]

#print adj_node
        
        
f = open("mutations.txt","w")
for a in adj_node:
    print adj_node[a]
    f.write(a+" "+" ".join([str(x) for y in adj_node[a] if y in gene_node for x in gene_node[y] ])+"\n")
        #map(lambda x: gene_node[x] if x in gene_node, adj[1:len(adj)])))
f.close()
