# MaCCS-BC

Maximum Coverage Connected Subgraph by Branch and Cut for driver pathways analysis

## Installation

Download the package as zip or via git clone. Eg:

```
git clone https://github.com/belzebuu/MaCCS-BC.git
```

Currently, MaCCS-BC consists of three Python scripts rather than being a
library module. The three scripts find a subgraph, verify it and execute
a permutation test, respectively.

## Usage

The data of the problem, that is, the genes, the interactions between
them, the mutations and the weights of the genes can be represented in
two alternative ways. In the extended way, the following set of files is
given:
1. a list of proteins with the corresponding genes carried by them
2. a list of knwon interactions between proteins.
3. a list of mutations related to a given type of cancer for each patient
4. a list of weights for each gene.  In the compressed format the first
two files are preprocessed and a single file representing the
interactions between genes is given.  See the documentation in the data
directory for more details on the input format.


Then, to find the maximum coverage connected subgraph, ie, a driver pathway
with extended data format:


``` bash
src/maccs-bc.py -g data/tmp/hippie_genes.txt -i data/tmp/hippie_edges.txt -m data/tmp/laml_pancancer.mm -k 10 -o sol.json
```
and with the compressed format:

``` bash
cd maccs-bc
src/maccs-bc.py -i data/interactions.txt -m data/mutations.txt -k 10 -o sol.json
```
Note that gene weights are optional. If a file containing a list of weights for each gene
is not given, then it is assumed that all genes have weight 1.


For a  list of options:

``` bash
src/maccs-bc.py -h
```

The set of genes composing the identified subgraph is written in the file sol.json.
To verify its properties:

``` bash
src/verifier.py -i data/interactions.txt -m data/mutations.txt -o sol.json
```

To execute a permutation test on the subgraph:

``` bash
src/permutation_test.py -i data/interactions.txt -m data/mutations.txt -k 10 -o sol.json -n 100
```
