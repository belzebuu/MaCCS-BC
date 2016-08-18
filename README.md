# MaCCS-BC

Maximum Coverage Connected Subgraph by Branch and Cut for driver pathways analysis

## Installation

Download the package as zip or via git clone. Eg:

```
git clone https://github.com/belzebuu/MaCCS-BC.git
```

Currently, MaCCS-BC consists of three Python scripts rather than being a library module. The three scripts find a subgraph, verify it and execute a permutation test, respectively.

## Usage

To find the set maximum coverage connected subgraph, ie, a driver pathway:

``` bash
cd maccs-bc
src/maccs-bc.py -i data/interactions.txt -m data/mutations.txt -k 10 -o sol.json
```

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
