# MaCCS-BC

Maximum Coverage Connected Subgraph by Branch and Cut for driver pathways analysis

## Installation

Download the package as zip or via git clone. Eg:

```
git clone https://github.com/belzebuu/MaCCS-BC.git
```

Currently, the package is not a Python module but a program to be used from command line.

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

To verify a solution:

``` bash
src/verifier.py -i data/interactions.txt -m data/mutations.txt -o sol.json
```
