- get back to the old instance format (three files). [Anna]

- make sure nothing in the instance format is as interpreted as integer
  numbers. [Anna]

- update all the documentation that we
  wrote. [Marco after Anna has done the previous two points]

- resolve incongruencies in preprocessing output: [Marco] eg:

```bash
Instance: 14930 genes, 194 patients, 1413 mutations, 179806 interactions
Preprocessing instance...
before, in component : 14897
degree 1, no cover: 2069
exists better node: 21
no cover, same nbs: 1874
small components:   0
all nodes to remove: 3964
Instance: 0 genes, 0 patients, 0 mutations, 0 interactions
single component
```





- Add more experiments on weighted model
  - Fabio provides new files with data
  - experiments similar to Tab 4 of WABI paper <- seek script
  
  
- decide k
```bash
  k-fold crossvalidation to decide k - Unweighted case 
  
  Assume 3000 patients. Use them for: Training, evaluation, assessment

  Divide data in 10 blocks of 300 patients

  - Training and evaluation by 10-fold cross validation

	For each block b:
		take block b out from data
		For k from 5 to 30 by intervals of h:
			solution[k,-b] = solve on data[-k] to optimality for k
			calculate coverage[k,b] on block b of solution[k,-b] 
	
	For k from 5 to 30 by intervals of h:
		coverage[k] = median over the blocks of coverage[k,b]
	
	k_best = arg min { coverage[k] }
		
  - Assessment
  For the best k determined above:
  Repeat 1000 times:
	  sample 90% of data
	  solve to optimality for k_best
  histogram of genes present in the solutions (frequency of each gene selected)

  Peeking avoided since we do not look at genes in selecting k...
```
- Add flow formulation
- Add better description heuristics

- Alternative objective function
- Jacquard Similarity analysis
- Graph of interdependencies
