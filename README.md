# Description

This repo contains source code for "Common NetworkS ReVealed" (CNSRV) presented in paper "Novel cross-species systems biology analyses reveal a deeply conserved transcriptional response to social challenge". Detailed explaination of the model please refer to the main paper and supplement [link].

# Code

### Compile:
1. cd makedir/
2. rm *
3. cmake ..
4. make

### Note:
1. You might need to update cmake version to the latest one
2. use the same steps to compile evaluate_CNSRV/

# Data:
See sample/data/

### Coexpression network:
1. Non-directed (or bi-directional) graph, if there is an edge between gene1 and gene2, then in the coexpression network file, there are two rows for this edge:

	gene1, gene2, 1 <br />
	gene2, gene1, 1

2. file format:
	1) 1st row of the file is species id
	2) 2nd row of the file is total number of genes (nodes) in the network of this species
	3) 3rd row to the end are coexpression edges

### Orthologous network:
1. Non-directed (or bi-directional) graph, if there is an edge between gene1 and gene2, then in the orthologous network file, there are two rows for this edge:

	gene1, gene2, 1 <br />
	gene2, gene1, 1

2. Orthologous edges are defined by orthologous mapping of genes between two species
3. file format:

	species1, species2, species1_gene_id, species2_gene_id

# To run:
$ ./makedir/CNSRV

### Run on sample dataset:
1. cd sample/run/
2. sh run.sh (please compiled CNSRV first)

### Evaluate results
1. cd sample/eval/
2. sh eval.sh (please compiled evaluate_CNSRV first)
3. Results of evaluation is a *.csv file
4. Script "./bin/combine_species_network.py" is used for combining coexpression networks of multiple species into one file
5. Metrics in evaluation are mainly In/Out density ratio of each cluster at both of coexpression and orthologous networks

### suggested parameters
1. Parameters used to generate clustering results in paper: 
	
	num_clusters: 20 <br />
	num_trials: 1 (run two independent trials and use the best one) <br />
	lambda: try range from 0 - 1, step size 0.1 <br />
	-t: 10

