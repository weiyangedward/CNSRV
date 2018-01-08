# Compile:
1. cd makedir/
2. rm *
3. cmake ..
4. make

### Note:
1. You might need to update cmake version to the latest one

# To run:
$ ./makedir/CNSRV

### Run on sample dataset:
1. cd output/
2. sh run.sh

# TODO:

### Prepare dataset
1. species coexpression network
2. orthologous genes between two species

### Evaluate clustering results:
1. Compute coexpression/orthologous In/Out density ratio at each output cluster for each species