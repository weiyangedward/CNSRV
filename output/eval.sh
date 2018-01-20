python ../bin/combine_species_network.py ../data/network0.txt ../data/network1.txt ../data/network2.txt > combined_networks

../evaluate_CNSRV/makedir/evaluate_CNSRV cluster.out combined_networks ../data/orth.txt eval_out.csv 10 3
