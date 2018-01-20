/*=======================================================================
 1. This code is used to test the noise cluster in a clustering results.
 Specifically, it aims to reason why the size of noise cluster is very
 small in the clustering result (noise size < 10).
 
 2. To justify this result, this code takes a node out of a cluster at
 the time, and see the change of final score, if the score decreases,
 then this node should not be put into the currently cluster and should
 be in the nosie cluster instead. If the clustering algorithm works
 correctly, this case should be very rare.
 
 3. Note that this code is not suitable for the new InOutRatio clustering
 algorithm, where the objective function is changed to:
 
 score =  sum_over_i(in_cluster_density(i) / ( in_cluster_density(i) + out_cluster_density(i))) * n + lamda * noise_term
 
 where lamda = sum_over_i(in_cluster_density(i) / ( in_cluster_density(i) + out_cluster_density(i))) / N
 
 ========================================================================*/

#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <list>
#include <unordered_set>
#include <stdio.h>
#include <string>

using namespace std;

#define EPSILON 0.0000001;

/* global variables */
int cluster_num = 0;
int species_num = 0;
int noise_cluster_flag = 1;
vector<vector<unordered_set<int> > > spe_clus; // arr to store genes per cluster per species
vector<int> total_gene_num(3, 0); // arr to store the total number of genes per species
string title; // "Best clustering..."
vector<unordered_map<int, vector<int> > >
        adj_list; // arr to store hash map for adjacency list per gene per spe for co-express edges
vector<vector<unordered_map<int, vector<int> > > >
        orth_adj_list; // adjacency list per gene per spe for orthologous edges

/* function declaration */
void in_out_cost(vector<vector<unordered_set<int> > > &spe_clus,
                 vector<unordered_map<int, vector<int> > > &adj_list,
                 vector<int> &total_gene_num,
                 char *output_file);

void read_clustering(char *clustering_file);

void read_network(char *network_file);


void in_out_cost(vector<vector<unordered_set<int> > > &spe_clus,
                 vector<unordered_map<int, vector<int> > > &adj_list,
                 vector<int> &total_gene_num,
                 char *output_file) {
    double ave_in_density, ave_out_density, ave_in_out_density_ratio;
    int count = 0;
    FILE *pFile2;
    pFile2 = fopen(output_file, "w");
    fprintf(pFile2, "Coexpression statistics:\n");
    fprintf(pFile2,
            "spe,cluster,#nodes,in-edge,out-edge,possible_in_edge,possible_out_edge,in-density,out-density,in/out-density_ratio,in/(in+out)_density_ratio\n");

    for (int i = 0; i < species_num; i++) {
        cout << "spe " << i << endl;
        for (int j = 0; j < cluster_num; j++) // 11 clusters
        {
            if (spe_clus[i][j].size() > 0) // if number of nodes > 0 in spe i cluster j
            {
                /*========== compute in/out density of cluster ============*/
                double in_edge = 0.0, out_edge = 0.0, num_node = (double) spe_clus[i][j].size();
                for (unordered_set<int>::iterator gene1 = spe_clus[i][j].begin();
                     gene1 != spe_clus[i][j].end(); gene1++) {
                    for (int g = 0; g < adj_list[i][*gene1].size(); g++) // for each gene connects to cur_gene
                    {
                        int gene2 = adj_list[i][*gene1][g];
                        if (spe_clus[i][j].find(gene2) != spe_clus[i][j].end()) { // if gene2 in cluster [i][j]
                            in_edge++;
                        } else {
                            out_edge++;
                        }
                    }
                }
                in_edge /= 2.0; // in_edge reduce to half because it was double counted in cluster
                double possible_in_edge = num_node * (num_node - 1) / 2.0;
                double possible_out_edge = num_node * (total_gene_num[i] - num_node);
                double in_density = in_edge / possible_in_edge;
                double out_density = out_edge / possible_out_edge;
                double in_out_density_ratio = in_density / out_density;
                double in_vs_inPlusOut_density_ratio = in_density / (in_density + out_density);
                double cost = (in_density / (in_density + out_density)) * num_node;

                fprintf(pFile2, "%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                        i, j, num_node,
                        in_edge, out_edge,
                        possible_in_edge, possible_out_edge,
                        in_density, out_density,
                        in_out_density_ratio, in_vs_inPlusOut_density_ratio);

                if (j != 0) // skip noise cluster
                {
                    count++;
                    ave_in_density += in_density;
                    ave_out_density += out_density;
                    ave_in_out_density_ratio += in_out_density_ratio;
                }
            }
        }
    }
    ave_in_density /= count;
    ave_out_density /= count;
    ave_in_out_density_ratio /= count;
//    fprintf(pFile2,
//            "\naverage_excluding_noise_cluster,ave_in_density,%g,ave_out_density,%g,ave_in_out_density_ratio,%g\n",
//            ave_in_density, ave_out_density, ave_in_out_density_ratio);


    /* orthologous in/out ratio */
    fprintf(pFile2, "Orthologous statistics:\n");
    fprintf(pFile2,
            "spe1,spe2,cluster,#nodes_in_spe1,#nodes_in_spe2,in-edge,out-edge1,out-edge2,possible_in_edge,possible_out_edge1,possible_out_edge2,in-density,out-density,in/out-density_ratio\n");
    for (int i = 0; i < species_num; i++) {
        for (int j = i + 1; j < species_num; j++) {
            for (int cluster = 0; cluster < cluster_num; cluster++) {
                double in_edge = 0.0, out_edge1 = 0.0, out_edge2 = 0.0,
                        num_node1 = (double) spe_clus[i][cluster].size(),
                        num_node2 = (double) spe_clus[j][cluster].size();
                if (spe_clus[i][cluster].size() > 0) {
                    for (unordered_set<int>::iterator gene1 = spe_clus[i][cluster].begin();
                         gene1 != spe_clus[i][cluster].end(); gene1++) {
                        fprintf(stderr, "spe %d, cluster %d, gene1 %d\n", i, cluster, *gene1);
                        for (int g = 0; g < orth_adj_list[i][j][*gene1].size(); g++) {
                            int gene2 = orth_adj_list[i][j][*gene1][g];
                            if (spe_clus[j][cluster].find(gene2) != spe_clus[j][cluster].end()) {
                                in_edge += 1.0;
                            } else {
                                out_edge1 += 1.0;
                            }
                        }
                    }

                    for (unordered_set<int>::iterator gene1 = spe_clus[j][cluster].begin();
                         gene1 != spe_clus[j][cluster].end(); gene1++) {
                        for (int g = 0; g < orth_adj_list[j][i][*gene1].size(); g++) {
                            int gene2 = orth_adj_list[j][i][*gene1][g];
                            if (spe_clus[i][cluster].find(gene2) == spe_clus[i][cluster].end()) {
                                out_edge2 += 1.0;
                            }
                        }
                    }
                    double possible_in_edge = num_node1 * num_node2;
                    double possible_out_edge1 = num_node1 * (total_gene_num[j] - num_node2);
                    double possible_out_edge2 = num_node2 * (total_gene_num[i] - num_node1);
                    double in_density = in_edge / possible_in_edge;
                    double out_density = ((out_edge1 / possible_out_edge1) + (out_edge2 / possible_out_edge2)) / 2;
                    if (out_density == 0)
                        out_density = EPSILON;
                    double in_out_density_ratio = in_density / out_density;
                    fprintf(pFile2, "%d,%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
                            i, j, cluster,
                            num_node1, num_node2,
                            in_edge, out_edge1, out_edge2,
                            possible_in_edge, possible_out_edge1, possible_out_edge2,
                            in_density, out_density, in_out_density_ratio);
                }
            }
        }
    }

    fclose(pFile2);
}

/*
 * store genes per species per cluster in spe_clus
 * */
void read_clustering(char *clustering_file) {
    ifstream CLUSTEROUT(clustering_file);
    if (!CLUSTEROUT) {
        cerr << "clustering output could not be opened for reading!" << endl;
        exit(1);
    }

    getline(CLUSTEROUT, title);

    // initialize spe_clus
    for (int i = 0; i < 3; i++) {
        vector<unordered_set<int> > spe;
        for (int j = 0; j < cluster_num; j++) {
            unordered_set<int> clus;
            spe.push_back(clus);
        }
        spe_clus.push_back(spe);
    }

    // read file CLUSTEROUT
    while (CLUSTEROUT) {
        string line = "";
        getline(CLUSTEROUT, line);
        if (line.compare("") != 0) {
            //        fprintf(stderr, "%s\n", line.c_str());
            istringstream split_line(line);
            string w1, w2, w3;
            int spe_id, gene_id, cluster_id;
            split_line >> w1 >> spe_id >> w2 >> gene_id >> w3 >> cluster_id;

            total_gene_num[spe_id]++;

            spe_clus[spe_id][cluster_id].insert(gene_id);
        }
    }
    CLUSTEROUT.close();
    printf("total number of genes: hb:%d, mm:%d, sb:%d\n", total_gene_num[0], total_gene_num[1], total_gene_num[2]);
}

void read_network(char *network_file) {
    ifstream NETWORK(network_file);
    if (!NETWORK) {
        cerr << "network file could not be opened for reading!" << endl;
        exit(1);
    }

    getline(NETWORK, title); // "Real network"

    for (int i = 0; i < species_num; i++) {
        unordered_map<int, vector<int> > spe_adj_list;
        adj_list.push_back(spe_adj_list);
    }

    // read file NETWORK
    while (NETWORK) {
        string line = "";
        getline(NETWORK, line);
        if (line.compare("") != 0) {
            istringstream split_line(line);
            int spe, gene_id1, gene_id2, edge;
            split_line >> spe >> gene_id1 >> gene_id2 >> edge;

            /* only add genes for id1 since edges are double counted */
            if (adj_list[spe].find(gene_id1) == adj_list[spe].end()) // if key does not exist in hash table
            {
                vector<int> tmp;
                tmp.push_back(gene_id2);
                adj_list[spe][gene_id1] = tmp;
            } else { // if key exists
                adj_list[spe].at(gene_id1).push_back(gene_id2);
            }
        }
    }
    NETWORK.close();
}

void read_orthologous(char *ortho_file) {
    ifstream ORTHO(ortho_file);
    if (!ORTHO) {
        cerr << "network file could not be opened for reading!" << endl;
        exit(1);
    }

    for (int i = 0; i < species_num; i++) {
        vector<unordered_map<int, vector<int> > > spe1_adj_list;
        for (int j = 0; j < species_num; j++) {
            unordered_map<int, vector<int> > spe2_adj_list;
            spe1_adj_list.push_back(spe2_adj_list);
        }
        orth_adj_list.push_back(spe1_adj_list);
    }

    while (ORTHO) {
        string line = "";
        getline(ORTHO, line);
        if (line.compare("") != 0) {
            istringstream split_line(line);
            int spe1, spe2, gene_id1, gene_id2;
            split_line >> spe1 >> spe2 >> gene_id1 >> gene_id2;

            if (orth_adj_list[spe1][spe2].find(gene_id1) ==
                orth_adj_list[spe1][spe2].end()) // if key does not exist in hash table
            {
                vector<int> tmp;
                tmp.push_back(gene_id2);
                orth_adj_list[spe1][spe2][gene_id1] = tmp;
            } else {// if key exists
                orth_adj_list[spe1][spe2][gene_id1].push_back(gene_id2);
            }
        }
    }
    ORTHO.close();
}


int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("USAGE:\n"
                       "\t%s clustering_output network_file orth_file eval_output_filename(.csv) cluster_num species_num"
                       "\nWhere:\n"
                       "\tclustering_output: clustering output from CNSRV\n"
                       "\tnetwork_file: coexpression edges\n"
                       "\torth_file: orthologous edges between two species\n"
                       "\teval_output_filename(.csv): filename of evaluation output\n"
                       "\tcluster_num: number of clusters\n"
                       "\tspecies_num: number of species\n"

                , argv[0]);
        exit(1);
    }

    cluster_num = stoi(argv[5]);
    species_num = stoi(argv[6]);
    // noise_cluster_flag = stoi(argv[7]);
    if (noise_cluster_flag == 1) {
        cluster_num += 1;
    } // include noise cluster index

    /* read clustering output file */
    fprintf(stderr, "read clustering output file\n");
    read_clustering(argv[1]);

    /* read network file */
    fprintf(stderr, "read network file\n");
    read_network(argv[2]);

    fprintf(stderr, "read orthologous file\n");
    read_orthologous(argv[3]);


    /* compute in-out cost per cluster */
    fprintf(stderr, "compute in-out cost per cluster\n");
    in_out_cost(spe_clus, adj_list, total_gene_num, argv[4]);

    return 0;
}
