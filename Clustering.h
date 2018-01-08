//
// Created by 杨巍 on 5/17/16.
//

#ifndef OCLODE_EFFICIENT_CLUSTER_H
#define OCLODE_EFFICIENT_CLUSTER_H


#include <unordered_set>
#include <unordered_map>
#include <map>
#include <sstream>
#include <math.h>

#include "SpeciesNetwork.h"
#include "Orthology.h"

#define MAXLOGSIZE 1000000
#define LINELEN 1024
#define HOWMANYITERTOSHOW 1000000
#define BADMOVETHRESHOLD 10000
#define ENDTEMPERATURE 100000
#define TIMESOFNODES 100


class Clustering
{
/* instances */
private:
    bool has_noise_cluster; // if clustering should include noise cluster
    int num_species;
    double lambda;
    int num_clusters;

    int *species_num_nodes; // number of nodes species co-express network
    int total_num_nodes;

    int undo_log_size;
    double current_coexpress_cost;
    double current_total_cost;
    double current_orth_cost;
    double ortho_normalize_constant;
    double coexpress_normalize_constant;

    int SA_counter;
    bool show_log;
    double max_temp;

    unordered_map<string, int> changed_nodes;
    vector< vector<int> > species_cluster_size;
    vector< vector<double> > species_cluster_inEdges;
    vector< vector<double> > species_cluster_outEdges;
    vector< vector<double> > species_cluster_coexpress_cost;

    int old_cluster_size_before_changed;
    double old_cluster_inEdges_before_changed;
    double old_cluster_outEdges_before_changed;
    double old_cluster_coexpress_cost_before_changed;

    int new_cluster_size_before_changed;
    double new_cluster_inEdges_before_changed;
    double new_cluster_outEdges_before_changed;
    double new_cluster_coexpress_cost_before_changed;

    double old_orth_cost;

private:
    struct LogEntry
    {
        int spc;
        int node;
        int oldstate;
    };

    LogEntry UndoLog[MAXLOGSIZE];

public:
    int **node_assigned_cluster;	// assignment of clusters for nodes in each species
    int **best_node_assigned_cluster; // best cluster assignment found
    vector< vector< unordered_set<int> > > spe_cluster_nodes; // 3D arr to store nodes per cluster per spe

    double best_total_cost;

/* methods */
public:

    Clustering(SpeciesNetwork **species_network,
               int num_species,
               Orthology *orthology,
               int num_clusters,
               double lambda,
               bool has_noise_cluster);

    ~Clustering();

    void learn_ground_state(SpeciesNetwork **species_network, Orthology *orthology);
    void pre_set(char *filename,Orthology *orthology, SpeciesNetwork **species_network);
    void set_max_temp(double mt);
    void Print(SpeciesNetwork **species_network);

private:
    void buffer_cluster_state(int size, double inEdges, double outEdges, double cost, bool is_new_state=false);
    void buffer_orth_cost(double cost);
    double get_orth_cost(Orthology *orthology);
    double get_coexpress_delta_cost(SpeciesNetwork **species_network);
    double get_coexpress_delta_cost_noNoiseCluster(SpeciesNetwork **species_network);
    double get_cluster_coexpress_delta_cost_noNoiseCluster(int spe, int cluster, int perturb_node, SpeciesNetwork **species_network, bool is_new_state=false);
    double get_spe_score(int spe, SpeciesNetwork **species_network);
    double get_coexpress_cost(SpeciesNetwork **species_network);
    double get_orth_delta_cost(Orthology *orthology, SpeciesNetwork **species_network, double curr_cost);

    int Perturb(SpeciesNetwork **species_network, int numchanges);
    void undo_perturb(int new_state);
    void update_best_clustering(double current_total_cost);
    vector<string> split(const string &s, char delim); // split string by delim into arr of strings

};


#endif //OCLODE_EFFICIENT_CLUSTER_H
