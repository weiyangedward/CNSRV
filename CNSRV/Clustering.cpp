#include "Clustering.h"

/**
 * Destructor
 */
Clustering::~Clustering()
{
    fprintf(stderr, "Clustering destructor\n");
    for (int spe=0; spe<num_species; spe++) {
        delete node_assigned_cluster[spe];
        delete best_node_assigned_cluster[spe];
        species_cluster_size.clear();
        species_cluster_inEdges.clear();
        species_cluster_outEdges.clear();
        species_cluster_coexpress_cost.clear();
        for (int clus = 0; clus< num_clusters; clus++) {
            spe_cluster_nodes[spe][clus].clear();
        }
        spe_cluster_nodes[spe].clear();
    }
    delete[] species_num_nodes;
    delete[] node_assigned_cluster;
    delete[] best_node_assigned_cluster;
    spe_cluster_nodes.clear();
    changed_nodes.clear();
}

/**
 * Constructor
 * @param species_network
 * @param num_species
 * @param orthology
 * @param num_clusters
 * @param lambda
 * @param has_noise_cluster
 */
Clustering::Clustering(SpeciesNetwork **species_network, int num_species, Orthology *orthology, int num_clusters, double lambda, bool has_noise_cluster)
        : num_species(num_species), num_clusters(num_clusters), lambda(lambda), has_noise_cluster(has_noise_cluster) {
    srand((int)time(NULL)); // random seed for local time
//    srand(1); // fix random seed

    show_log = false; // change manually


    SA_counter = 0;
    undo_log_size = 0;

    /* init total number of nodes in a species network */
    species_num_nodes = new int[num_species];
    for (int spe=0; spe<num_species; spe++)
    {
        species_num_nodes[spe] = species_network[spe]->num_nodes;
        total_num_nodes += species_num_nodes[spe];
    }

    /* init cluster size, inDens, outDens */
    species_cluster_size = vector< vector<int> >(num_species,vector<int>(num_clusters,0));
    species_cluster_inEdges = vector< vector<double> >(num_species,vector<double>(num_clusters,0));
    species_cluster_outEdges = vector< vector<double> >(num_species,vector<double>(num_clusters,0));
    species_cluster_coexpress_cost = vector< vector<double> >(num_species,vector<double>(num_clusters,0));
    /* init 3D arr to store nodes per cluster per species */
    for (int spe = 0; spe<num_species; spe++)
    {
        vector< unordered_set<int> > tmp_cluster_node;
        for (int clus = 0; clus< num_clusters; clus++)
        {
            unordered_set<int> tmp_node;
            tmp_cluster_node.push_back(tmp_node);
        }
        spe_cluster_nodes.push_back(tmp_cluster_node);
    }

    /* randomly assign cluster to nodes */
    node_assigned_cluster = new int *[num_species];
    best_node_assigned_cluster = new int *[num_species];
    for (int spe=0; spe<num_species; spe++)
    {
        node_assigned_cluster[spe] = new int[species_num_nodes[spe]];
        best_node_assigned_cluster[spe] = new int[species_num_nodes[spe]];
        for (int node=0; node<species_num_nodes[spe]; node++)
        {
            // randomly assign clusters to each node
            int rand_index = 0;
            // noise cluster included, cluster (0,k) will be assigned
            if (has_noise_cluster)
                rand_index = rand() % num_clusters;
                // no noise cluster -> no cluster0, only cluster(1,k) will be assigned
            else
                rand_index = (rand() % (num_clusters-1))+1;

            node_assigned_cluster[spe][node] = rand_index; // [spe][node] = cluster
            spe_cluster_nodes[spe][rand_index].insert(node); // [spe][cluster].emplace(node,1), maybe to keep this only (and remove node_assigned_cluster) is enough???
        }
    }

    // score the init cluster
    current_coexpress_cost = 0.0;
    current_orth_cost = 0.0;

    /**
     * a normalize constant to scale ortho cost to
     * the same as coexpress cost = #nodes in all species
     * given cost for each cluster is [0,1]
     * */
//    ortho_normalize_constant = (double)num_species / ((double)num_species * ((double)num_species-1) / 2.0);
    ortho_normalize_constant = (double)num_species / ((double)num_species * ((double)num_species-1) / 2.0);

    coexpress_normalize_constant = (double)total_num_nodes / ((double)num_clusters * (double)num_species / 2);

    if (lambda > 0 && lambda < 1) /* has orthologous term */
    {
        current_coexpress_cost = get_coexpress_cost(species_network);
        current_orth_cost = get_orth_cost(orthology);
        current_total_cost = (1-lambda) * current_coexpress_cost * coexpress_normalize_constant + lambda * current_orth_cost * ortho_normalize_constant;
    }
    else if (lambda == 0)
    { /* no orthologous term */
        current_coexpress_cost = get_coexpress_cost(species_network);
        current_total_cost = current_coexpress_cost * coexpress_normalize_constant;
    }
    else if (lambda == 1)
    {
        current_orth_cost = get_orth_cost(orthology);
        current_total_cost = current_orth_cost * ortho_normalize_constant;
    }
    fprintf(stderr, "orth_cost %g, ortho_normalize_constant %g, coexpress_normalize_constant %g, coexpression cc = %g, current_coexpress_cost %g\n", current_orth_cost, ortho_normalize_constant, coexpress_normalize_constant, 1-lambda, current_coexpress_cost);

    // init best_total_cost to curr cost
    best_total_cost = current_total_cost;
}


/* print final clustering result to output */
void Clustering::Print(SpeciesNetwork **species_network)
{
    for (int spe=0; spe<num_species; spe++)
    {
        for (unordered_map<string, int>::iterator x= species_network[spe]->geneName_to_uniqId.begin(); x!=species_network[spe]->geneName_to_uniqId.end(); x++)
        {
            string name = x->first;
            int id = x->second;
//            int id = species_network[i]->geneName_to_uniqId[name];
            printf("Species %d\tGene %s\tCluster %d\n",spe,name.c_str(),best_node_assigned_cluster[spe][id]);
        }
    }
}

/* set max tempreture for simulated annealing */
void Clustering::set_max_temp(double mt) { max_temp = mt; }

/* use grouth truth to init cluster labels */
void Clustering::pre_set(char *filename, Orthology *orthology, SpeciesNetwork **species_network)
{
    FILE *F = fopen(filename, "r");
    if (F==NULL) return;

    char line[LINELEN];
    fgets(line, LINELEN, F); // skip the single header line
    while (fgets(line, LINELEN, F))
    {
        char spc[LINELEN];
        char id[LINELEN];
        char clusterid[LINELEN];
        char d1[LINELEN], d2[LINELEN], d3[LINELEN];
        // format of line is "Species <spc> Gene <gene> Cluster <cluster>"
        sscanf(line, "%s %s %s %s %s %s",(char *)&d1, (char *)&spc, (char *)&d2, (char *)&id, (char *)&d3, (char *)&clusterid);

        // convert spc id to numerics
        if (orthology->speciesName_to_id.find(spc) == orthology->speciesName_to_id.end()) {
            printf("Error: seed-clustering file %s has species name that wasnt seen in network files\n", filename);
            exit(1);
        }
        int spcid = orthology->speciesName_to_id[spc];

        // convert gene id to numerics
        if (species_network[spcid]->geneName_to_uniqId.find(id) == species_network[spcid]->geneName_to_uniqId.end()) {
            continue; // this may happen because clusters file includes information about genes that are not there in
            // ... the network file (gene node with degree 0). It doesnt matter where (which cluster) we place such genes.
        }
        int geneid = species_network[spcid]->geneName_to_uniqId[id];

        // record information
        node_assigned_cluster[spcid][geneid] = atoi(clusterid);
        fprintf(stderr,"Initial clustering has species %d gene %s assigned to cluster %d\n", spcid, id, node_assigned_cluster[spcid][geneid]);
    }
    fclose(F);
    current_coexpress_cost = get_coexpress_cost(species_network);
    current_total_cost = current_coexpress_cost + lambda * get_orth_cost(orthology);
}


/*==============================================================================
 simulate annealing to find the best clustering

 1) for tempreture from max to min (max/1000),
    with step size = 0.9
        for 50*size_of_network
            compute delta_cost = new_cost - old_cost
            if delta_cost >= 0, accept

            new move with probability = exp(-delta_cost/temp)
            if delta_cost < 0, accept new move
    break for loop if reject new move for 500 times
 2) update best clustering when current clustering is better with a margine of 1
 ==============================================================================*/
void Clustering::learn_ground_state(SpeciesNetwork **species_network, Orthology *orthology)
{
    fprintf(stderr,"Called LearnGroundState\n");
    double tempmax = max_temp; // set start tempreture
    double tempmin = tempmax/ENDTEMPERATURE; // set end tempreture

    // count total number of nodes in all species
    int num_total_nodes = 0;
    for (int i=0; i<num_species; i++)
        num_total_nodes += species_num_nodes[i];

    // counter for number of iterations that no better score is seen
    int cost_not_change = 0;

    /* simulated annealing (SA), each step = 0.9,
     loop until tempreture is 1000 times smaller */
    for(double temp=tempmax; temp >= tempmin; temp *= 0.9)
    {
        if (cost_not_change > BADMOVETHRESHOLD) break;

        // a "sweep" of the network does about 20 x as many changes as there are nodes overall
        for (int i=0; i < TIMESOFNODES*num_total_nodes; i++)
        {
            // counter for SA iterations, used to determine when to print log
            SA_counter++;
            /* stop for loop if more than 500 iteration without change of cost,
             here the counter 'cost_not_change' is the same as outer loop,
             which means once the inner loop break, the outer will also break */
            if (cost_not_change > BADMOVETHRESHOLD) break;

            // record the old cost
            double old_total_cost = current_total_cost;
            double old_coexpress_cost = current_coexpress_cost;

            int new_state = Perturb(species_network, 1);

            // compute new cost of coexpress network
            double coexpress_delta_cost = 0.0;
            double new_total_cost = 0.0;
            double new_coexpress_cost = 0.0;
            double delta_orth_cost = 0.0;

//            double new_coexpress_cost = get_coexpress_cost(species_network);

            /* ====================================
             * compute new cost adding ortho term
             * only compute get_orth_cost or
             * get_coexpress_delta_cost when needed
             * ====================================*/

            /*
             * delta cost is used
             * */
            if (lambda > 0 && lambda < 1)
            {
                if (has_noise_cluster) {
                    coexpress_delta_cost = get_coexpress_delta_cost(species_network);
                } else {
                    coexpress_delta_cost = get_coexpress_delta_cost_noNoiseCluster(species_network);
                }

                new_coexpress_cost = old_coexpress_cost + coexpress_delta_cost;
                delta_orth_cost = get_orth_delta_cost(orthology, species_network, current_orth_cost);

                current_orth_cost += delta_orth_cost;

                new_total_cost = (1-lambda) * new_coexpress_cost * coexpress_normalize_constant + lambda * current_orth_cost * ortho_normalize_constant;

            }
            else if (lambda == 0)
            {
                if (has_noise_cluster) {
                    coexpress_delta_cost = get_coexpress_delta_cost(species_network);
                } else {
                    coexpress_delta_cost = get_coexpress_delta_cost_noNoiseCluster(species_network);
                }
                new_coexpress_cost = old_coexpress_cost + coexpress_delta_cost;
                new_total_cost = new_coexpress_cost * coexpress_normalize_constant;
            }
            else if (lambda == 1)
            {
                delta_orth_cost = get_orth_delta_cost(orthology, species_network, current_orth_cost);

                current_orth_cost += delta_orth_cost;
                new_total_cost = current_orth_cost * ortho_normalize_constant;
            }

            /*
             * not using delta cost
             * */
//            if (lambda > 0 && lambda < 1)
//            {
//                new_coexpress_cost = get_coexpress_cost(species_network);
//
//                delta_orth_cost = get_orth_delta_cost(orthology, species_network, current_orth_cost);
//
//                current_orth_cost += delta_orth_cost;
//
//                new_total_cost = (1-lambda) * new_coexpress_cost + lambda * current_orth_cost * ortho_normalize_constant;
//
//            }
//            else if (lambda == 0)
//            {
//                new_coexpress_cost = get_coexpress_cost(species_network);
//                new_total_cost = new_coexpress_cost;
//            }
//            else if (lambda == 1)
//            {
//                delta_orth_cost = get_orth_delta_cost(orthology, species_network, current_orth_cost);
//
//                current_orth_cost += delta_orth_cost;
//                new_total_cost = current_orth_cost * ortho_normalize_constant;
//            }



            if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
                fprintf(stderr, "orth_cost %g, ortho_normalize_constant %g, coexpress_normalize_constant %g, coexpress cost %g, lambda %g\n", current_orth_cost, ortho_normalize_constant, coexpress_normalize_constant, new_coexpress_cost, lambda);

            double delta_total_cost = new_total_cost - old_total_cost; // compute delta cost that includes orth term

//            if (SA_counter % HOWMANYITERTOSHOW == 0)
//                fprintf(stderr, "delta cost %g, old cost %g\n", delta_total_cost, old_total_cost);

            /* bad move, new_total_cost >= old_total_cost */
            if (delta_total_cost >= 0)
            {
                /* reject bad move
                 when a random double (0~1) >= e^(-delta / t)
                 if new cost is more, i.e., worse,
                 delta is more positive, thus exp(-delta/T) is closer to 0,
                 thus reject probability is larger */
                if (double(rand())/RAND_MAX >= exp(-delta_total_cost/temp))
                {
                    // prev log
                    if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
                        fprintf(stderr, "REJECTED MOVE\t%g to %g at t %g, prob %g\n", old_total_cost, new_total_cost, temp, 1-exp(-delta_total_cost/temp));

                    // undo pertubation
                    undo_perturb(new_state);
                    // increment counter no change of cost
                    cost_not_change++;
                }
                    /* accept bad move
                     with probability = exp(-delta_total_cost/temp) */
                else
                {
                    current_total_cost = new_total_cost;
                    current_coexpress_cost = new_coexpress_cost;

                    if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
                        fprintf(stderr, "BAD MOVE\t%g to %g at t %g, prob %g\n", old_total_cost, new_total_cost, temp, exp(-delta_total_cost/temp)); // prev log

                    // reset counter for no change of cost
                    cost_not_change = 0;
                    // reset counter if accept move
                    undo_log_size = 0;
                }
            }
                /* good moves, accept it.
                 when delta_total_cost < 0, if new_total_cost < old_total_cost */
            else
            {
                current_total_cost = new_total_cost;
                current_coexpress_cost = new_coexpress_cost;

                if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
                    fprintf(stderr, "GOOD MOVE\t%g to %g at t %g\n", old_total_cost, new_total_cost, temp); // prev log, print out good move

                /* any improvement less than this is not counted as an improvement,
                 here set to 0, all improvements are counted */
                if (delta_total_cost < -0.01)
                    cost_not_change = 0;

                // reset log size
                undo_log_size = 0;

                /* update best_total_cost and clustering
                 * curr cost < best_total_cost - 1,
                 * here 1 is a margine of improvement */
                if (current_total_cost < best_total_cost - 0.1)
                    update_best_clustering(current_total_cost);

            }
            // print size of noise cluster
            if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
            {
                for (int spe = 0; spe<num_species; spe++)
                {
                    fprintf(stderr, "spe %d, noise size %lu\n", spe, spe_cluster_nodes[spe][0].size());
                }
                fprintf(stderr, "C = %g\n", current_total_cost); // prev log, print out curr cost
            }
        }
    }
}

void Clustering::update_best_clustering(double current_total_cost)
{
    // update best_total_cost
    best_total_cost = current_total_cost;
    // update best clustering

    for (unordered_map<string, int>::iterator x=changed_nodes.begin(); x!=changed_nodes.end(); x++)
    {
        int new_state = x->second;
        vector<string> spe_node = split(x->first, '.');
        best_node_assigned_cluster[stoi(spe_node[0])][stoi(spe_node[1])] = new_state;
    }
    changed_nodes.clear();
}

/*============================================
 * compute total delta cost (no noise clsuter)
 * ===========================================*/
double Clustering::get_coexpress_delta_cost_noNoiseCluster(SpeciesNetwork **species_network)
{
    double delta_cost = 0.0;

    int perturb_node = UndoLog[undo_log_size-1].node;
    int spe = UndoLog[undo_log_size-1].spc;
    int old_cluster = UndoLog[undo_log_size-1].oldstate;
    int new_cluster = node_assigned_cluster[spe][perturb_node];

    double old_cluster_deltaCost = get_cluster_coexpress_delta_cost_noNoiseCluster(spe, old_cluster, perturb_node, species_network, false);

    spe_cluster_nodes[spe][old_cluster].erase(perturb_node);
    spe_cluster_nodes[spe][new_cluster].insert(perturb_node);

    double new_cluster_deltaCost = get_cluster_coexpress_delta_cost_noNoiseCluster(spe, new_cluster, perturb_node, species_network, true);
    delta_cost = old_cluster_deltaCost + new_cluster_deltaCost;
//    fprintf(stderr, "delta_cost %g\n", delta_cost);
    return delta_cost;
}

/*======================================================
 * compute coexpress delta cost for the old/new cluster
 * with respect to perturbed node (no noise cluster)
 * ====================================================*/
double Clustering::get_cluster_coexpress_delta_cost_noNoiseCluster(int spe, int cluster, int perturb_node, SpeciesNetwork **species_network, bool is_new_state)
{
    double cluster_delta_score = 0.0, in_edge_changed = 0.0, out_edge_changed = 0.0, after_changed_cluster_size = 0.0;

    /* buffer cluster state for recover if undo_perturb */
    if (!is_new_state)
    {
        buffer_cluster_state(species_cluster_size[spe][cluster],
                             species_cluster_inEdges[spe][cluster],
                             species_cluster_outEdges[spe][cluster],
                             species_cluster_coexpress_cost[spe][cluster], false);
        after_changed_cluster_size = species_cluster_size[spe][cluster] - 1;
    } else {
        buffer_cluster_state(species_cluster_size[spe][cluster],
                             species_cluster_inEdges[spe][cluster],
                             species_cluster_outEdges[spe][cluster],
                             species_cluster_coexpress_cost[spe][cluster], true);
        after_changed_cluster_size = species_cluster_size[spe][cluster] + 1;
    }

    if (after_changed_cluster_size < 0)
    {
        fprintf(stderr, "after_changed_cluster_size < 0, spe %d, cluster %d\n", spe, cluster);
        exit(1);
    }

    /* count in, out edges */
    double possible_in_edge = after_changed_cluster_size * (after_changed_cluster_size - 1) / 2;
    double possible_out_edge = after_changed_cluster_size * (species_num_nodes[spe] - after_changed_cluster_size);
    for (int i=0; i<species_network[spe]->node_adjList[perturb_node].size(); i++)
    {
        int node = species_network[spe]->node_adjList[perturb_node][i];
        if (spe_cluster_nodes[spe][cluster].find(node) != spe_cluster_nodes[spe][cluster].end())
        {
            in_edge_changed += 2;
            /* perturb_node moves into new_cluster,
             * so old_in_edge decrease decrease/compensate out_edge_changed */
            if (is_new_state)
                out_edge_changed -= 1;
            /* perturb_node moves out of old_cluster, so
             * old_in_edge becomes new_out_edge and thus decrease/compensate out_edge_changed */
            if (!is_new_state)
                out_edge_changed -= 1;
        } else {
            out_edge_changed += 1;
        }
    }

    /* compute score */
    double new_in_edge, new_out_edge, new_local_density, new_relative_density, new_score, new_degree, new_x_log_x;
    if (after_changed_cluster_size >= 1) // cluster is not empty, then count in, out edges
    {
        if (!is_new_state)
        {
            new_in_edge = species_cluster_inEdges[spe][cluster] - in_edge_changed;
            new_out_edge = species_cluster_outEdges[spe][cluster] - out_edge_changed;
        } else {
            new_in_edge = species_cluster_inEdges[spe][cluster] + in_edge_changed;
            new_out_edge = species_cluster_outEdges[spe][cluster] + out_edge_changed;
        }

        if (new_in_edge >= 1) // cluster has in edge, then compute inDens, outDens, score
        {
            new_degree = new_in_edge + new_out_edge;
            new_x_log_x = (new_degree / (double)species_network[spe]->num_edges) * log(new_degree / (double)species_network[spe]->num_edges);
//            new_local_density = new_in_edge / possible_in_edge;

            new_relative_density = new_in_edge / new_degree;
            new_score = - new_x_log_x * new_relative_density;
        } else {
            new_in_edge = 0.0;
            new_local_density = 0.0;
            new_relative_density = 0.0;
            new_score = 0.0;
        }
    } else {
        new_in_edge = 0.0;
        new_out_edge = 0.0;
        new_local_density = 0.0;
        new_relative_density = 0.0;
        new_score = 0.0;
    }
    /* compute delta score after purturb */
    cluster_delta_score = new_score - species_cluster_coexpress_cost[spe][cluster];

    /* update cluster info */
    species_cluster_size[spe][cluster] = after_changed_cluster_size;
    species_cluster_inEdges[spe][cluster] = new_in_edge;
    species_cluster_outEdges[spe][cluster] = new_out_edge;
    species_cluster_coexpress_cost[spe][cluster] = new_score;

//    fprintf(stderr, "spe %d, cluster %d, cluster size %g, indensity %g, outdensity %g, new_out_edge %g, possible_out_edge %g, new_in_edge %g\n", spe, cluster, after_changed_cluster_size, new_inDens, new_outDens, new_out_edge, possible_out_edge, new_in_edge);

    return - cluster_delta_score;
}

/*============================================
 delta cost after perturb of a node in species

 1) if either the prev or new cluster of this node is noise cluster (0)
 then, compute the old and new score (delete this node from old cluster
 and add to new cluster) of this species, define delta_score = new - old,
 and return -delta_score

 2) if both of prev and new cluster of this node are not noise cluster,
 then do the same thing as above to simplify computation steps. (it is possible
 to only score the new and old cluster, but that needs to take care of
 adaptive nosie cost term)
 ============================================*/
double Clustering::get_coexpress_delta_cost(SpeciesNetwork **species_network)
{
    double delta_score = 0.0; // delta_cost = -delta_score

    int perturb_node = UndoLog[undo_log_size-1].node;
    int spe = UndoLog[undo_log_size-1].spc;
    int old_cluster = UndoLog[undo_log_size-1].oldstate;
    int new_cluster = node_assigned_cluster[spe][perturb_node]; // perturb assigns a new cluster to node_assigned_cluster, that's why we need to keep it

    double old_score = get_spe_score(spe,species_network);

    spe_cluster_nodes[spe][old_cluster].erase(perturb_node);
    spe_cluster_nodes[spe][new_cluster].insert(perturb_node);

    double new_score = get_spe_score(spe,species_network);

    delta_score = new_score - old_score;
    return -delta_score;
}

/*===========================
 score clustering of a species
 1) used as a helper function for 'get_coexpress_delta_cost' to score
    new and old clusterings
 ============================*/
double Clustering::get_spe_score(int spe, SpeciesNetwork **species_network)
{
    double spe_score = 0.0;

    unsigned long noise_cluster_size = spe_cluster_nodes[spe][0].size();
    // ignore noise cluster here
    for (int clus = 1; clus < num_clusters; clus++)
    {
        int in_edges = 0;
        int out_edges = 0;
        unsigned long cluster_size = spe_cluster_nodes[spe][clus].size();

        for (unordered_set<int>::iterator node= spe_cluster_nodes[spe][clus].begin(); node!=spe_cluster_nodes[spe][clus].end(); node++)
        {
            for (int i=0; i<species_network[spe]->node_adjList[*node].size(); i++)
            {
                int adj_gene = species_network[spe]->node_adjList[*node][i];
                // adj_gene in the same cluster as cur gene
                if (spe_cluster_nodes[spe][clus].find(adj_gene)  != spe_cluster_nodes[spe][clus].end())
                {
                    in_edges++;
                }
                else
                {
                    // adj_gene not in noise cluster
                    if (spe_cluster_nodes[spe][0].find(adj_gene) == spe_cluster_nodes[spe][0].end())
                    {
                        out_edges++;
                    }
                }
            }
        }
        in_edges /=2;
        // cluster with size <= 1 will have possible_in_edge =0
        if (cluster_size >=2 && (in_edges > 0 || out_edges > 0))
        {
            unsigned long possible_in_edge = (cluster_size) * (cluster_size - 1) / 2;
            unsigned long possible_out_edge = ( species_num_nodes[spe] - (noise_cluster_size) - (cluster_size) ) * (cluster_size);

            double in_density = double(in_edges) /possible_in_edge;
            double out_density = double(out_edges) / possible_out_edge;

            double inOutRatio = (in_density / (in_density + out_density)) * double(cluster_size);

            spe_score += inOutRatio;
//            fprintf(stderr, "spe %d, cluster size %d, indensity %g, outdensity %g\n", spe, cluster_size, in_density, out_density);
        }
    }
//    double noise_cost = (double)noise_cluster_size * (spe_score / species_num_nodes[spe]);
    double noise_cost = 0.5 * (double)noise_cluster_size;
    spe_score +=noise_cost;
//    if (SA_counter % HOWMANYITERTOSHOW == 0)
//        fprintf(stderr, "spe %d, noise cluster size %d, spe_score %g\n", spe, noise_cluster_size, spe_score);
    return spe_score;
}


/*==============================
 score a clustering of a network
 1) loop through each node in a speices and assign them to corresponding clusters,
    count in-cluster edges, out-cluster edges, in-cluster nodes
 2) loop through each cluster of a species, compute size(cluster) * In-density/(In-density + Out-density)
 3) for each species, sum up InOutRatio score and compute adaptive noise score
 4) return -score as cost
 ==============================*/
double Clustering::get_coexpress_cost(SpeciesNetwork **species_network)
{
    // init total score = 0
    double score = 0.0;

    for (int spe = 0; spe < num_species; spe++)
    {
        unordered_map<int, double> total_cluster_outEdge; // #out-cluster edges
        unordered_map<int, double> total_cluster_inEdge; // #in-cluster edges
        map<int, double> total_cluster_nodes; // #nodes in cluster, ordered
        int noise_clus_size = 0; // #nodes in noise cluster

        /*===== cost in normal clusters =====*/
        /*
         * loop through node_i in spe
         * */
        for (int i=0; i<species_num_nodes[spe]; i++)
        {
            // skip noise cluster
            if (node_assigned_cluster[spe][i] != 0)
            {
                // initialize in-cluster total edge = 0
                total_cluster_inEdge[node_assigned_cluster[spe][i]] += 0.0;
                total_cluster_outEdge[node_assigned_cluster[spe][i]] += 0.0;
                // increment #nodes in clus
                total_cluster_nodes[node_assigned_cluster[spe][i]] +=1.0;
                /*
                 * loop through node_j connects to node_i
                 * so that the out_edge can be double counted, which is needed
                 */
                for (int k=0; k<species_network[spe]->node_adjList[i].size(); k++)
                {
                    int j = species_network[spe]->node_adjList[i][k];
                    // skip noise cluster
                    if (node_assigned_cluster[spe][j] != 0)
                    {
                        // i and j in the same cluster
                        if (node_assigned_cluster[spe][i] == node_assigned_cluster[spe][j])
                        {
                            total_cluster_inEdge[node_assigned_cluster[spe][i]]++;

                        }
                            // i and j are in diff cluster
                        else
                        {
                            total_cluster_outEdge[node_assigned_cluster[spe][i]] ++;
                        }
                    }
                }
            }
            else
            { // if cur cluster is noise cluster
                // increment #nodes in noise cluster
                noise_clus_size ++;
            }
        }

        /*
         * compute Schaeffer score
         * */
        double spe_score = 0.0;
        for (map<int, double>::iterator x=total_cluster_nodes.begin(); x!=total_cluster_nodes.end(); x++)
        {
            int cluster = x->first;
            double cluster_size =x->second;
            // #in-cluster edges > 0, do not need to worry about 'noise cluster' since it does not take into account above
            double possible_in_edge = 0.0;
            double x_log_x = 0.0, relative_density = 0.0, total_degree = 0.0;
            if (total_cluster_inEdge[cluster] > 0)
            {
//                total_cluster_inEdge[cluster] = total_cluster_inEdge[cluster]/2;
                total_cluster_inEdge[cluster] = total_cluster_inEdge[cluster];
                total_degree = total_cluster_inEdge[cluster] + total_cluster_outEdge[cluster];
                x_log_x = (total_degree / (double)species_network[spe]->num_edges) * log(total_degree / (double)species_network[spe]->num_edges);
                relative_density = total_cluster_inEdge[cluster] / total_degree;
                // cost of a spc = in-density * size / (in-density + out-density)
                species_cluster_coexpress_cost[spe][cluster] = - x_log_x * relative_density;
                spe_score += species_cluster_coexpress_cost[spe][cluster];

//                fprintf(stderr, "cluster: %d, size: %d, in_edge: %d, out_edge: %d, possible_in_edge: %d\n", cluster, cluster_size, total_cluster_inEdge[cluster], total_cluster_outEdge[cluster], possible_in_edge);
            }
            /* init cluster size, inEdge, outEdge */
            species_cluster_size[spe][cluster] = cluster_size;
            species_cluster_inEdges[spe][cluster] = total_cluster_inEdge[cluster];
            species_cluster_outEdges[spe][cluster] = total_cluster_outEdge[cluster];

//            fprintf(stderr, "cluster: %d, size: %g, in_edge: %g, out_edge: %g, possible_in_edge: %g\n", cluster, cluster_size, total_cluster_inEdge[cluster], total_cluster_outEdge[cluster], possible_in_edge);
        }
        // add spe_score to total score
        score += spe_score;

        /* cost of noise cluster = #nodes in noise cluster * (spc_cost / size)
         * adaptive weight for noise cluster
         */
//        score += (double)noise_clus_size * (spe_score / species_num_nodes[spe]);

        score += (double)noise_clus_size * 0.5; // 0.5 is the mean in/(in+out) density in a random graph
//        fprintf(stderr, "spe %d, noise cluster size %d, spe_score/#nodes %g\n", spe, noise_clus_size, (spe_score / species_num_nodes[spe]));
    }

    // return -score as cost
    return - score;
}

double Clustering::get_orth_delta_cost(Orthology *orthology, SpeciesNetwork **species_network, double curr_cost)
{
    buffer_orth_cost(curr_cost);
    double delta_score = 0.0;
    int perturb_node = UndoLog[undo_log_size-1].node;
    int spe = UndoLog[undo_log_size-1].spc;
    int old_cluster = UndoLog[undo_log_size-1].oldstate;
    int new_cluster = node_assigned_cluster[spe][perturb_node]; // perturb assigns a new cluster to node_assigned_cluster, that's why we need to keep it

//    fprintf(stderr, "===== get_orth_delta_cost spe %d, node %s ======\n", spe, (species_network[spe]->uniqId_to_geneName[perturb_node]).c_str());
    for (unordered_set<string>::iterator x=orthology->spe_orth_node_adjList[spe][perturb_node].begin(); x!=orthology->spe_orth_node_adjList[spe][perturb_node].end(); x++)
    {
        vector<string> spe_gene = split(*x,'.');
        int spe2 = stoi(spe_gene[0]);
        int gene2 = stoi(spe_gene[1]);

//        fprintf(stderr,"adj nodes: %s in spe %d\n", (species_network[spe2]->uniqId_to_geneName[gene2]).c_str(), spe2);

        double edge_score = (orthology->weighted_orthEdge[spe][spe2][perturb_node] + orthology->weighted_orthEdge[spe2][spe][gene2]) / 2.0;

        if (node_assigned_cluster[spe2][gene2] == old_cluster)
        {
            delta_score -= edge_score;
//            fprintf(stderr,"old_cluster orth nodes: %d in spe %d\n", gene2, spe2);
        }
        else if (node_assigned_cluster[spe2][gene2] == new_cluster)
        {
            delta_score += edge_score;
//            fprintf(stderr,"new_cluster orth nodes: %d in spe %d\n", gene2, spe2);
        }
    }
    return - delta_score;
}

void Clustering::buffer_orth_cost(double cost)
{
    old_orth_cost = cost;
}

/*===============================
 get_orth_cost of all pairs of species
 1) for  from 1 to total spe
        for spe2 from +1 to total spe
            for each ortholog edge
                if both end nodes are in the same cluster, orth_score++
                if only one end node is in noise cluster, noise_score++
 2) return -score as cost
 ===============================*/
double Clustering::get_orth_cost(Orthology *orthology)
{
    double orthterms = 0.0; // orth cost without noise clusters

//    double total_noise_nodes = 0.0;
    double orth_noise_term = 0.0; // orth cost of noise clusters
//    double total_gene_among_spc = 0.0; // total number of genes among three spes
    double total_ortho_noise = 0.0;

    for (int spe1=0; spe1 < num_species; spe1++)
    {
//        total_gene_among_spc += species_num_nodes[spc1]; // sum up total genes
//        total_noise_nodes += spe_cluster_nodes[spc1][0].size(); // sum up noise nodes

        for (int spe2=spe1+1; spe2 < num_species; spe2++) // avoid double-counting ortholog edges
        {
            /*=============
             write a helper funciton to compute the cost for a pair of species
             ==============*/
            for (unordered_map< Orthology::NodePair*, int, Orthology::NodePairHasher, Orthology::eqNodePair >::iterator x= orthology->orth_edges[spe1][spe2].begin(); x!=orthology->orth_edges[spe1][spe2].end(); x++)
            {
                Orthology::NodePair *np = (Orthology::NodePair *)(x->first);
                int n1 = np->id1; // node1
                int n2 = np->id2; // node2



                // count edges where both end nodes are in the same cluster, not counting nodes in noise cluster
                if (node_assigned_cluster[spe1][n1] == node_assigned_cluster[spe2][n2] && node_assigned_cluster[spe1][n1] != 0 && node_assigned_cluster[spe2][n2] != 0)
                {
                    double weight = (orthology->weighted_orthEdge[spe1][spe2][n1] + orthology->weighted_orthEdge[spe2][spe1][n2]) / 2.0;
                    orthterms += weight;
//                    fprintf(stderr, "spe%d, spe%d, gene%d, gene%d, weight %g\n", spe1, spe2, n1, n2, weight);

                    /*===========
                     1. change this obj function to (((_clus_size + spe2_clus_size) / 2) * in-density) / (in-density + out-density)
                     2. add noise version later
                     ============*/
                }
                /* For noise nodes, only count ortholog edges that cross noise cluster. When both end of an ortho edge are in noise cluster, chances are these pair of nodes have co-expression edges, and so should be put in normal clusters.
                 This rewards:
                 1) all 3 ortho edge within a normal cluster
                 2) one end of ortho edge in noise cluster, and the other in normal cluster, max reward happens when end node of the other two species are in the same normal cluster: 2(cross noise) + 1(within normal) = 3
                 3) when 1 and 2 become a tie, then In/Out-density ratio will decide if a noise should be moved out of noise cluster
                 */
//                if (node_assigned_cluster[spc1][n1] != node_assigned_cluster[spe2][n2] && (node_assigned_cluster[spc1][n1] == 0 || node_assigned_cluster[spe2][n2] == 0))
                if ((node_assigned_cluster[spe1][n1] == 0 || node_assigned_cluster[spe2][n2] == 0))
                {
                    total_ortho_noise += (orthology->weighted_orthEdge[spe1][spe2][n1] + orthology->weighted_orthEdge[spe1][spe2][n2]) / 2.0;
                }
            }
        }
    }
    orth_noise_term = 0.5 * total_ortho_noise;

    // this noise term can also be scaled by number of orth edges between two species
    orthterms += orth_noise_term;

//    if (SA_counter % HOWMANYITERTOSHOW == 0)
//    fprintf(stderr, "total ortho noise %g\n", total_ortho_noise);
    //    fprintf(stderr, "cc = %g\tcccost = %g\n",lambda, lambda*orthterms); // prev log, print orth term
    return -orthterms;
}


/*================================================================
 undo_perturb
 1) recover perturbed node to its old cluster
 2) for spe_cluster_nodes, delete perturbed node from new cluster,
    add perturbed node back to old cluster
 ==============================================================*/
void Clustering::undo_perturb(int new_state)
{
    for (int i = undo_log_size; i>0; i--)
    {
        int spe = UndoLog[i-1].spc;
        int node = UndoLog[i-1].node;
        int old_state = UndoLog[i-1].oldstate;
        node_assigned_cluster[spe][node] = old_state;
        spe_cluster_nodes[spe][new_state].erase(node);
        spe_cluster_nodes[spe][old_state].insert(node);
        changed_nodes[to_string(static_cast<long long>(spe)) + "." + to_string(static_cast<long long>(node))] = old_state;

        species_cluster_size[spe][old_state] = old_cluster_size_before_changed;
        species_cluster_inEdges[spe][old_state] = old_cluster_inEdges_before_changed;
        species_cluster_outEdges[spe][old_state] = old_cluster_outEdges_before_changed;
        species_cluster_coexpress_cost[spe][old_state] = old_cluster_coexpress_cost_before_changed;

        species_cluster_size[spe][new_state] = new_cluster_size_before_changed;
        species_cluster_inEdges[spe][new_state] = new_cluster_inEdges_before_changed;
        species_cluster_outEdges[spe][new_state] = new_cluster_outEdges_before_changed;
        species_cluster_coexpress_cost[spe][new_state] = new_cluster_coexpress_cost_before_changed;

        current_orth_cost = old_orth_cost;
    }

    undo_log_size = 0;
}

/*===========================
 Perturb
 1) randomly select a species
 2) randomly select a node
 3) randomly select a new cluster
 4) assign this to the new cluster
 ==========================*/
int Clustering::Perturb(SpeciesNetwork **species_network, int numchanges = 1)
{
    int newstate = 0;
    for (int i=0; i<numchanges; i++)
    {
        // choose a random species
        int spc = rand() % num_species;
        // choose a random node
        int node = rand() % (species_num_nodes[spc]); // node is a uniq_ID instead of gene_ID
        int oldstate = node_assigned_cluster[spc][node];

        // choose a random cluster
        while (1)
        {
            // MAXLOGSIZE = 1000
            if (undo_log_size > MAXLOGSIZE) {
                printf("Error: Undo Log Size Limit reached\n");
                exit(1);
            }

            // choose a randome new cluster
            int rand_index = 0;
            if (has_noise_cluster)
            {
                rand_index = rand() % num_clusters; // if noise cluster included
            }
            else
            {
                rand_index = (rand() % (num_clusters-1))+1; // no cluster0
            }
            newstate = rand_index;

            /* if new cluster is the same as old,
            jump to the end of loop, not increase undo_log_size */
            int min_clus_size = spe_cluster_nodes[spc][newstate].size();

            // not assign node to a cluster if it is too small
//            if (newstate == oldstate or min_clus_size <= species_num_nodes[spc]/100) continue;
            if (newstate == oldstate) continue;

            // assign node to new cluster
            node_assigned_cluster[spc][node] = newstate;
            changed_nodes[to_string(static_cast<long long>(spc)) + "." + to_string(static_cast<long long>(node))] = newstate;

            if (SA_counter % HOWMANYITERTOSHOW == 0 or show_log)
//                fprintf(stderr,"Perturb spe %d, node %d, old %d, new %d\n",spc,node,oldstate,newstate); // prev log
                fprintf(stderr,"Perturb spe %d, node %d, old %d, new %d\n",spc,node,oldstate,newstate); // prev log


            // record the change for later undo
            UndoLog[undo_log_size].spc = spc;
            UndoLog[undo_log_size].node = node;
            UndoLog[undo_log_size].oldstate = oldstate;
            undo_log_size++;
            // stop the while loop
            break;
        }
    }
    return newstate;
}

vector<string> Clustering::split(const string &s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> tokens;
    while (getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

void Clustering::buffer_cluster_state(int size, double inEdges, double outEdges, double cost, bool is_new_state)
{
    if (is_new_state)
    {
        new_cluster_size_before_changed = size;
        new_cluster_inEdges_before_changed = inEdges;
        new_cluster_outEdges_before_changed = outEdges;
        new_cluster_coexpress_cost_before_changed = cost;
    }
    else
    {
        old_cluster_size_before_changed = size;
        old_cluster_inEdges_before_changed = inEdges;
        old_cluster_outEdges_before_changed = outEdges;
        old_cluster_coexpress_cost_before_changed = cost;
    }
}
