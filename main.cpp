#include "Clustering.h"
#include "SpeciesNetwork.h"
#include "Orthology.h"
#include <string.h>

using namespace std;


int main (int argc, char * const argv[]) {

    if (argc < 7)
    {
        printf("usage: %s "
                       "<num_clusters> "
                       "<num_trials> "
                       "<coupling_constant> "
                       "<orth_file> "
                       "<num_species> "
                       "<spc1nw> "
                       "<spc2nw> ... "
                       "-s:string <start_clustering_file> "
                       "-t:double <max_temp> "
                       "-n <clustering with noise cluster> \n", argv[0]);
        exit(1);
    }

    clock_t begin = clock();
    /*============ read files and parameters ============*/
    int argbase = 1; // count for argv
    int num_clusters = atoi(argv[argbase++]); // num clusters
    num_clusters +=1; // add one more cluster: 'noise cluster 0'
    int num_trials = atoi(argv[argbase++]); // number of trials
    double coupling_constant = atof(argv[argbase++]); // coupling_constant
    char *orth_file = argv[argbase++]; // ortholog files
    int num_species = atoi(argv[argbase++]); // number of species

    /* init species co-expression network */
    SpeciesNetwork **species_network = new SpeciesNetwork*[num_species];
    for (int i=0; i<num_species; i++)
    {
        char *spe_network_file = argv[argbase++];
//        species_network[i] = read_network_file(spc); // init network for a species
        species_network[i] = new SpeciesNetwork(spe_network_file);
    }



    /* read parameters: startclusfn '-s' and max_temp '-t', noise_cluster '-n', weight for noise term '-lambda' */
    char startclusfn[1024]; startclusfn[0] = 0;
    double max_temp = -1;
    bool has_noise_cluster = false;
    double lambda = 1.0;
    while (argbase < argc)
    {
        // ground-truth clutering
        if (!strcmp(argv[argbase], "-s"))
        {
            argbase++;
            strcpy(startclusfn, argv[argbase++]);
            continue;
        }
        // starting tempreture
        if (!strcmp(argv[argbase], "-t"))
        {
            argbase++;
            max_temp = atof(argv[argbase++]);
            continue;
        }
        // has noise cluster
        if (!strcmp(argv[argbase], "-n"))
        {
            argbase++;
            has_noise_cluster = true;
            continue;
        }
    }



    // create ortholog obj using (ortholog file, network arr, number of spe)
//    Orthology *orth = read_orthology_file(orth_file, species_network, num_species);

    /* init orthologous */
    Orthology *orthology = new Orthology(orth_file, species_network, num_species);

    /*======= clustering ==========*/
    // create arr of size num_trials for clustering obj
    Clustering **clustering = new Clustering *[num_trials];

    for (int i=0; i<num_trials; i++)
    {
        // init clustering for a trial
        clustering[i] = new Clustering(species_network,
                                       num_species,
                                       orthology,
                                       num_clusters,
                                       coupling_constant,
                                       has_noise_cluster);

        // use ground-truth clustering label as start point
        if (startclusfn[0]!=0)
        {
            fprintf(stderr,"Using seed clustering file %s\n",startclusfn);
            clustering[i]->pre_set(startclusfn, orthology, species_network);
        }

        // if starting tempreture used, set max_temp
        if (max_temp >= 0) clustering[i]->set_max_temp(max_temp);

        // simulated annealing to find parameters with best cost
        clustering[i]->learn_ground_state(species_network, orthology);
    }
    /*==== go through all trails, and record the best cost ====*/
    double best_cost = 0.0;
    int best_trial = -1;

    for (int i=0; i<num_trials; i++)
    {
        if (clustering[i]->best_total_cost < best_cost) // cost is negative
        {
            best_trial = i;
            best_cost = clustering[i]->best_total_cost;
        }
    }
    /*====== print best clustering found =======*/
    printf("Best clustering, with cost %g is:\n", best_cost);
    clustering[best_trial]->Print(species_network);

    /*====== print running time =========*/
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    fprintf(stderr,"Total time: %f s\n",elapsed_secs);

    /* free memory */
    for (int i=0; i<num_trials; i++)
    {
        delete clustering[i];
    }
    delete[] clustering;
    for (int i=0; i<num_species; i++)
    {
        delete species_network[i];
    }
    delete[] species_network;
    delete orthology;
    return 0;
}
