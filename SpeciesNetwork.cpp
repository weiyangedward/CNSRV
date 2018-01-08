//
// Created by 杨巍 on 5/18/16.
//

#include "SpeciesNetwork.h"

/**
 * Destructor
 */
SpeciesNetwork::~SpeciesNetwork()
{
    fprintf(stderr, "SpeciesNetwork destructor\n");
    node_adjList.clear();
    geneName_to_uniqId.clear();
}

/**
 * default Constructor
 */
SpeciesNetwork::SpeciesNetwork()
{ }

/**
 * Constructor
 * @param network_file
 */
SpeciesNetwork::SpeciesNetwork(char *network_file)
{
    unique_node = 0;
    num_edges = 0;
    FILE *F = fopen(network_file, "r"); // open file
    if (F == NULL)
    {
        fprintf(stderr, "cannot open file %s", network_file);
        exit(1);
    }

    char line[LINELEN]; // buffer for each line

    fgets(line, LINELEN, F); // first line: species name
    sscanf(line, "%s", (char *) &species_name);


    fgets(line, LINELEN, F); // second line: number of nodes
    sscanf(line, "%d", &num_nodes);


    /* init network adjacency list */
    for (
            int i = 0; i < num_nodes; i++
            )
    {
        vector<int> adj_genes;
        node_adjList.push_back(adj_genes);
    }

    // a unique ID assigned each gene so that the all gene IDs will be mapped to 1,2,3,...,n
    /*====== read the rest of network file ======*/
    while (fgets(line, LINELEN, F))
    {
//        char *id1 = new char[LINELEN]; // buff for gene1
//        char *id2 = new char[LINELEN]; // buff for gene2
        char   id1[LINELEN];
        char   id2[LINELEN];
        double edge_weight; // edge weight
        // edge_weight = 1
        sscanf(line, "%s\t%s\t%lf", id1, id2, &edge_weight);
        string str_id1(id1); // copy id1 to a str
        string str_id2(id2); // copy id2 to a str

        /*======== convert gene_ID to uniq_ID =========*/
        hash_gene_id(str_id1, network_file); // convert gene1
        hash_gene_id(str_id2, network_file); // convert gene2

        /*====== build network using adj matrix and adj list =======*/
//        network[geneName_to_uniqId[str_id1]][geneName_to_uniqId[str_id2]] = edge_weight; // store edge weight between gene1 and gene2 in network adj matrix

        node_adjList[geneName_to_uniqId[str_id1]].push_back(geneName_to_uniqId[str_id2]); // store gene2 in the adj list of gene1, only do this for gene1 assuming that network file has both edges of: gene1->gene2 and gene2->gene1
        num_edges ++;
    }
    fclose(F); // close file

}


bool SpeciesNetwork::is_edge(int i, int j)
{
    if (i >= num_nodes || j >= num_nodes)
    {
        return false;
    }

    for (
            int k = 0; k < node_adjList[i].size(); k++
            )
    {
        int gene = node_adjList[i][k];
        if (gene == j)
        {
            return true;
        }
    }
    return false;
}


void SpeciesNetwork::hash_gene_id(string &gene_ID, char *network_file)
{
    if (geneName_to_uniqId.find(gene_ID) == geneName_to_uniqId.end()) // if gene1 is not seen
    {
        fprintf(
                stderr, "Assigning id %d to node %s in species %s\n", unique_node, gene_ID.c_str(), network_file
               ); // print out which gene is assigned what unique ID
//        uniqId_to_geneName[unique_node] = gene_ID; // create a new entry for (uniq_ID, gene_ID) in hashmap
        geneName_to_uniqId[gene_ID] = unique_node; // create a new entry for (gene_ID, gene_ID) in hashmap
        unique_node++; // increase uniq_ID by 1

        if (unique_node > num_nodes) // error: if uniq_ID > number of nodes
        {
            printf(
                    "Error: In file %s, read more unique node ids (%d) than the %d expected\n", network_file,
                    unique_node, num_nodes
                  );
            exit(1);
        }
    }

}
