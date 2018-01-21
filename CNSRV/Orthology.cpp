//
// Created by 杨巍 on 5/18/16.
//

#include "Orthology.h"

/**
 * Destructor
 */
Orthology::~Orthology()
{
    fprintf(stderr, "Orthology destructor\n");
    for (int i=0;i<num_species;i++)
    {
        delete speciesId_to_name[i];
        for (int j=0;j<num_species;j++)
        {
            orth_edges[i][j].clear();
            weighted_orthEdge[i][j].clear();
        }
        delete[] orth_edges[i];
        delete[] weighted_orthEdge[i];
    }
    delete[] speciesId_to_name;
    delete[] orth_edges;
    delete[] weighted_orthEdge;
    speciesName_to_id.clear();
}

/**
 * Default Constructor
 */
Orthology::Orthology() { }

/**
 * Constructor
 * @param orth_file
 * @param species_network
 * @param num_species
 */
Orthology::Orthology(char *orth_file, SpeciesNetwork **species_network, int num_species) {
    fprintf(stderr, "Orthology....\n");
    this->num_species = num_species;
    orth_edges = new unordered_map<NodePair*, int, NodePairHasher, eqNodePair> *[num_species];
    weighted_orthEdge = new unordered_map<int, double> *[num_species];

    // pair up species for orthologous edges
    for (int i=0; i<num_species; i++)
    {
        orth_edges[i] = new unordered_map<NodePair *, int, NodePairHasher, eqNodePair> [num_species];
        weighted_orthEdge[i] = new unordered_map<int, double> [num_species];

        /* init adjList for orth_edges */
        vector< unordered_set<string> > tmp_orth_node_adjList;
        for (int j=0;j<species_network[i]->num_nodes;j++)
        {
            unordered_set<string> tmp_node_adjList;
            tmp_orth_node_adjList.push_back(tmp_node_adjList);
        }
        spe_orth_node_adjList.push_back(tmp_orth_node_adjList);
    }

    // create numerical ids for species names
    speciesId_to_name = new char *[num_species];
    for (int i=0; i<num_species; i++)
    {
        speciesId_to_name[i] = new char[LINELEN];
        strcpy(speciesId_to_name[i], species_network[i]->species_name);
        speciesName_to_id[(char *)species_network[i]->species_name] = i;
    }
    // read the file and populate table
    FILE *F = fopen(orth_file, "r");

    if (F == NULL)
    {
        fprintf(stderr, "cannot open file %s", orth_file);
        exit(1);
    }

    unordered_set<string> seen_line;
    char line[LINELEN];
    while (fgets(line, LINELEN, F))
    {
        char spe1[LINELEN];
        char spe2[LINELEN];
        char gene_id1[LINELEN];
        char gene_id2[LINELEN];
        sscanf(line, "%s %s %s %s",(char *)&spe1, (char *)&spe2, (char *)&gene_id1, (char *)&gene_id2);
        /* char to string */

        pair<int, int> species = get_spe_pair(spe1, spe2,orth_file);
        pair<int, int> genes = get_gene_pair(species_network, gene_id1, gene_id2, species);
        if (genes.first == -1 and genes.second == -1) continue;

        Orthology::NodePair *np1 = Orthology::get_node_pair(genes.first, genes.second);
        orth_edges[species.first][species.second][np1] = 1; // assign an ortho edge for g1,g2 in spe1,spe2, this double counted edges, should only count one direction??

        Orthology::NodePair *np2 = Orthology::get_node_pair(genes.second, genes.first);
        orth_edges[species.second][species.first][np2] = 1;

        string key1 = to_string(static_cast<long long>(species.first)) + "." + to_string(static_cast<long long>(species.second)) + "." + to_string(static_cast<long long>(genes.first)) + "." + to_string(static_cast<long long>(genes.second));
        if (seen_line.find(key1) == seen_line.end())
        {
            weighted_orthEdge[species.first][species.second][genes.first] +=
                    1.0; // sum edges for g1 in spe1,spe2, why only g1???
            seen_line.insert(key1);
        }

        string key2 = to_string(static_cast<long long>(species.second)) + "." + to_string(static_cast<long long>(species.first)) + "." + to_string(static_cast<long long>(genes.second)) + "." + to_string(static_cast<long long>(genes.first));
        if (seen_line.find(key2) == seen_line.end())
        {
            weighted_orthEdge[species.second][species.first][genes.second] += 1.0;
            seen_line.insert(key2);
        }

//        if (species.first == 2 and stoi(species_network[2]->uniqId_to_geneName[genes.first]) == 108)
//            fprintf(stderr, "----- orth edge to spe2 gene108: spe %d, gene %d, line %s ------\n", species.second, stoi(species_network[species.second]->uniqId_to_geneName[genes.second]), line);
        /* add orth edges to adjList
         * do not count species.second because orth file is already bi-direction*/
        spe_orth_node_adjList[species.first][genes.first].insert(to_string(static_cast<long long>(species.second)) + "." + to_string(static_cast<long long>(genes.second)));
        spe_orth_node_adjList[species.second][genes.second].insert(to_string(static_cast<long long>(species.first)) + "." + to_string(static_cast<long long>(genes.first)));
    }
    fclose(F);

    normalize_edge_weight();
}



Orthology::NodePair *Orthology::get_node_pair(int gene1, int gene2)
{
    NodePair *np = new NodePair;
    np->id1 = gene1;
    np->id2 = gene2;
    return np;
}


void Orthology::normalize_edge_weight()
{
    for (int spe1=0; spe1 < num_species; spe1++)
    {
        for (int spe2=0; spe2 < num_species; spe2++) // one direction
        {
//            unordered_map<int, double>::iterator iter;
//            for ( iter = o->weighted_orthEdge[spc1][spc2].begin(); iter != o->weighted_orthEdge[spc1][spc2].end(); ++iter)
            for (unordered_map<int, double>::iterator x=weighted_orthEdge[spe1][spe2].begin(); x!=weighted_orthEdge[spe1][spe2].end(); x++)
            {
                int gene = x->first, count = x->second;
                double orth_effect = 1.0 / count;
                weighted_orthEdge[spe1][spe2][gene] = orth_effect;
//                fprintf(stderr, "spe%d, spe%d, gene%d, weight %g\n", spe1, spe2, gene, weighted_orthEdge[spe1][spe2][gene]);
            }
        }
    }
}

pair<int, int> Orthology::get_gene_pair(SpeciesNetwork **species_network, char *g1, char *g2, pair<int, int>& species) {
    string str_gene_id1(g1);
    string str_gene_id2(g2);
//    if (stoi(str_gene_id1) == 108 and species.first == 2)
//        fprintf(stderr, "spe%d, spe%d, g%d, g%d\n", species.first, species.second, stoi(str_gene_id1), stoi(str_gene_id2));

    // skip genes that are not found in co-expression network
    if (species_network[species.first]->geneName_to_uniqId.find(str_gene_id1) == species_network[species.first]->geneName_to_uniqId.end()) {
        return make_pair(-1,-1);
    }

    if (species_network[species.second]->geneName_to_uniqId.find(str_gene_id2) == species_network[species.second]->geneName_to_uniqId.end()) {
        return make_pair(-1,-1);
    }

    // get converted geneID in (1,N)
    int gene1_id = species_network[species.first]->geneName_to_uniqId[str_gene_id1];
    int gene2_id = species_network[species.second]->geneName_to_uniqId[str_gene_id2];
    return make_pair(gene1_id, gene2_id);
}

pair<int, int> Orthology::get_spe_pair(char *spe1, char *spe2, char *orth_file)
{
    string str_spe1(spe1);
    string str_spe2(spe2);

    // convert spc ids to numerics
    if (speciesName_to_id.find(str_spe1) == speciesName_to_id.end())
    {
        fprintf(stderr, "Error: orthology file %s has species name that wasn't seen in network files\n", orth_file);
        exit(1);
    }

    if (speciesName_to_id.find(str_spe2) == speciesName_to_id.end()) {
        fprintf(stderr, "Error: orthology file %s has species name that wasnt seen in network files\n", orth_file);
        exit(1);
    }

    int species1_id = speciesName_to_id[str_spe1];
    int species2_id = speciesName_to_id[str_spe2];
    return make_pair(species1_id, species2_id);
}
