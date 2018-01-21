//
// Created by 杨巍 on 5/18/16.
//

#ifndef OCLODE_EFFICIENT_ORTHOLOGY_H
#define OCLODE_EFFICIENT_ORTHOLOGY_H

#include <cstdio>
#include <ctype.h>
#include <unordered_map>
#include <string>
#include <utility>
#include <sstream>
#include <unordered_set>
#include <ext/hash_map>
#include <string.h>

#include "SpeciesNetwork.h"

using namespace std;

#define LINELEN 1024

class Orthology {


public:
    struct NodePair
    {
        int id1;
        int id2;
    };

    struct eqNodePair
    {
        bool operator()(NodePair* n1, NodePair* n2) const
        {
            return (n1->id1 == n2->id1 && n1->id2 == n2->id2);
        }
    };


    class NodePairHasher {
    public:
        size_t operator()(const NodePair *r) const
        {
            char line[LINELEN];
            sprintf(line,"%d,%d",r->id1, r->id2);
            return h(line);
        };

    private:
        __gnu_cxx::hash<char*> h;
    };


private:
    int num_species;
    char **speciesId_to_name; // mapping speciesID -> speciesName

public:
    unordered_map< NodePair*, int, NodePairHasher, eqNodePair > **orth_edges;
    unordered_map<int, double> **weighted_orthEdge; // hash table to store 1/d(#ortho genes from spe1 to spe2)
    unordered_map<string, int> speciesName_to_id; // mapping speciesName -> speciesID
    vector< vector< unordered_set<string> > > spe_orth_node_adjList; // string: spe.gene

/* methods */
public:
    Orthology();
    Orthology(char *orth_file, SpeciesNetwork **species_network, int num_species);
    ~Orthology();

private:
    NodePair* get_node_pair(int gene1, int gene2);
    void normalize_edge_weight();
    pair<int,int> get_spe_pair(char *spe1, char *spe2,char *orth_file);
    pair<int, int> get_gene_pair(SpeciesNetwork **species_network, char *g1, char *, pair<int,int>& species);
};


#endif //OCLODE_EFFICIENT_ORTHOLOGY_H
