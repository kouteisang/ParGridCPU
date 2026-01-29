#ifndef FCSyncMem_H
#define FCSyncMem_H

#include "Graph/MultilayerGraph.h"
#include "FCCoreTree.h"


class FCSyncMem
{
private:
    /* data */
public:
    FCSyncMem(/* args */);
    ~FCSyncMem();


    static void Execute(MultilayerGraph &mg, FCCoreTree &tree);

    static void constructCoreSync(coreNodeP *node, uint k, uint lmd, uint n_vertex, uint n_layer, bool* valid, int** degs, int total, bool serial, coreNodeP* father);

    static void PeelSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, coreNodeP* father);

    static void PathSerialSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, coreNodeP* father);

    static void BuildSubFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, int &total, coreNodeP* father);
// PathByK(mg, root->degs, root->k+1, root->lmd, leftChild, root->valid, root->total);
    static void PathByK(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int &total, coreNodeP* father);


    // The following are mix solution
    // static void PathByKMix(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, vector<bool>& valid, int &total);
    // static void PeelSyncMix(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, vector<bool>& valid, bool serial, int &total);
    // static void ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree); 
    // static void BuildSubFCTreeSyncMix(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, vector<bool>& valid, int &total);
    // static void BuildSubcFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, vector<bool>& valid, uint* invalid, uint* cnts);
};




#endif