#ifndef FCSyncLeft_H
#define FCSyncLeft_H

#include "Graph/MultilayerGraph.h"
#include "FCCoreTree.h"
#include "Util/UtilStats.h"


class FCSyncLeft
{
private:
    /* data */
public:
    FCSyncLeft(/* args */);
    ~FCSyncLeft();


    static void Execute(MultilayerGraph &mg, FCCoreTree &tree);

    static void constructCoreSync(coreNodeP *node, uint k, uint lmd, uint n_vertex, uint n_layer, bool* valid, int** degs, int total, bool serial);

    static void PeelSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats);

    static void PathSerialSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats);

    static void BuildSubFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, int &total, UtilStats& stats);
// PathByK(mg, root->degs, root->k+1, root->lmd, leftChild, root->valid, root->total);
    static void PathByK(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int &total, UtilStats& stats);


    // The following are mix solution
    static void PathByKMix(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int &total, UtilStats& stats);
    static void PeelSyncMix(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats);
    static void ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree); 
    static void BuildSubFCTreeSyncMix(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, int &total, UtilStats& stats);
    // static void BuildSubcFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts);
};




#endif