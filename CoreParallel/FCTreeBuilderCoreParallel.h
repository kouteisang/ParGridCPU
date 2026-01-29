// #ifndef FCTreeBuilderCoreParallel_H
// #define FCTreeBuilderCoreParallel_H

// #include "Graph/MultilayerGraph.h"
// #include "FCCoreTree.h"

// class FCTreeBuilderCoreParallel
// {
// private:
//     /* data */
// public:
//     FCTreeBuilderCoreParallel(/* args */);
//     ~FCTreeBuilderCoreParallel();

//     static void Execute(MultilayerGraph &mg, FCCoreTree &tree);
    
//     static void BuildSubFCTree(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts);
 
//     static void PathSerial(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

//     static void PathByK(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

//     // Core Parallel to Peel invalid vertex
//     static uint PeelInvalidInParallelByCount(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e, bool serial);

//     static uint PeelInvalidInParallelByCheck(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e, bool serial);

//     // useful common method
//     static bool check(uint **degs, uint u, uint k, uint lmd, uint n_layers);

//     static void constructCore(coreNodeP *node, uint k, uint lmd, uint new_e, uint n_vertex, uint n_layer, bool* valid, uint* invalid, uint** degs, uint* cnts, bool serial, bool wds);

//     // Mix strategy
//     static void ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree);

//     static void BuildSubFCTreeMix(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts);

//     static void PathSerialMix(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

// };



// #endif