// #ifndef FCTreeBuilderCoreParallelByK_H
// #define FCTreeBuilderCoreParallelByK_H


// #include "Graph/MultilayerGraph.h"
// #include "FCCoreTree.h"

// class FCTreeBuilderCoreParallelByK
// {
// private:
//     /* data */
// public:
//     FCTreeBuilderCoreParallelByK(/* args */);
//     ~FCTreeBuilderCoreParallelByK();

//     static void Execute(MultilayerGraph &mg, FCCoreTree &tree);
    
//     static void BuildSubFCTree(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts);
 
//     static void PathSerial(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

//     static void PathBylmd(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

//     // Core Parallel to Peel invalid vertex
//     static uint PeelInvalidInParallelByCount(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e, bool serial);

//     // useful common method
//     static bool check(int **degs, uint u, uint k, uint lmd, uint n_layers);

//     static void constructCore(coreNodeP *node, uint k, uint lmd, uint new_e, uint n_vertex, uint n_layer, bool* valid, uint* invalid, uint** degs, uint* cnts, bool serial);

//     // Mix strategy
//     static void ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree);

//     static void BuildSubFCTreeMix(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts);

//     static void PathByLmdTask(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e);

// };


// #endif
