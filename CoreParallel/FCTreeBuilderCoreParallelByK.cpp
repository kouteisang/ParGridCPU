// #include "FCTreeBuilderCoreParallelByK.h"


// FCTreeBuilderCoreParallelByK::~FCTreeBuilderCoreParallelByK(){
// }


// void FCTreeBuilderCoreParallelByK::constructCore(coreNodeP *node, uint k, uint lmd, uint new_e, uint n_vertex, uint n_layer, bool* valid, uint* invalid, uint** degs, uint* cnts, bool serial){
    
//     node->k = k;
//     node->lmd = lmd;
//     node->length = n_vertex - new_e;
//     node->e = new_e;
//     node->valid = new bool[n_vertex];
     

//     // Plan A store the valid array
//     memcpy(node->valid, valid, sizeof(bool) * n_vertex);
    

//     // TODO: Plan B store only the valid vertex index

//     // If serial store the deg info and the invalid info
//     if(serial){
//         node->degs = new int*[n_vertex];
//         for(uint v = 0; v < n_vertex; v ++){
//             node->degs[v] = new int[n_layer];
//             memcpy(node->degs[v], degs[v], n_layer * sizeof(uint));
//         }

//         node->invalid = new uint[n_vertex];
//         memcpy(node->invalid, invalid, n_vertex * sizeof(uint));
//         node->cnts = new uint[n_vertex];
//         memcpy(node->cnts, cnts, n_vertex * sizeof(uint));
//     }

//     // cout << "k = " << node->k << " lmd = " << node->lmd << " length = " << node->length << " new_e = " << node->e << endl;

// }


// uint FCTreeBuilderCoreParallelByK::PeelInvalidInParallelByCount(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e, bool serial){
//     uint n_vertex = mg.GetN(); // number of vertex
//     uint n_layers = mg.getLayerNumber();
//     uint cnt = 0;
//     uint s = e, new_e = 0;
//     uint **adj_lst;

//     // Debug for test
//     // omp_set_num_threads(1);

// #pragma omp parallel private(cnt) shared(cnts)
// {
//     uint tid = omp_get_thread_num();
//     uint maxthreads = omp_get_num_threads();

//     for(int v = tid; v < n_vertex; v += maxthreads){
//         cnt = 0;
//         if(valid[v] != 1){
//             cnts[v] = 0;
//             continue; // only process the valid vertex
//         } 
//         for(int l = 0; l < n_layers; l ++){
//             if(degs[v][l] >= k){
//                 cnt += 1;
//             }
//         }
//         cnts[v] = cnt;


//         if(cnts[v] < lmd){
//             valid[v] = 0;
//             int index;
//             #pragma omp atomic capture
//             index = e++;
//             invalid[index] = v;
//         }
//     }
// }
//     new_e = e;
//     // cout << "new_e = " << new_e << endl;   
//     uint vv;
//     while(s < e){
// #pragma omp parallel private(vv, adj_lst) shared(s, e, new_e, valid, invalid, degs, cnts, mg)
// {
//         #pragma omp for
//         for(uint j = s; j < e; j ++){
//             vv = invalid[j];
//             for(uint l = 0; l < n_layers; l ++){
//                 adj_lst = mg.GetGraph(l).GetAdjLst();
//                 for(uint i = 1; i <= adj_lst[vv][0]; i ++){
//                     uint u = adj_lst[vv][i]; // the neighbourhood
//                     if(valid[u] != 1) continue; // only process if u is valid
//                     // minus one and return the old value
//                     auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
//                     if(originDeg == k){
//                        auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
//                        if(originCnt == lmd && valid[u] == 1){
//                             #pragma omp critical
//                             {
//                                 invalid[new_e] = u; 
//                                 new_e ++;
//                             }
//                             valid[u] = 0; 
//                             cnts[u] = 0;
//                        } 
//                     }
//                 }
//             }
//         }
//         #pragma omp barrier
// }

//         s = e;
//         e = new_e;
//     }
//     if(n_vertex - new_e > 0){
//         constructCore(node, k, lmd, new_e, n_vertex, n_layers, valid, invalid, degs, cnts, serial);
//     }else{
//         node = nullptr;
//     }
   

//     return new_e;
// }


// void FCTreeBuilderCoreParallelByK::PathBylmd(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e){
//     uint n_vertex = mg.GetN(); // number of vertex
//     uint n_layers = mg.getLayerNumber(); // number of layer
    

//     uint new_e = PeelInvalidInParallelByCount(mg, degs, k, lmd, node, valid, invalid, cnts, e, false);

//     // means the (k, lambda)-constaint has the valid vertex
//     if(n_vertex - new_e > 0){
//         lmd += 1;
//         coreNodeP* rightChild = new coreNodeP();
//         node->right = rightChild;
//         PathBylmd(mg, degs, k, lmd, rightChild, valid, invalid, cnts, new_e);
//     }

// } 

// void FCTreeBuilderCoreParallelByK::PathSerial(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e){

//    uint n_vertex = mg.GetN(); // number of vertex
//    uint n_layers = mg.getLayerNumber(); // number of layer
   
//    uint new_e = PeelInvalidInParallelByCount(mg, degs, k, lmd, node, valid, invalid, cnts, e, true);

//    if(k == 10){
//         omp_set_num_threads(1);
//    }

//    // means the (k, lambda)-constaint has the valid vertex
//    if(n_vertex - new_e > 0){
//         k += 1;
//         coreNodeP* leftChild = new coreNodeP();
//         node->left = leftChild;
//         PathSerial(mg, degs, k, lmd, leftChild, valid, invalid, cnts, new_e);
//    }
// }

// void FCTreeBuilderCoreParallelByK::BuildSubFCTree(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts){

//     uint k = klmd[0];
//     uint lmd = klmd[1];
//     uint n_layers = mg.getLayerNumber();

//     uint n_vertex = mg.GetN(); // number of vertex

//     auto start_time_serial = omp_get_wtime(); 
//     PathSerial(mg, degs, k, lmd, node, valid, invalid, cnts, 0);
//     auto end_time_serial = omp_get_wtime(); 

//     double elapsed_time_serial = end_time_serial - start_time_serial;
//     std::cout << "Core Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";


//     // coreNodeP* root = tree.getNode();

//     // while(root != nullptr && root->k != 0){

//     //         coreNodeP* rightChild = new coreNodeP();
//     //         root->right = rightChild;
//     //         PathBylmd(mg, root->degs, root->k, root->lmd+1, rightChild, root->valid, root->invalid, root->cnts, root->e); 

//     //         root = root->left;
//     // }

// }

// void FCTreeBuilderCoreParallelByK::Execute(MultilayerGraph &mg, FCCoreTree &tree){

//     coreNodeP* node = tree.getNode();
        
//     uint count = 0;
//     uint n_vertex = mg.GetN(); // number of vertex
//     uint n_layers = mg.getLayerNumber();
//     uint **degs, **adj_list;
//     bool* valid = new bool[n_vertex]; // 1 is valid
//     uint* cnts = new uint[n_vertex];
//     uint* invalid = new uint[n_vertex];

//     degs = new uint*[n_vertex];


//     // Parallel init the degree and valid part
//     #pragma omp parallel
//     {
//         #pragma omp for schedule(static)
//         for(int v = 0; v < n_vertex; v ++){
//              degs[v] = new uint[n_layers];
//              valid[v] = true; // 1 is valid
//         } 

//         #pragma omp for schedule(static) collapse(2)
//         for(int v = 0; v < n_vertex; v ++){
//             // degs[v] = new uint[n_layers];
//             for(int l = 0; l < n_layers; l ++){
//                 degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
//             }
//         }
//     }



//     uint klmd[2];
//     klmd[0] = 1; // k
//     klmd[1] = 1; // lmds

//     BuildSubFCTree(tree, mg, degs, klmd, node, valid, invalid, cnts);
    
//     // Free the memory
//     for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
//     delete[] degs;

// }



// // Mix strategy

// void FCTreeBuilderCoreParallelByK::PathByLmdTask(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts, uint e){
//     uint n_vertex = mg.GetN(); // number of vertex
//     uint n_layers = mg.getLayerNumber(); // number of layer
    

// #pragma omp parallel
// {

//     #pragma omp single
//     {
//         uint new_e = PeelInvalidInParallelByCount(mg, degs, k, lmd, node, valid, invalid, cnts, e, false);
//         // means the (k, lambda)-constaint has the valid vertex
//         if(n_vertex - new_e > 0){
//             lmd += 1;
//             coreNodeP* rightChild = new coreNodeP();
//             node->right = rightChild;
//             PathBylmd(mg, degs, k, lmd, rightChild, valid, invalid, cnts, new_e);
//         }
//     }
// }

// } 

// void FCTreeBuilderCoreParallelByK::BuildSubFCTreeMix(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, uint* invalid, uint* cnts){

//     uint k = klmd[0];
//     uint lmd = klmd[1];
//     uint n_layers = mg.getLayerNumber();

//     uint n_vertex = mg.GetN(); // number of vertex

//     auto start_time_serial = omp_get_wtime(); 
//     PathSerial(mg, degs, k, lmd, node, valid, invalid, cnts, 0);
//     auto end_time_serial = omp_get_wtime(); 

//     double elapsed_time_serial = end_time_serial - start_time_serial;
//     std::cout << "Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";

// #pragma omp parallel
// {
//     #pragma omp single
//     {
//         coreNodeP* root = tree.getNode();
//         while(root != nullptr && root->k != 0){
//                 #pragma omp task
//                 {
//                     coreNodeP* rightChild = new coreNodeP();
//                     root->right = rightChild;
//                     PathByLmdTask(mg, root->degs, root->k, root->lmd+1, rightChild, root->valid, root->invalid, root->cnts, root->e); 
//                 }

//                 root = root->left;
//         }
//         #pragma omp taskwait
//     }
// }

// }



// void FCTreeBuilderCoreParallelByK::ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree){
    
//     coreNodeP* node = tree.getNode();
        
//     uint count = 0;
//     uint n_vertex = mg.GetN(); // number of vertex
//     uint n_layers = mg.getLayerNumber();
//     uint  **adj_list;
//     int **degs;
//     bool* valid = new bool[n_vertex]; // 1 is valid
//     uint* cnts = new uint[n_vertex];
//     uint* invalid = new uint[n_vertex];

//     degs = new int*[n_vertex];


//     // Parallel init the degree and valid part
//     #pragma omp parallel
//     {
//         #pragma omp for schedule(static)
//         for(int v = 0; v < n_vertex; v ++){
//              degs[v] = new int[n_layers];
//              valid[v] = true; // 1 is valid
//         } 

//         #pragma omp for schedule(static) collapse(2)
//         for(int v = 0; v < n_vertex; v ++){
//             // degs[v] = new uint[n_layers];
//             for(int l = 0; l < n_layers; l ++){
//                 degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
//             }
//         }
//     }



//     uint klmd[2];
//     klmd[0] = 1; // k
//     klmd[1] = 1; // lmds

//     BuildSubFCTreeMix(tree, mg, degs, klmd, node, valid, invalid, cnts);
    
//     // Free the memory
//     for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
//     delete[] degs; 
// }