#include "CoreParallelNew.h"


CoreParallelNew::~CoreParallelNew(){

}

void CoreParallelNew::constructCore(coreNodeP *node, uint k, uint lmd,  uint n_vertex, uint n_layer, uint* valid, int** degs, uint* cnts, int invalidIndex, bool serial){
    
    node->k = k;
    node->lmd = lmd;
    node->valid2 = new uint[n_vertex];
    node->invalidIndex = invalidIndex;
     

    // Plan A store the valid array
    memcpy(node->valid2, valid, sizeof(uint) * n_vertex);

    // int numCnt = 0;

    // for(int i = 0; i < n_vertex; i ++){
    //     if(node->valid2[i] == 0){
    //         numCnt ++;
    //     }
    // }
    // if(lmd == 10 && k == 3){
    //     cout << "k = " << k << " lmd = " << lmd << " length = " << numCnt << endl;
    // }

    // If serial store the deg info and the invalid info
    if(serial){
        node->degs = new int*[n_vertex];
        for(uint v = 0; v < n_vertex; v ++){
            node->degs[v] = new int[n_layer];
            memcpy(node->degs[v], degs[v], n_layer * sizeof(uint));
        }

        node->cnts = new uint[n_vertex];
        memcpy(node->cnts, cnts, n_vertex * sizeof(uint));

    }

}

bool CoreParallelNew::PeelInvalidInParallel(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, uint* valid, uint* cnts, int invalidIndex, bool serial){
    uint n_layers = mg.getLayerNumber(); // number of layer
    uint n_vertex = mg.GetN(); // number of vertex
    uint **adj_lst;
    uint newInvalidIndex;
    uint cnt;
    bool flag = false;

#pragma omp parallel private(cnt) shared(cnts, flag)
{
    uint tid = omp_get_thread_num();
    uint maxthreads = omp_get_num_threads();

    for(int v = tid; v < n_vertex; v += maxthreads){
        cnt = 0;
        if(valid[v] != 0){
            cnts[v] = 0;
            continue; // only process the valid vertex
        } 
        for(int l = 0; l < n_layers; l ++){
            if(degs[v][l] >= k){
                cnt += 1;
            }
        }

        cnts[v] = cnt;


        if(cnts[v] < lmd){
            valid[v] = invalidIndex;
            flag = true;
        }
    }
}

    // while flag = true
    while(flag){
        flag = false;
        newInvalidIndex = invalidIndex + 1; 
        
        #pragma omp parallel private(adj_lst) shared(valid, degs, cnts, mg)
        {
                #pragma omp for
                for(int v = 0; v < n_vertex; v ++){

                    if(valid[v] == invalidIndex){

                        for(int l = 0; l < n_layers; l ++){
                            adj_lst = mg.GetGraph(l).GetAdjLst();
                            for(uint i = 1; i <= adj_lst[v][0]; i ++){
                                uint u = adj_lst[v][i];
                                if(valid[u] != 0){
                                    continue;
                                }
                                auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
                                if(originDeg == k){
                                    auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
                                    if(originCnt == lmd && valid[u] == 0){
                                        flag = true;
                                        valid[u] = newInvalidIndex; 
                                        cnts[u] = 0;
                                    } 
                                }
                            }
                        }

                    }

                }
                #pragma omp barrier
        }
        invalidIndex = newInvalidIndex; 
    }

bool isContinue = false;
#pragma omp parallel
{ 
    #pragma omp for
    for(int v = 0; v < n_vertex; v ++){
        if(valid[v] == 0){
           isContinue = true;
            #pragma omp cancel for 
        }
    }
    
}

    if(isContinue){
        constructCore(node, k, lmd, n_vertex, n_layers, valid, degs, cnts, invalidIndex, serial);
    }else{
        node->k = 0;
        node->lmd = 0;
    }

    return isContinue;

}

void CoreParallelNew::PathSerial(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, uint* valid, uint* cnts, int invalidIndex){

    uint n_layers = mg.getLayerNumber(); // number of layer
    uint n_vertex = mg.GetN(); // number of vertex

    if(k == 1 && lmd == 1){
        invalidIndex = 1;
    }
    
    bool isContinue = PeelInvalidInParallel(mg, degs, k, lmd, node, valid, cnts, invalidIndex, true);

    if(isContinue){
        lmd += 1;
        if(lmd <= n_layers){
            coreNodeP* rightChild = new coreNodeP();
            node->right = rightChild;
            PathSerial(mg, degs, k, lmd, rightChild, valid, cnts, node->invalidIndex);
        }else{
            node->right = nullptr;
        } 
    }else{
        node->right = nullptr;
    }

}

void CoreParallelNew::PathLeft(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, uint* valid, uint* cnts, int invalidIndex){
    uint n_layers = mg.getLayerNumber(); // number of layer
    uint n_vertex = mg.GetN(); // number of vertex

    bool isContinue = PeelInvalidInParallel(mg, degs, k, lmd, node, valid, cnts, invalidIndex, false);
    
    if(isContinue){
        k += 1;
        // if(k <= 10){
            coreNodeP* leftChild = new coreNodeP();
            node->left = leftChild;
            PathLeft(mg, degs, k, lmd, leftChild, valid, cnts, node->invalidIndex);
        // }
    }else{
        node->left = nullptr;
    } 
} 

void CoreParallelNew::BuildSubFCTree(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, uint* valid, uint* cnts){
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint n_layers = mg.getLayerNumber();

    uint n_vertex = mg.GetN(); // number of vertex

    coreNodeP* root = tree.getNode();
    root = new coreNodeP();
    root->invalidIndex = 1;

    auto start_time_serial = omp_get_wtime(); 
    PathSerial(mg, degs, k, lmd, root, valid, cnts, root->invalidIndex);
    auto end_time_serial = omp_get_wtime(); 
    double elapsed_time_serial = end_time_serial - start_time_serial;
    std::cout << "Core Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";


    
    // cout << root->right << endl;
    // while(root != nullptr && root->k != 0){
        coreNodeP* leftChild = new coreNodeP(); 
        root->left = leftChild;
        PathLeft(mg, root->degs, root->k+1, root->lmd, leftChild, root->valid2, root->cnts, root->invalidIndex);
        root = root->right;
    // }
}


void CoreParallelNew::Execute(MultilayerGraph &mg, FCCoreTree &tree){

    coreNodeP* node = tree.getNode();

    uint count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **adj_list;
    int **degs;
    uint* valid = new uint[n_vertex]; // 1 is valid
    uint* cnts = new uint[n_vertex];

    degs = new int*[n_vertex];


    // Parallel init the degree and valid part
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
             degs[v] = new int[n_layers];
             valid[v] = 0; // 0 is valid, others are not valid
        } 

        #pragma omp for schedule(static) collapse(2)
        for(int v = 0; v < n_vertex; v ++){
            // degs[v] = new uint[n_layers];
            for(int l = 0; l < n_layers; l ++){
                degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
            }
        }
    }



    uint klmd[2];
    klmd[0] = 1; // k
    klmd[1] = 1; // lmds

    BuildSubFCTree(tree, mg, degs, klmd, valid, cnts);

    // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;

}
