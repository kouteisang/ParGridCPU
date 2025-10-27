#include "FCSyncLeft.h"

FCSyncLeft::FCSyncLeft(/* args */){
}

FCSyncLeft::~FCSyncLeft(){
}


void FCSyncLeft::constructCoreSync(coreNodeP *node, uint k, uint lmd, uint n_vertex, uint n_layer, bool* valid, uint** degs, int total, bool serial){
        
    node->k = k;
    node->lmd = lmd;
    node->valid = new bool[n_vertex];
    node->total = total;
    node->length = n_vertex - total;
     
    // Plan A store the valid array
    memcpy(node->valid, valid, sizeof(bool) * n_vertex);

    // int count = 0;
    // for(int i = 0; i < n_vertex; i ++){
    //     if(node->valid[i] == 1){
    //         count ++;
    //     }
    // }
   
    // cout << "node->k = " << node->k << " node->lmd = " << node->lmd << " node->length = " << node->length << endl;

    if(serial){
        node->degs = new uint*[n_vertex];
        for(uint v = 0; v < n_vertex; v ++){
            node->degs[v] = new uint[n_layer];
            memcpy(node->degs[v], degs[v], n_layer * sizeof(uint));
        }
    }

}

void FCSyncLeft::PeelSync(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); 
    int *cnts = new int[n_vertex];
    // memset(cnts, 0, sizeof(int)*n_vertex);
    
    
    #pragma omp parallel shared(cnts, valid, degs)  num_threads(40)
    {

        uint **adj_lst;
        int buff_size = n_vertex;
        int *buff = (int *)malloc(buff_size*sizeof(int));

        int start = 0, end = 0;
        int cnt = 0;
        int chunk_size = n_vertex/10;

        #pragma omp for schedule(dynamic, 512)
        // #pragma omp for schedule(guided, 2048) nowait
        for(int v = 0; v < n_vertex; v ++){
            cnt = 0;
            if(valid[v] == 0){
                 cnts[v] = 0;
                 continue; // only process the valid vertex
            } 
            for(int l = 0; l < n_layers; l ++){
                cnt += (degs[v][l] >= k); 
            }

            if(cnt < lmd){
                cnts[v] = 0;
                valid[v] = 0;
                buff[end ++] = v;
            }else{
                cnts[v] = cnt;   
            }
        }

        #pragma omp barrier
        // printf("Thread %d removed node end = %d\n", omp_get_thread_num(), end);

        while(start < end){
            int vv = buff[start];
            start ++;
            for(uint l = 0; l < n_layers; l ++){
                adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint i = 1; i <= adj_lst[vv][0]; i ++){
                    uint u = adj_lst[vv][i]; // the neighbourhood
                    if(valid[u] == 0) continue; // only process if u is valid
                    //  // minus one and return the old value
                    auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
                    if(originDeg == k){
                        auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
                        // valid[u] == 1){
                        if(originCnt == lmd && __sync_bool_compare_and_swap(&valid[u], 1, 0)){
                            cnts[u] = 0;
                            buff[end++] = u;
                       }
                    }
                }
            }
        }
        __sync_fetch_and_add(&total, end);
        free(buff);
    }

    if(total < n_vertex){
        constructCoreSync(node, k, lmd, n_vertex, n_layers, valid, degs, total, serial);
    }
    // cout << "total = " << total << endl;
    delete[] cnts;
   
}


void FCSyncLeft::PathSerialSync(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total){

    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    
    PeelSync(mg, degs, k, lmd, node, valid, serial, total);
 
    // means the (k, lambda)-constaint has the valid vertex
     if(node->length > 0){
         lmd += 1;
         if(lmd <= n_layers){
             coreNodeP* rightChild = new coreNodeP();
             node->right = rightChild;
             PathSerialSync(mg, degs, k, lmd, rightChild, valid, serial, total);
         }else{
             node->right = nullptr;
         }
     }
 }
 

void FCSyncLeft::PathByK(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int& total){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    

    PeelSync(mg, degs, k, lmd, node, valid, false, total);

    // means the (k, lambda)-constaint has the valid vertex
    if(node->length > 0){
        k += 1;
        coreNodeP* leftChild = new coreNodeP();
        node->left = leftChild;
        PathByK(mg, degs, k, lmd, leftChild, node->valid, total);
    }

} 

void FCSyncLeft:: BuildSubFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, bool* valid, int& total){
    
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint n_layers = mg.getLayerNumber();

    uint n_vertex = mg.GetN(); // number of vertex

    auto start_time_serial = omp_get_wtime(); 
    PathSerialSync(mg, degs, k, lmd, node, valid, 1, total);
    auto end_time_serial = omp_get_wtime(); 

    double elapsed_time_serial = end_time_serial - start_time_serial;
    std::cout << "Core Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";


    coreNodeP* root = tree.getNode();


    while(root != nullptr && root->k != 0){

        coreNodeP* leftChild = new coreNodeP();
        root->left = leftChild;
        PathByK(mg, root->degs, root->k+1, root->lmd, leftChild, root->valid, root->total); 

        root = root->right;
    }

}

void FCSyncLeft::Execute(MultilayerGraph &mg, FCCoreTree &tree){
    coreNodeP* node = tree.getNode();
    int count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **degs, **adj_list;
    bool* valid = new bool[n_vertex]; // 1 is valid
    int total = 0;

    degs = new uint*[n_vertex];

        // Parallel init the degree and valid part
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
                degs[v] = new uint[n_layers];
            //  valid[v] = true; // 1 is valid
        } 

        #pragma omp for schedule(static) collapse(2)
        for(int v = 0; v < n_vertex; v ++){
            // degs[v] = new uint[n_layers];
            for(int l = 0; l < n_layers; l ++){
                degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
            }
        }
    }

    memset(valid, true, sizeof(bool) * n_vertex);

    uint klmd[2];
    klmd[0] = 1; // k
    klmd[1] = 1; // lmds

    BuildSubFCTreeSync(tree, mg, degs, klmd, node, valid, total);
    
     // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;


}


// ========== The following are mix strategy==========


void FCSyncLeft::PeelSyncMix(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); 
    int *cnts = new int[n_vertex];
    memset(cnts, 0, sizeof(int)*n_vertex);
    
    
    #pragma omp parallel shared(cnts, valid, degs) num_threads(3)
    {

        uint **adj_lst;
        int buff_size = n_vertex;
        int *buff = (int *)malloc(buff_size*sizeof(int));

        int start = 0, end = 0;
        int cnt = 0;
        int chunk_size = (n_vertex) / (80);

        #pragma omp for schedule(dynamic, chunk_size)
        for(int v = 0; v < n_vertex; v ++){
            cnt = 0;
            if(valid[v] == 0){
                cnts[v] = 0;
                continue; // only process the valid vertex
            } 
            for(int l = 0; l < n_layers; l ++){
                cnt += (degs[v][l] >= k); 
            }

            if(cnt < lmd){
                cnts[v] = 0;
                valid[v] = 0;
                buff[end ++] = v;
            }else{
                cnts[v] = cnt;   
            }

        }

        #pragma omp barrier
        // printf("Thread %d removed node end = %d\n", omp_get_thread_num(), end);

        while(start < end){
            int vv = buff[start];
            start ++;
            for(uint l = 0; l < n_layers; l ++){
                adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint i = 1; i <= adj_lst[vv][0]; i ++){
                    uint u = adj_lst[vv][i]; // the neighbourhood
                    if(valid[u] == 0) continue; // only process if u is valid
                    //  // minus one and return the old value
                    auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
                    if(originDeg == k){
                        auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
                        if(originCnt == lmd && __sync_bool_compare_and_swap(&valid[u], 1, 0)){
                            cnts[u] = 0;
                            buff[end++] = u;
                       }
                    }
                }
            }
        }
        __sync_fetch_and_add(&total, end);
        free(buff);
    }

    if(total < n_vertex){
        constructCoreSync(node, k, lmd, n_vertex, n_layers, valid, degs, total, serial);
    }

    delete[] cnts;
   
}

void FCSyncLeft::PathByKMix(MultilayerGraph &mg, uint **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int& total){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    

    PeelSyncMix(mg, degs, k, lmd, node, valid, false, total);

    // means the (k, lambda)-constaint has the valid vertex
    if(node->length > 0){
        k += 1;
        coreNodeP* leftChild = new coreNodeP();
        node->left = leftChild;
        PathByKMix(mg, degs, k, lmd, leftChild, node->valid, total);
    }

} 

void FCSyncLeft::BuildSubFCTreeSyncMix(FCCoreTree &tree, MultilayerGraph &mg, uint **degs, uint *klmd, coreNodeP* node, bool* valid, int& total){
    
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint n_layers = mg.getLayerNumber();

    uint n_vertex = mg.GetN(); // number of vertex

    auto start_time_serial = omp_get_wtime(); 
    PathSerialSync(mg, degs, k, lmd, node, valid, 1, total);
    auto end_time_serial = omp_get_wtime(); 

    double elapsed_time_serial = end_time_serial - start_time_serial;
    std::cout << "Core Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";


    coreNodeP* root = tree.getNode();

    #pragma omp parallel
    {
        #pragma omp single
        {
            while(root != nullptr && root->k != 0){
                coreNodeP* thisNode = root; // 👈 临时变量，捕捉当前 root

                #pragma omp task firstprivate(thisNode)
                {
                    coreNodeP* leftChild = new coreNodeP();
                    thisNode->left = leftChild;
                    PathByKMix(mg, thisNode->degs, thisNode->k+1, thisNode->lmd, leftChild, thisNode->valid, thisNode->total); 
                }
                root = root->right;
            }
        }
        #pragma omp taskwait
    }


}

void FCSyncLeft::ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree){

    coreNodeP* node = tree.getNode();
    int count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **degs, **adj_list;
    bool* valid = new bool[n_vertex]; // 1 is valid
    int total = 0;

    degs = new uint*[n_vertex];

        // Parallel init the degree and valid part
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
                degs[v] = new uint[n_layers];
            //  valid[v] = true; // 1 is valid
        } 

        #pragma omp for schedule(static) collapse(2)
        for(int v = 0; v < n_vertex; v ++){
            // degs[v] = new uint[n_layers];
            for(int l = 0; l < n_layers; l ++){
                degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
            }
        }
    }

    memset(valid, true, sizeof(bool) * n_vertex);

    uint klmd[2];
    klmd[0] = 1; // k
    klmd[1] = 1; // lmds

    BuildSubFCTreeSyncMix(tree, mg, degs, klmd, node, valid, total);
    
     // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;


}