#include "FCSyncLeft.h"

FCSyncLeft::FCSyncLeft(/* args */){
}

FCSyncLeft::~FCSyncLeft(){
}


void FCSyncLeft::constructCoreSync(coreNodeP *node, uint k, uint lmd, uint n_vertex, uint n_layer, bool* valid, int** degs, int total, bool serial){
        
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
        node->degs = new int*[n_vertex];
        for(uint v = 0; v < n_vertex; v ++){
            node->degs[v] = new int[n_layer];
            memcpy(node->degs[v], degs[v], n_layer * sizeof(int));
        }
    }

}

void FCSyncLeft::PeelSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); 
    int *cnts = new int[n_vertex];
    // memset(cnts, 0, sizeof(int)*n_vertex);
    
    double sum_core = 0.0;
    double t_par0 = omp_get_wtime();
    int P = 0;


    #pragma omp parallel shared(cnts, valid, degs)  reduction(+:sum_core) num_threads(32)
    {

        P = omp_get_num_threads();
        double local_core = 0.0;
        uint **adj_lst;
        int buff_size = n_vertex;
        int *buff = (int *)malloc(buff_size*sizeof(int));

        int start = 0, end = 0;
        int cnt = 0;
        int chunk_size = n_vertex/10;

        double t0 = omp_get_wtime();
        #pragma omp for schedule(dynamic, 512)
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
        double t1 = omp_get_wtime();
        local_core += (t1 - t0);

        #pragma omp barrier
        // printf("Thread %d removed node end = %d\n", omp_get_thread_num(), end);

        double t2 = omp_get_wtime();
        while(start < end){
            int vv = buff[start];
            start ++;
            for(uint l = 0; l < n_layers; l ++){
                adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint i = 1; i <= adj_lst[vv][0]; i ++){
                    uint u = adj_lst[vv][i]; // the neighbourhood
                    // if(valid[u] == 0) continue; // only process if u is valid
                    auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
                    if(originDeg == k){
                        auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
                        if(originCnt == lmd && __sync_bool_compare_and_swap(&valid[u], 1, 0)){
                            buff[end++] = u;
                       }
                    }
                }
            }
        }
        __sync_fetch_and_add(&total, end);
        free(buff);
        double t3 = omp_get_wtime();
        local_core += (t3 - t2);
        sum_core += local_core;

    }
    // cout << "P here is " << P << endl;
    double t_par1 = omp_get_wtime();
    double T_par_wall = t_par1 - t_par0;    
    
    stats.total_core     += sum_core;
    stats.total_capacity += (double)P * T_par_wall;
    stats.calls++;


    if(total < n_vertex){
        constructCoreSync(node, k, lmd, n_vertex, n_layers, valid, degs, total, serial);
    }
    // cout << "total = " << total << endl;
    delete[] cnts;
   
}


void FCSyncLeft::PathSerialSync(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats){

    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    
    PeelSync(mg, degs, k, lmd, node, valid, serial, total, stats);
 
    // means the (k, lambda)-constaint has the valid vertex
     if(node->length > 0){
         lmd += 1;
         if(lmd <= n_layers){
             coreNodeP* rightChild = new coreNodeP();
             node->right = rightChild;
             PathSerialSync(mg, degs, k, lmd, rightChild, valid, serial, total, stats);
         }else{
             node->right = nullptr;
         }
     }
 }
 

void FCSyncLeft::PathByK(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int& total, UtilStats& stats){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    

    PeelSync(mg, degs, k, lmd, node, valid, false, total, stats);

    // means the (k, lambda)-constaint has the valid vertex
    if(node->length > 0){
        k += 1;
        coreNodeP* leftChild = new coreNodeP();
        node->left = leftChild;
        PathByK(mg, degs, k, lmd, leftChild, node->valid, total, stats);
    }

} 

void FCSyncLeft:: BuildSubFCTreeSync(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, int& total, UtilStats& stats){
    
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint n_layers = mg.getLayerNumber();

    uint n_vertex = mg.GetN(); // number of vertex

    auto start_time_serial = omp_get_wtime(); 
    PathSerialSync(mg, degs, k, lmd, node, valid, 1, total, stats);
    auto end_time_serial = omp_get_wtime(); 

    double elapsed_time_serial = end_time_serial - start_time_serial;
    std::cout << "Core Parallel Serial part Elapsed time: " << elapsed_time_serial << " seconds\n";


    coreNodeP* root = tree.getNode();


    while(root != nullptr && root->k != 0){

        coreNodeP* leftChild = new coreNodeP();
        root->left = leftChild;
        PathByK(mg, root->degs, root->k+1, root->lmd, leftChild, root->valid, root->total, stats); 

        root = root->right;
    }

}

void FCSyncLeft::Execute(MultilayerGraph &mg, FCCoreTree &tree){
  
    UtilStats peel_util;

    coreNodeP* node = tree.getNode();
    int count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint  **adj_list;
    int **degs;
    bool* valid = new bool[n_vertex]; // 1 is valid
    int total = 0;

    degs = new int*[n_vertex];

        // Parallel init the degree and valid part
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
                degs[v] = new int[n_layers];
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

    BuildSubFCTreeSync(tree, mg, degs, klmd, node, valid, total, peel_util);
    
    double U_overall = peel_util.total_core / peel_util.total_capacity;

    cout << "U_overall = " << U_overall*100 << " % " << endl;

     // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;


}


// ========== The following are mix strategy==========


void FCSyncLeft::PeelSyncMix(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, bool serial, int &total, UtilStats& stats){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); 
    int *cnts = new int[n_vertex];
    memset(cnts, 0, sizeof(int)*n_vertex);
    
    double sum_core = 0.0;
    double t_par0 = omp_get_wtime();
    int P = 0;
    
    // printf("Hello I am here\n");

    #pragma omp parallel shared(cnts, valid, degs) num_threads(8)
    {
        P = omp_get_num_threads();
        double local_core = 0.0;

        uint **adj_lst;
        int buff_size = n_vertex;
        int *buff = (int *)malloc(buff_size*sizeof(int));

        int start = 0, end = 0;
        int cnt = 0;
        int chunk_size = (n_vertex) / (80);
        
        double t0 = omp_get_wtime();
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
        double t1 = omp_get_wtime();
        local_core += (t1 - t0);

        #pragma omp barrier
        // printf("Thread %d removed node end = %d\n", omp_get_thread_num(), end);
        double t2 = omp_get_wtime();
        while(start < end){
            int vv = buff[start];
            start ++;
            for(uint l = 0; l < n_layers; l ++){
                adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint i = 1; i <= adj_lst[vv][0]; i ++){
                    uint u = adj_lst[vv][i]; // the neighbourhood
                    // if(valid[u] == 0) continue; // only process if u is valid
                    //  // minus one and return the old value
                    auto originDeg = __sync_fetch_and_sub(&degs[u][l], 1);
                    if(originDeg == k){
                        auto originCnt = __sync_fetch_and_sub(&cnts[u], 1);
                        if(originCnt == lmd && __sync_bool_compare_and_swap(&valid[u], 1, 0)){
                            buff[end++] = u;
                       }
                    }
                }
            }
        }
        __sync_fetch_and_add(&total, end);
        free(buff);

        double t3 = omp_get_wtime();
        local_core += (t3 - t2);
        sum_core += local_core;
    }
    
    double t_par1 = omp_get_wtime();
    double T_par_wall = t_par1 - t_par0;    

    // cout << "P = " << P << endl;
    
    stats.total_core     += sum_core;
    stats.total_capacity += (double)P * T_par_wall;
    stats.calls++;

    if(total < n_vertex){
        constructCoreSync(node, k, lmd, n_vertex, n_layers, valid, degs, total, serial);
    }

    delete[] cnts;
   
}

void FCSyncLeft::PathByKMix(MultilayerGraph &mg, int **degs, uint k, uint lmd, coreNodeP* node, bool* valid, int& total, UtilStats& stats){
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber(); // number of layer
    

    PeelSyncMix(mg, degs, k, lmd, node, valid, false, total, stats);

    // means the (k, lambda)-constaint has the valid vertex
    if(node->length > 0){
        k += 1;
        coreNodeP* leftChild = new coreNodeP();
        node->left = leftChild;
        PathByKMix(mg, degs, k, lmd, leftChild, node->valid, total, stats);
    }

} 

void FCSyncLeft::BuildSubFCTreeSyncMix(FCCoreTree &tree, MultilayerGraph &mg, int **degs, uint *klmd, coreNodeP* node, bool* valid, int& total, UtilStats& stats){
    
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint n_layers = mg.getLayerNumber();

    uint n_vertex = mg.GetN(); // number of vertex

    auto start_time_serial = omp_get_wtime(); 
    PathSerialSync(mg, degs, k, lmd, node, valid, 1, total, stats);
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
                    PathByKMix(mg, thisNode->degs, thisNode->k+1, thisNode->lmd, leftChild, thisNode->valid, thisNode->total, stats); 
                }
                root = root->right;
            }
        }
        #pragma omp taskwait
    }


}

void FCSyncLeft::ExecuteMix(MultilayerGraph &mg, FCCoreTree &tree){

    UtilStats peel_util;

    coreNodeP* node = tree.getNode();
    int count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **adj_list;
    int **degs;
    bool* valid = new bool[n_vertex]; // 1 is valid
    int total = 0;

    degs = new int*[n_vertex];

        // Parallel init the degree and valid part
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
                degs[v] = new int[n_layers];
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

    BuildSubFCTreeSyncMix(tree, mg, degs, klmd, node, valid, total, peel_util);
    
    double U_overall = peel_util.total_core / peel_util.total_capacity;
    cout << "U_overall = " << U_overall*100 << " % " << endl;

     // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;


}