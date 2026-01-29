#include "CoreIndex.h"

CoreIndex::CoreIndex(/* args */){
}

CoreIndex::~CoreIndex(){
}

struct  MaxDensity{
    uint resk;
    uint reslmd;
    float density;
};



struct LayerDegree2{
    uint id;
    float degree;
};

bool cmp(const LayerDegree2 &a, const LayerDegree2 &b){
    return a.degree > b.degree;
}

bool cmp2(const MaxDensity &a, const MaxDensity &b){
    return a.density > b.density;
}


bool reverseSort(uint a, uint b) {
    return a > b; 
}

void print(std::vector<std::vector<uint>> B, uint n_vertex){
    cout << "==========Round==========" << endl;
   for(uint k = 1; k <= n_vertex; k ++){
        cout << "k = " << k << ": ";
        for(uint j = 0; j < B[k].size(); j ++){
            cout << B[k][j] << " ";
        }
        cout << endl;
    }
}

uint findLmdthLargest(uint *deg, uint lmd, uint n_layer){

    uint *deg_sort = new uint[n_layer];
    memcpy(deg_sort, deg, sizeof(uint) * n_layer);


    sort(deg_sort, deg_sort+n_layer, reverseSort);

    int res = deg_sort[lmd-1];
    delete[] deg_sort;

    return res;
}



uint* get_index2(MultilayerGraph &mg, uint **degs, uint lmd, ll_uint *id2vtx){
     

    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN(); 
    uint *i_v = new uint[n_vertex];

    std::vector<std::set<uint>> B(n_vertex+1);



    for(uint v = 0; v < n_vertex; v ++){
        i_v[v] = findLmdthLargest(degs[v], lmd, n_layer); // This is where the problem happens
        B[i_v[v]].insert(v);
    }

    
    for(uint k = 0; k <= n_vertex; k ++){
        while(!B[k].empty()){ 
            auto v = B[k].begin();
            uint vv = *v;
            B[k].erase(vv);
            i_v[vv] = k;
            std::set<uint> N;
            for(uint l = 0; l < n_layer; l ++){
                uint **adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint u = 1; u <= adj_lst[vv][0]; u ++){
                    uint nb_u = adj_lst[vv][u]; // node v's neighbourhood in layer l
                    if(i_v[nb_u] <= k) continue;
                    degs[nb_u][l] = degs[nb_u][l] - 1;
                    if(degs[nb_u][l] + 1 == i_v[nb_u]){
                        N.insert(nb_u);
                    }
                }
            }
            
            for(const uint& u : N){
                B[i_v[u]].erase(u);
                i_v[u] = findLmdthLargest(degs[u], lmd, n_layer); 
                B[std::max(i_v[u], k)].insert(u);
            }

        }

    }

    return i_v;
}


uint* get_index(MultilayerGraph &mg, uint **degs, uint lmd, ll_uint *id2vtx){
     

    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN(); 
    uint *i_v = new uint[n_vertex];
    uint *core = new uint[n_vertex];

    std::vector<std::set<uint>> B(n_vertex+1);

    bool visit[n_vertex];
    memset(visit, false, sizeof(bool) * n_vertex);


    for(uint v = 0; v < n_vertex; v ++){
        i_v[v] = findLmdthLargest(degs[v], lmd, n_layer); // This is where the problem happens
        B[i_v[v]].insert(v);
    }

    
    for(uint k = 0; k <= n_vertex; k ++){
        while(!B[k].empty()){ 
            auto v = B[k].begin();
            uint vv = *v;
            B[k].erase(vv);
            core[vv] = k;
            visit[vv] = true;
            std::set<uint> N;
            for(uint l = 0; l < n_layer; l ++){
                uint **adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint u = 1; u <= adj_lst[vv][0]; u ++){
                    uint nb_u = adj_lst[vv][u]; // node v's neighbourhood in layer l
                    if(i_v[nb_u] <= k) continue;
                    if(visit[nb_u]) continue;
                    degs[nb_u][l] = degs[nb_u][l] - 1;
                    if(degs[nb_u][l] + 1 == i_v[nb_u]){
                        N.insert(nb_u);
                    }
                }
            }
            
            for(const uint& u : N){
                B[i_v[u]].erase(u);
                i_v[u] = findLmdthLargest(degs[u], lmd, n_layer); 
                B[std::max(i_v[u], k)].insert(u);
            }

        }

    }

    return core;
}



void CoreIndex::Execute(MultilayerGraph &mg, ll_uint *id2vtx){
    cout << "I am the CoreIndex" << endl;
    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN();

    cout << "n_layer = " << n_layer << endl;
    cout << "n_vertex = " << n_vertex << endl;

    uint **degs;
    degs = new uint*[n_vertex];
   
    uint** final_res = new uint*[n_layer+1];

    for(int i = 0; i <= n_layer; i ++){
        final_res[i] = new uint[n_vertex];
    }


    // int cnt = 0;
    for(uint lmd = 1; lmd <= n_layer; lmd ++){

        for(int v = 0; v < n_vertex; v ++){
             degs[v] = new uint[n_layer];
        } 

         for(int v = 0; v < n_vertex; v ++){
            for(int l = 0; l < n_layer; l ++){
                degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
            }
        }

        uint *core = get_index2(mg, degs, lmd, id2vtx);
        memcpy(final_res[lmd], core, n_vertex * sizeof(uint));

        // if(lmd == 3){
        //     for(int vv = 0; vv < n_vertex; vv ++){
        //         cnt += final_res[lmd][vv] >= 4;
        //     }
        // }
    
        // Free core array
        delete[] core;
    }

    // printf("cnt = %d\n", cnt);
}



void CoreIndex::ExecuteParallel(MultilayerGraph &mg, ll_uint *id2vtx){
    
    UtilStats peel_util;

    cout << "I am the ParCoreIndex" << endl;
    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN();

    cout << "n_layer = " << n_layer << endl;
    cout << "n_vertex = " << n_vertex << endl;

    uint **degs;
    degs = new uint*[n_vertex];
   
    uint** final_res = new uint*[n_layer+1];

    for(int i = 0; i <= n_layer; i ++){
        final_res[i] = new uint[n_vertex];
    }

    for(int v = 0; v < n_vertex; v ++){
        degs[v] = new uint[n_layer];
    } 

    for(int v = 0; v < n_vertex; v ++){
        for(int l = 0; l < n_layer; l ++){
            degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
        }
    }

    double sum_core = 0.0;
    double t_par0 = omp_get_wtime();
    int P = 0;

    // int cnt = 0;
    #pragma omp parallel for reduction(+:sum_core)
    for(uint lmd = 1; lmd <= n_layer; lmd ++){

        P = omp_get_num_threads();
        double t0 = omp_get_wtime();
        double local_core = 0.0;
        // Copy the degs degree
        uint **degs_copy = new uint*[n_vertex];
        for(uint v = 0; v < n_vertex; v ++){
            degs_copy[v] = new uint[n_layer];
            memcpy(degs_copy[v], degs[v], n_layer * sizeof(uint));
        }
        // for(int v = 0; v < n_vertex; v ++){
        //      degs[v] = new uint[n_layer];
        // } 

        //  for(int v = 0; v < n_vertex; v ++){
        //     for(int l = 0; l < n_layer; l ++){
        //         degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
        //     }
        // }

        uint *core = get_index2(mg, degs_copy, lmd, id2vtx);
        memcpy(final_res[lmd], core, n_vertex * sizeof(uint));

        // if(lmd == 3){
        //     for(int vv = 0; vv < n_vertex; vv ++){
        //         cnt += final_res[lmd][vv] >= 4;
        //     }
        // }
    
        // Free core array
        delete[] core;
        double t1 = omp_get_wtime();
        local_core += (t1 - t0);
        sum_core += local_core;
    }
    
    double t_par1 = omp_get_wtime();
    double T_par_wall = t_par1 - t_par0;

    peel_util.total_core     += sum_core;
    peel_util.total_capacity += (double)P * T_par_wall;
    peel_util.calls++;

    cout << "P = " << P << endl;
    double U_overall = peel_util.total_core / peel_util.total_capacity;
    cout << "U_overall = " << U_overall*100 << " % " << endl;

    // printf("cnt = %d\n", cnt);
}


// The following two algorithms are getting weight denset subgraph algorithm

void calculateWeightDenestSubgraphCommon(uint k, uint lmd, MultilayerGraph &mg, std::unordered_set<uint> valid_node_set, MaxDensity &maximum_density, float beta, uint &res_k, uint &res_lmd){
    
    uint n_layers = mg.getLayerNumber();
    uint n_vertex = mg.GetN();
    uint n_node = valid_node_set.size();

    LayerDegree2 *layer_degree = new LayerDegree2[n_layers];

    for(uint l = 0; l < n_layers; l ++){
        uint **adj = mg.GetGraph(l).GetAdjLst();
        float n_edge = 0;
        for (const auto& vv:valid_node_set){
            for(uint u = 1; u <= adj[vv][0]; u ++){
                auto uu = adj[vv][u];
                if(valid_node_set.find(uu) != valid_node_set.end()){
                    n_edge ++;
                }
            }
        }
        layer_degree[l].id = l;
        layer_degree[l].degree = n_edge/2;
    }

    float *maximum_average_degree = new float[n_layers];
    sort(layer_degree, layer_degree+n_layers, cmp);

    for(uint l = 0; l < n_layers; l ++){
        maximum_average_degree[l] = layer_degree[l].degree / n_node;
        float layer_density = maximum_average_degree[l] * pow((l+1), beta); 
        if(layer_density >= maximum_density.density){
            maximum_density.density = layer_density;
            maximum_density.resk = k;
            maximum_density.reslmd = lmd;
        }
    }
    
    delete[] layer_degree;

}

uint* get_wds(MultilayerGraph &mg, uint **degs, uint lmd, ll_uint *id2vtx, MaxDensity &maximum_density, float beta, uint &res_k, uint &res_lmd){
     

    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN(); 
    uint *i_v = new uint[n_vertex];
    uint *core = new uint[n_vertex];

    std::vector<std::set<uint>> B(n_vertex+1);

    bool visit[n_vertex];
    memset(visit, false, sizeof(bool) * n_vertex);
    
    // valid set
    std::unordered_set<uint> valid_node_set;

    for(uint v = 0; v < n_vertex; v ++){
        i_v[v] = findLmdthLargest(degs[v], lmd, n_layer); // This is where the problem happens
        B[i_v[v]].insert(v);
        valid_node_set.insert(v);
    }



    for(uint k = 0; k <= n_vertex; k ++){
        std::unordered_set<uint> t_node_set;
        while(!B[k].empty()){ 
            auto v = B[k].begin();
            uint vv = *v;
            B[k].erase(vv);
            core[vv] = k;
            visit[vv] = true;
            t_node_set.insert(vv);
            std::set<uint> N;
            for(uint l = 0; l < n_layer; l ++){
                uint **adj_lst = mg.GetGraph(l).GetAdjLst();
                for(uint u = 1; u <= adj_lst[vv][0]; u ++){
                    uint nb_u = adj_lst[vv][u]; // node v's neighbourhood in layer l
                    if(i_v[nb_u] <= k) continue;
                    if(visit[nb_u]) continue;
                    degs[nb_u][l] = degs[nb_u][l] - 1;
                    if(degs[nb_u][l] + 1 == i_v[nb_u]){
                        N.insert(nb_u);
                    }
                }
            }
            
            for(const uint& u : N){
                B[i_v[u]].erase(u);
                i_v[u] = findLmdthLargest(degs[u], lmd, n_layer); 
                B[std::max(i_v[u], k)].insert(u);
            }

        }

        if(k >= 1 && !valid_node_set.empty()){
            calculateWeightDenestSubgraphCommon(k, lmd, mg, valid_node_set, maximum_density, beta, res_k, res_lmd);
        }

        for(const auto& ele : t_node_set){
           valid_node_set.erase(ele); 
        }

    }

    return core;
}

void CoreIndex::WdsCoreIndex(MultilayerGraph &mg, ll_uint *id2vtx, float beta){
    cout << "I am the WDS CoreIndex" << endl;
    uint n_layer = mg.getLayerNumber();
    uint n_vertex = mg.GetN();

    cout << "n_layer = " << n_layer << endl;
    cout << "n_vertex = " << n_vertex << endl;

    uint **degs;
    degs = new uint*[n_vertex];

    // degs_sort = new uint*[n_vertex];
    
    // Parallel init the degree information
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for(int v = 0; v < n_vertex; v ++){
             degs[v] = new uint[n_layer];
        } 

        #pragma omp for schedule(static) collapse(2)
        for(int v = 0; v < n_vertex; v ++){
            for(int l = 0; l < n_layer; l ++){
                degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
            }
        }
    }

    // float maximum_density = 0.0f;
    uint res_k = 0;
    uint res_lmd = 0;

    MaxDensity* maximum_density = new MaxDensity[n_layer+1];
    


    // for(uint lmd = 1; lmd <= n_layer; lmd ++){

        
    // }

    #pragma omp parallel for
    for(uint lmd = 1; lmd <= n_layer; lmd ++){

        maximum_density[lmd].reslmd = 0;
        maximum_density[lmd].resk = 0;
        maximum_density[lmd].density = 0;

        // Copy the degs degree
        uint **degs_copy = new uint*[n_vertex];
        for(uint v = 0; v < n_vertex; v ++){
            degs_copy[v] = new uint[n_layer];
            memcpy(degs_copy[v], degs[v], n_layer * sizeof(uint));
        }

        get_wds(mg, degs_copy, lmd, id2vtx, maximum_density[lmd], beta, res_k, res_lmd);
    }

    sort(maximum_density, maximum_density + n_layer, cmp2);

    cout << "maximum_density = " << maximum_density[0].density << endl;
    cout << "res_k = " << maximum_density[0].resk << endl;
    cout << "res_lmd = " << maximum_density[0].reslmd << endl;
     
}
