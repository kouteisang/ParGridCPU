#include "FCTreeBuilderLeft.h"

FCTreeBuilderLeft::~FCTreeBuilderLeft(){
}

// ========== Common method ==========

bool FCTreeBuilderLeft::check(int **degs, uint u, uint *klmd, uint n_layers){
    uint k = klmd[0];
    uint lmd = klmd[1];
    uint cnt = 0;
    for(int l = 0; l < n_layers; l ++){
        if(degs[u][l] >= k){
            cnt += 1;
            if(cnt >= lmd){
                return true;
            }
        }
    }
    return false;
}

uint FCTreeBuilderLeft::peel(MultilayerGraph &mg, int **degs, uint *klmd, uint *core, uint *pos, uint s, uint e){

    uint n_layers = mg.getLayerNumber();
    uint old_s = s;
    uint old_e = e;
    uint **adj_lst;
    uint v, u;
    uint pos_u;
    bool flag = true;

    while(s < e){
        old_e = e;


        for(uint i = 0; i < n_layers; i ++){
            adj_lst = mg.GetGraph(i).GetAdjLst(); 
            for(uint j = s; j < old_e; j ++){
                v = core[j];
                for(uint l = 1; l <= adj_lst[v][0]; l ++){
                    u = adj_lst[v][l];
                    degs[u][i] --;

                    if(degs[u][i] < klmd[0] && pos[u] >= e){
                        flag = check(degs, u, klmd, n_layers);
                        if(!flag){ // node u does not satisfy the condition
                            pos_u = pos[u];
                            core[pos_u] = core[e];
                            core[e] = u;

                            pos[core[pos_u]] = pos_u;
                            pos[u] = e;
                            e ++;
                        }
                    }
                }
            }
 
        }
        

        s = old_e;
    }

    return e;

}

void FCTreeBuilderLeft::restore(MultilayerGraph &mg, int **degs, uint *core, uint old_e, uint new_e){

    uint n_layers = mg.getLayerNumber();
    uint **adj_lst;
    uint u, v;

    // [old_e, new_e) is the elements that change from active -> inactive during this time.
    for(int i = 0; i < n_layers; i ++){
        adj_lst = mg.GetGraph(i).GetAdjLst();
        for(int j = old_e; j < new_e; j ++){   
            v = core[j];
            for(int l = 1; l <= adj_lst[v][0]; l ++){
                u = adj_lst[v][l];
                degs[u][i] ++;
            }
        }
    } 
}

void FCTreeBuilderLeft::PrintCoreInfor(uint *klmd, uint *core, uint new_e, uint n_vertex, int** degs, uint n_layer, long long& sum_num_of_edge_fc){

    // cout << "( ";
    // cout << "k=" << klmd[0] << ",";
    // cout << "lmd=" << klmd[1];
    // cout << ") len  = ";
    // cout << n_vertex - new_e;
    

    if(klmd[1] == 1){
        long long num_of_edge_fc = 0;
        for(int v = 0; v < n_vertex; v ++){
            for(int l = 0; l < n_layer; l ++){
                num_of_edge_fc += degs[v][l];
            }
        }
        cout << "k = " << klmd[0] << ", ";
        cout << "lmd = " << klmd[1] << ", ";
        cout << "num_of_edge_fc = " << num_of_edge_fc/2 << endl;
        sum_num_of_edge_fc += num_of_edge_fc/2;
    }

    // cout << endl;
}

void FCTreeBuilderLeft::constructCore(uint *klmd, uint *core, uint new_e, uint n_vertex, coreNode *node, coreNode *father){

    // The node infor
 
        
    node->k = klmd[0];
    node->lmd = klmd[1];
    node->length = n_vertex - new_e;
    node->core = new uint[node->length];
    memcpy(node->core, core + new_e, (n_vertex - new_e) * sizeof(uint));
    
    // if(father != nullptr){
    //     delete[] father->core;
    //     father->core = nullptr;
    //     if(father->length > node->length){
    //         cout << "father->length = "<< father->length << endl;
    //         cout << "father->length - node->length = "<< father->length - node->length << endl; 
    //         father->core = new uint[father->length - node->length];
    //         father->length = father->length - node->length; 
    //     }
    // }

    // if (father != nullptr) {
    //     // 先判断是否真的需要重建 core 数组
    //     uint old_len = father->length;
    //     uint new_len = old_len > node->length ? (old_len - node->length) : 0;
    
    //     if (father->core != nullptr) {
    //         delete[] father->core;      // 确保释放原内存
    //         father->core = nullptr;
    //     }
    
    //     if (new_len > 0) {
    //         father->core = new uint[new_len];  // 重新分配更小数组
    //         father->length = new_len;
    
    //         // ⚠️ 如果你打算访问这段内存，一定要初始化它（例如置0）
    //         // memset(father->core, 0, new_len * sizeof(uint));
    //     } else {
    //         father->core = nullptr;
    //         father->length = 0;
    //     }
    // }
}


// ========== The following build the tree structure to actually store the data ==========

void FCTreeBuilderLeft::Execute(MultilayerGraph &mg, FCTree &tree){
    coreNode* node = tree.getNode();

        
    uint count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **adj_list;
    int **degs;

    degs = new int*[n_vertex];

    for(int v = 0; v < n_vertex; v ++){
        degs[v] = new int[n_layers];
        for(int l = 0; l < n_layers; l ++){
            degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
        }
    }

    uint* core = new uint[n_vertex];
    uint* pos = new uint[n_vertex];

    for(int i = 0; i < n_vertex; i ++){
        core[i] = i;
        // vertex i is in the i-th position
        pos[i] = i;
    }

    uint klmd[2];
    klmd[0] = 1; // k
    klmd[1] = 1; // lmd

    long long sum_num_of_edge_fc = 0;
    uint e = 0;
    BuildSubFCTree(mg, degs, core, pos, klmd, e, count, node, nullptr,sum_num_of_edge_fc);

    cout << "sum_num_of_edge_fc = " <<  sum_num_of_edge_fc << endl;
    // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;

    cout << "count = " << count << endl;

   
}

void FCTreeBuilderLeft::BuildSubFCTree(MultilayerGraph &mg, int **degs, uint *core, uint *pos, uint *klmd, uint e, uint &count, coreNode* node, coreNode* father, long long& sum_num_of_edge_fc){

    uint s = e;
    uint old_e = e;
    uint new_e = 0;

    uint n_vertex = mg.GetN();
    uint n_layer = mg.getLayerNumber();
    uint old_pos_v;

    for(int i = s; i < n_vertex; i ++){
        uint v = core[i];
        uint cnt = 0;
        bool flag = false;
        for(int j = 0; j < n_layer; j ++){
            if(degs[v][j] >= klmd[0]){
                cnt ++;
                if(cnt >= klmd[1]){
                    flag = true;
                    break;
                }
            }
        }
        if(!flag){
            old_pos_v = pos[v];
            core[old_pos_v] = core[e];
            core[e] = v;

            pos[core[old_pos_v]] = old_pos_v;
            pos[v] = e ++;

        }
    }

    if(s != e){
        new_e = peel(mg, degs, klmd, core, pos, s, e);
    }else if(s == e){
        new_e = e;
    }

    if(n_vertex - new_e > 0){
        count ++;
        PrintCoreInfor(klmd, core, new_e, n_vertex, degs, mg.getLayerNumber(), sum_num_of_edge_fc); 
        constructCore(klmd, core, new_e, n_vertex, node, father);
    }else{
        node->k = 0;
        node->lmd = 0;
        node->length = 0; 
    }

    if(n_vertex - new_e > 0){

        coreNode* leftChild = new coreNode();
        node->left = leftChild;

        klmd[0] += 1;
        // BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count, leftChild);
        BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count, leftChild, node,sum_num_of_edge_fc);  
        klmd[0] -= 1;
    }

    if(n_vertex - new_e > 0 && klmd[0] == 1){

        coreNode* rightChild = new coreNode();
        node->right = rightChild;

        klmd[1] += 1;
        // BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count, rightChild);
        BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count, rightChild, node, sum_num_of_edge_fc);
        klmd[1] -= 1; 

    }

    restore(mg, degs, core, old_e, new_e); 

    s = old_e;
    e = old_e;   
}

// ========== The following without building the tree structure to actually store the data ==========

void FCTreeBuilderLeft::Execute(MultilayerGraph &mg){
    
    uint count = 0;
    uint n_vertex = mg.GetN(); // number of vertex
    uint n_layers = mg.getLayerNumber();
    uint **adj_list;
    int **degs;

    degs = new int*[n_vertex];

    for(int v = 0; v < n_vertex; v ++){
        degs[v] = new int[n_layers];
        for(int l = 0; l < n_layers; l ++){
            degs[v][l] = mg.GetGraph(l).GetAdjLst()[v][0];
        }
    }

    uint* core = new uint[n_vertex];
    uint* pos = new uint[n_vertex];

    for(int i = 0; i < n_vertex; i ++){
        core[i] = i;
        // vertex i is in the i-th position
        pos[i] = i;
    }

    uint klmd[2];
    klmd[0] = 1; // k
    klmd[1] = 1; // lmd

    uint e = 0;
    long long sum_num_of_edge_fc = 0;

    BuildSubFCTree(mg, degs, core, pos, klmd, e, count, sum_num_of_edge_fc);

    // Free the memory
    for (uint i = 0; i < n_vertex; i++) delete[] degs[i];
    delete[] degs;

    cout << "count = " << count << endl;

}

void FCTreeBuilderLeft::BuildSubFCTree(MultilayerGraph &mg, int **degs, uint *core, uint *pos, uint *klmd, uint e, uint &count, long long& sum_num_of_edge_fc){

    uint s = e;
    uint old_e = e;
    uint new_e = 0;

    uint n_vertex = mg.GetN();
    uint n_layer = mg.getLayerNumber();
    uint old_pos_v;

    for(int i = s; i < n_vertex; i ++){
        uint v = core[i];
        uint cnt = 0;
        bool flag = false;
        for(int j = 0; j < n_layer; j ++){
            if(degs[v][j] >= klmd[0]){
                cnt ++;
                if(cnt >= klmd[1]){
                    flag = true;
                    break;
                }
            }
        }
        if(!flag){
            old_pos_v = pos[v];
            core[old_pos_v] = core[e];
            core[e] = v;

            pos[core[old_pos_v]] = old_pos_v;
            pos[v] = e ++;

        }
    }

    if(s != e){
        new_e = peel(mg, degs, klmd, core, pos, s, e);
    }else if(s == e){
        new_e = e;
    }

    if(n_vertex - new_e > 0){
        count ++;
        PrintCoreInfor(klmd, core, new_e, n_vertex, degs, mg.getLayerNumber(), sum_num_of_edge_fc); 
        klmd[0] += 1;
        BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count,sum_num_of_edge_fc); 
        klmd[0] -= 1;
    }

    if(n_vertex - new_e > 0 && klmd[0] == 1){
        klmd[1] += 1;
        BuildSubFCTree(mg, degs, core, pos, klmd, new_e, count,sum_num_of_edge_fc);
        klmd[1] -= 1; 

    }

    restore(mg, degs, core, old_e, new_e); 

    s = old_e;
    e = old_e;   
}