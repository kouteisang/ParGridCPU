#include "Graph/MultilayerGraph.h"
#include "utils.h"
#include "Algorithm/ParGridGPU.cuh"
#include "Algorithm/ParGridGPU_block.cuh"
#include "Algorithm/ParGridGPU_k.cuh"




enum Algorithm{
    llist = 1,
    klist = 2,
    llist_block = 3,
};

int main(int argc, char* argv[]){

    cudaSetDevice(1);

    string dataset = "example";
    int order = 0;
    int alg = 1;

    for(int i = 1; i < argc; i ++){
        string arg = argv[i];
        if (arg == "-d" && i + 1 < argc) {
            dataset = argv[++i];
        }else if(arg == "-a" && i+1 < argc){
            alg = std::stoi(argv[++i]);
        }
    }


    // load the dataset
    MultilayerGraph mg;
    // relative path is relative to the execute file
    mg.LoadFromFile("./dataset/"+dataset+"/");
    mg.SetGraphOrder(0);
    uint *order_list = mg.GetOrder();
    for(int i = 0; i < mg.getLayerNumber(); i ++){
        cout << order_list[i] << " ";;
    }
    cout << endl;
    mg.PrintStatistics();

    int n_vertex = mg.GetN();
    int n_layer = mg.getLayerNumber();

    int *degs;
    degs = new int[n_vertex * n_layer];
    for(int v = 0; v < n_vertex; v ++){
        for(int l = 0; l < n_layer; l ++){
            degs[v * n_layer + l] = mg.GetGraph(l).GetAdjLst()[v][0];
        }
    }
   

    G_pointers data_pointers;
    data_pointers.n_vtx = n_vertex;
    data_pointers.n_layer = n_layer;

    // // memory alloc adj list for each layer
    std::vector<int*> d_adj_list_ptrs(n_layer);

    int* h_offset;
    h_offset = new int[n_layer * (n_vertex+1)];
    int cnt = 0;

    int *num_edge;
    num_edge = new int[n_layer];

    for(int l = 0; l < n_layer; l ++){

        std::vector<int> h_adj_list;
        int totalsum = 0;
        h_offset[cnt ++] = totalsum;

        uint** adj_lst = mg.GetGraph(l).GetAdjLst();
        for(int v = 0; v < n_vertex; v ++){
            int offset = mg.GetGraph(l).GetAdjLst()[v][0];
            totalsum += offset;
            h_offset[cnt ++] = totalsum;
            for(int nb = 1; nb <= offset; nb ++){
                h_adj_list.push_back(adj_lst[v][nb]);
            }
        }
        
        // for(int uu = 0; uu < h_adj_list.size(); uu ++){
        //     cout << h_adj_list[uu] << " ";
        // }

        int* d_list;
        int len = h_adj_list.size();
        num_edge[l] = len;
        chkerr(cudaMalloc(&d_list, len * sizeof(int)));
        cudaMemcpy(d_list, h_adj_list.data(), sizeof(int) * len, cudaMemcpyHostToDevice);
        d_adj_list_ptrs[l] = d_list;
    }



    malloc_graph_gpu_memory(data_pointers, degs, h_offset, d_adj_list_ptrs);

    cudaEvent_t start, stop; // Calculate time
    cudaEventCreate(&start); // Calculate time
    cudaEventCreate(&stop);  // Calculate time
    
    cudaEventRecord(start, 0);

    switch (alg)
    {
        case Algorithm::llist:
            gpu_baseline_de(data_pointers, degs);
            break;
        case Algorithm::klist:
            gpu_baseline_de_klist(data_pointers, degs);
            break;
        case Algorithm::llist_block:
            gpu_baseline_block(data_pointers, degs);
            break;
        default:
            break;
    }

    



    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    float gpu_time = 0;
    cudaEventElapsedTime(&gpu_time, start, stop);
    std::cout << "GPU time = " << gpu_time*1.0/1000 << " s" << std::endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);


    // ll_uint *id2vtx = new ll_uint[mg.GetN()];
    // mg.LoadId2VtxMap(id2vtx);



}