
#ifndef MULTILAYERGRAPH_H
#define MULTILAYERGRAPH_H


#include "Graph.h"
#include <cstring> 

const uint DEFAULT_EDGE_BUF_SIZE = 5000;


class MultilayerGraph{
private:
    uint n{0};
    uint n_layers{0};
    uint *order{nullptr};
    Graph *graph_layers{nullptr};
    string map_file;



    static uint LoadLayer(const std::string &graph_file, edge *&edge_buf, std::unordered_map<ll_uint, uint> &vtx2id, std::basic_ofstream<char> &map_file_out);

public:


    MultilayerGraph() = default;
    ~MultilayerGraph();

    void LoadFromFile(const std::string &input_path);
    void GetGraphFile(const std::string &input_path,  std::vector<std::string> &graph_files);
    void SetGraphOrder(int o);

    void LoadId2VtxMap(ll_uint * id2vtx);
    void LoadVtx2IdMap(unordered_map<ll_uint, uint> &vtx2id);

    static MultilayerGraph* Load(const string &file);
    void PrintStatistics();

    inline Graph &GetGraph(uint i) {
        return graph_layers[order[i]];
    }

    inline uint getLayerNumber(){
        return n_layers;
    }

    [[nodiscard]] inline uint GetN() const {
        return n;
    }

    [[nodiscard]] inline uint* GetOrder() const {
        return order;
    }

    
};

#endif