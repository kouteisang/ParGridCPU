


#include "MultilayerGraph.h"


MultilayerGraph::~MultilayerGraph(){
    delete[] graph_layers;    
}

void MultilayerGraph::SetGraphOrder(int o){
    if (o == 0){
        
    }else if(o == 1){
      
    }else if(o == 2){

    }
}

void MultilayerGraph::GetGraphFile(const std::string &input_path, std::vector<std::string> &graph_files){

    std::string conf_file = input_path + "mlg.conf";
    std::string line;

    auto fin = std::ifstream(conf_file);
    while (fin.peek() != EOF) {
        getline(fin, line);

        if (!line.empty()) {
            graph_files.emplace_back(line);
        }
    }
    fin.close();

}

void MultilayerGraph::LoadFromFile(const std::string &input_path){
    std::vector<std::string> graph_files;
    std::vector<edge *> edge_buf;
    std::vector<uint> edge_size;

    map_file = input_path + "__vtx_map.txt";
    std::unordered_map<ll_uint, uint> vtx2id;
    std::basic_ofstream<char> map_file_out;

    GetGraphFile(input_path, graph_files);

    n_layers = graph_files.size(); 
    //How many layer we have for this graph
    graph_layers = new Graph[n_layers];

    order = new uint[n_layers];
    for(int i = 0; i < n_layers; i ++){
        order[i] = i;
    }


    edge_buf.resize(n_layers);
    edge_size.resize(n_layers);

    map_file_out = ofstream(map_file);
    for(uint i = 0; i < n_layers; i ++){
        edge_size[i] = LoadLayer(input_path + graph_files[i], edge_buf[i], vtx2id, map_file_out);
    } 
    map_file_out.close();

    n = vtx2id.size();
    for (uint i = 0; i < n_layers; i++) {
        graph_layers[i].BuildFromEdgeLst(edge_buf[i], n, edge_size[i]);
        delete[] edge_buf[i];
    }


}

uint MultilayerGraph::LoadLayer(const std::string &graph_file, edge *&edge_buf, unordered_map<ll_uint, uint> &vtx2id, std::basic_ofstream<char> &map_file_out){
    
    uint uid, vid;
    uint edge_buf_size, num_of_vtx, num_of_edge;
    uint u, v;

    edge *tmp_edge_buf;

    std::basic_ifstream<char> graph_in;
    graph_in = std::ifstream(graph_file);

    num_of_vtx = (uint) vtx2id.size();

    // load graph
    num_of_edge = 0;
    edge_buf_size = DEFAULT_EDGE_BUF_SIZE;
    edge_buf = new edge[edge_buf_size];
    
    while (graph_in.peek() != EOF) {

        graph_in >> u >> v;
        // remove self-loop
        if (u != v) {
            auto iter1 = vtx2id.find(u);
            if (iter1 != vtx2id.end()) {
                uid = iter1->second;
            } else {
                uid = num_of_vtx++;
                vtx2id.emplace(u, uid);
                map_file_out << u << " " << uid << endl;
            }

            auto iter2 = vtx2id.find(v);
            if (iter2 != vtx2id.end()) {
                vid = iter2->second;
            } else {
                vid = num_of_vtx++;
                vtx2id.emplace(v, vid);
                map_file_out << v << " " << vid << endl;
            }

            if (num_of_edge + 2 > edge_buf_size) {
                edge_buf_size = edge_buf_size << 1;
                tmp_edge_buf = new edge[edge_buf_size];

                memcpy(tmp_edge_buf, edge_buf, num_of_edge * sizeof(edge));
                delete[] edge_buf;
                edge_buf = tmp_edge_buf;
            }

            edge_buf[num_of_edge++] = edge(uid, vid);
            edge_buf[num_of_edge++] = edge(vid, uid);
        }
    }

    // for(int i = 0; i < num_of_edge; i ++){
    //     cout << "s  = " << edge_buf[i].s << " t = " << edge_buf[i].t << endl;
    // }

    graph_in.close();
    return num_of_edge;

} 

void MultilayerGraph::PrintStatistics() {
    uint max_m = 0, sum_m = 0, m;

    for (uint i = 0; i < n_layers; i++) {
        m = graph_layers[i].getM() >> 1;
        max_m = std::max(max_m, m);
        sum_m += m;
    }

    std::cout << "|L| = " << n_layers << ", |V| = " << n << ", |E| = " << sum_m << ", max_|E_i| = " << max_m << std::endl;
}

void MultilayerGraph::LoadId2VtxMap(ll_uint *id2vtx) {
    ll_uint u;
    uint uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        id2vtx[uid] = u;
    }
    fin.close();
}

void MultilayerGraph::LoadVtx2IdMap(unordered_map<ll_uint, uint> &vtx2id) {
    ll_uint u;
    uint uid;

    ifstream fin(map_file);
    while (fin.good() && !fin.eof()) {
        fin >> u >> uid;
        vtx2id[u] = uid;
    }
    fin.close();
}

MultilayerGraph* MultilayerGraph::Load(const string &file) {
    auto mg = new MultilayerGraph();
    mg->LoadFromFile(file);
    return mg;
}
