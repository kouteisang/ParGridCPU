#include "Graph.h"

Graph::~Graph(){

}

void Graph::BuildFromGraph(uint n_, uint m_, uint **adj_lst_, uint *adj_lst_buf_){
    n = n_;
    m = m_;
    adj_lst = adj_lst_;
    adj_lst_buf = adj_lst_buf_;
}



// use CSR to store the graph
void Graph::BuildFromEdgeLst(edge *edge_buf, uint num_of_vtx, uint num_of_edge){

    uint i, j, pj;

    n = num_of_vtx;
    // adj_lst is a pointer pointer, therefore it is 2D. The first dimension is the number of vertex
    // The second dimension is the number of outgoing edges I think.
    adj_lst = new uint *[n]; 
    adj_lst_buf = new uint [num_of_edge + n];

    std::sort(edge_buf, edge_buf + num_of_edge); // 对边的大小进行排序


    i = 0; // index of edge buf;
    j = 0; // index of edge lst buf;
   
    for(uint v = 0; v < n; v ++){
        adj_lst[v] = &adj_lst_buf[j];
        if(edge_buf[i].s > v || i >= num_of_edge) adj_lst_buf[j++] = 0;
        else{
            pj = j++;
            adj_lst_buf[j++] = edge_buf[i++].t;
            while(i < num_of_edge && edge_buf[i].s == v){
                if (edge_buf[i].t == edge_buf[i - 1].t) i++;  // duplicated edge
                else adj_lst_buf[j++] = edge_buf[i++].t;
            }
            adj_lst_buf[pj] = j - pj - 1;
            m += adj_lst_buf[pj];
            if (adj_lst_buf[pj] > max_deg) max_deg = adj_lst_buf[pj];   
        }
    } 


}