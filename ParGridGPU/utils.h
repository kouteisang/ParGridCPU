#ifndef UTILS_H
#define UTILS_H

#include "header.h"


void malloc_graph_gpu_memory(G_pointers &p, int* h_deg, int* h_offset, std::vector<int*> d_adj_list_ptrs){

    chkerr(cudaMalloc((&p.adj), sizeof(int*) * p.n_layer));
    cudaMemcpy(p.adj, d_adj_list_ptrs.data(), sizeof(int*) * p.n_layer, cudaMemcpyHostToDevice);

    // deg
    chkerr(cudaMalloc(&(p.deg), p.n_layer * p.n_vtx * sizeof(int)));
    chkerr(cudaMemcpy(p.deg, h_deg, p.n_layer * p.n_vtx * sizeof(int), cudaMemcpyHostToDevice));

    // tdeg
    chkerr(cudaMalloc(&(p.t_deg), p.n_layer * p.n_vtx * sizeof(int)));

    // visit
    chkerr(cudaMalloc((&p.visit), p.n_vtx * sizeof(int)));
    cudaMemset(p.visit, 0, p.n_vtx * sizeof(int));

    chkerr(cudaMalloc((&p.precount), p.n_vtx * sizeof(int)));
    cudaMemset(p.precount, 0, p.n_vtx * sizeof(int));

    // offset
    chkerr(cudaMalloc((&p.offset), (p.n_layer)*(p.n_vtx+1)* sizeof(int)));
    chkerr(cudaMemcpy(p.offset, h_offset, (p.n_layer)*(p.n_vtx+1) * sizeof(int), cudaMemcpyHostToDevice));
    
}

#endif