#include "ParGridGPU.cuh"

__global__ void trykernel(int** adj ,int* deg, int* offset, int n_layer, int n_vtx){

    // printf("Hello world\n");
    // for(int v = 0; v < n_vtx; v ++){
    //     for(int l = 0; l < n_layer; l ++){
    //         printf("%d ,", deg[v*n_layer+l]);
    //     }
    //     printf("\n\n");
    // }
    // for(int l = 0; l < n_layer; l ++){
    //     for(int v = 0; v < n_vtx; v ++){
    //         int begin = l * (n_vtx+1) + v;
    //         int end =  l * (n_vtx+1) + v + 1;
    //         printf("%d, %d\n", offset[begin], offset[end]);
    //     }
    // }

    // printf("\n");
    // printf("\n");
    // for(int l = 0; l  < n_layer; l ++){
    //     int len = offset[(n_vtx+1)*(l+1) - 1];
    //     // printf("len = %d", len);
    //     int* adj_l = adj[l];
    //     for(int v = 0; v < len; v ++){
    //         printf("%d,", adj_l[v]);
    //     }
    //     printf("\n");
    // }

    // // for(int l = 0; l  < n_layer; l ++){
    //     int len = 18;
    //     int* adj_l = adj[0];
    //     for(int v = 0; v < 18; v ++){
    //         printf("%d,", adj_l[v]);
    //     }
    //     printf("\n\n\n");
    // // }

}


__global__ void scan(int* global_buffer, int* buf_count, int* deg, int* precount, int* visit, int n_vtx, int n_layer, int k, int l){
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ int sh_buf_count;
    __shared__ int* t_global_buffer;

    if(threadIdx.x == 0){
        sh_buf_count = 0;
        t_global_buffer = global_buffer + blockIdx.x * BUFFER_SIZE;
    }
    __syncthreads();

    for(int v = tid; v < n_vtx; v += BLK_DIM * BLK_NUMS){
        if(visit[v] == 1) continue;
        int count = 0;
        int start = v*n_layer;
        int end = (v+1)*n_layer;
        for(int d = start; d < end; d ++){
            count += (deg[d] >= k);
        }
        if(count < l){
            visit[v] = 1; // visit 设置为false
            count = 0;  // count 设置为 0
            int pos = atomicAdd(&sh_buf_count, 1);
            t_global_buffer[pos] = v;
            // printf("v = %d\n", v);
        }
        precount[v] = count;
    }
    __syncthreads();

    if(threadIdx.x == 0){
        buf_count[blockIdx.x] = sh_buf_count;
    
    }
}

__global__ void update(int* global_buffer, int* buf_count, int *deg, int **adj, int* offset, int* precount, int* visit, int n_vtx, int n_layer, int k, int lambda, int* global_count){
    
    __shared__ int start, end;
    __shared__ int* t_global_buffer;

    int warp_per_block = blockDim.x / WARP_SIZE;
    int warp_id = threadIdx.x / WARP_SIZE;
    int lane_id = threadIdx.x % WARP_SIZE;
    int start_prime, end_prime;
    if(threadIdx.x == 0){
        t_global_buffer = global_buffer + blockIdx.x * BUFFER_SIZE;
        start = 0;
        end = buf_count[blockIdx.x]; // The end position of the buffer
        // printf("id = %d, end = %d\n", blockIdx.x, end);
    } 

    __syncthreads();

    while(true){
        __syncthreads();
        // printf("end = %d\n", end);
        if(start >= end) break; // All the thread break the iteration
        start_prime = start + warp_id; // Get the vertex id position
        end_prime = end; // Get the last position of the vertex id
        __syncthreads();
        if(start_prime >= end_prime) continue; // The vertex position is larger than the number of valid vertices in the buffer
        if(threadIdx.x == 0){
            start = min(start + warp_per_block, end); // update the start position
        }
        int v = t_global_buffer[start_prime]; // Get the vertex id

        for(int l = 0; l < n_layer; l ++){
            __syncwarp();
            int offset_start = offset[l * (n_vtx+1) + v]; // offset of v 
            int offset_end = offset[l * (n_vtx+1) + v + 1]; // offset of v
            int* adj_l = adj[l];
            while (true){
                __syncwarp();
                if(offset_start >= offset_end) break;
                int uid = offset_start + lane_id;
                offset_start = offset_start + WARP_SIZE; // update the offset position, each thread maintain its own offset_start
                if(uid >= offset_end) continue; // This vertex does not has so many neighbouthood
                int u = adj_l[uid]; // v's out-neighbouthood u
                if(visit[u] == 1) continue;
                // printf("l = %d, v = %d, u = %d\n", l, v, u);
                int originDeg = atomicSub(&deg[u*n_layer + l], 1);
                if (originDeg == k){
                    // printf("u = %d, originDeg = %d, precount = %d\n", u, originDeg, precount[u]);
                    int originCnt = atomicSub(&precount[u], 1);
                    // printf("originCnt = %d\n", originCnt); 
                    if(originCnt == lambda && visit[u] == 0){
                        // printf("u = %d\n", u);
                        visit[u] = 1;
                        precount[u] = 0;
                        int end_pos = atomicAdd(&end, 1);
                        t_global_buffer[end_pos] = u;
                    }
                }
            }

        }
    }

    if(threadIdx.x == 0 && end > 0){
        atomicAdd(global_count, end);
    }


}

void gpu_baseline_de(G_pointers &p, int* dges){
    // printf("Here?\n");
    // trykernel<<<1, 1>>>(p.adj, p.deg, p.offset, p.n_layer, p.n_vtx);
    // cudaDeviceSynchronize(); 

    int* global_count = 0;
    chkerr(cudaMalloc(&global_count, sizeof(int)));

    int* buf_count;
    chkerr(cudaMalloc(&buf_count, sizeof(int) * BLK_NUMS));
    cudaMemset(buf_count, 0, sizeof(int) * BLK_NUMS);

    int* global_buffer;
    chkerr(cudaMalloc(&global_buffer, sizeof(int) * BLK_NUMS * BUFFER_SIZE));
    

    int n_layer = p.n_layer;
    int n_vtx = p.n_vtx;

    int k = 0;
    int count = 0;

    for(int l = 1; l <= n_layer; l ++ ){
        k = 1;
        cudaMemset(global_count, 0, sizeof(int));
        cudaMemset(p.precount, 0, sizeof(int)*n_vtx);
        cudaMemset(p.visit, 0, p.n_vtx * sizeof(int)); // flag = false means has not visited
        chkerr(cudaMemcpy(p.t_deg, p.deg, p.n_vtx * p.n_layer * sizeof(int), cudaMemcpyDeviceToDevice));
        count = 0;
        while(count < n_vtx){ 
            cudaMemset(buf_count, 0, sizeof(int) * BLK_NUMS);
            scan<<<BLK_NUMS, BLK_DIM>>>(global_buffer, buf_count, p.t_deg, p.precount, p.visit,n_vtx, n_layer, k, l);
            update<<<BLK_NUMS, BLK_DIM>>>(global_buffer, buf_count, p.t_deg, p.adj, p.offset, p.precount, p.visit, n_vtx, n_layer, k, l, global_count);
            chkerr(cudaMemcpy(&count, global_count, sizeof(int), cudaMemcpyDeviceToHost));
            // printf("count = %d\n", count);
            // printf("l = %d, k = %d, valid = %d\n", l, k, n_vtx - count);
            if(count < n_vtx){
                k ++;
            }else if(count >= n_vtx){
               break;
            }
        }
        if(k == 1 && count == n_vtx) break;
    }
}

// 0 1 2
// 0 3 5