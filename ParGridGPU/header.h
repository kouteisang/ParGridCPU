#ifndef ALL_SYS_HEADER
#define ALL_SYS_HEADER

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <cuda_runtime.h>


using std::string;
using std::ofstream;
using std::unordered_map;
using std::endl;
using std::cout;
using std::ifstream;
typedef long long unsigned ll_uint;

 
#define BLK_NUMS 256
#define BLK_DIM 1024
#define BUFFER_SIZE 1000000
#define WARP_SIZE 32


typedef struct G_pointers {
    int** adj;
    int* deg;
    int* t_deg;
    int* offset;
    int* visit;
    int* precount;
    int n_vtx;
    int n_layer;
} G_pointers;//graph related

inline void chkerr(cudaError_t code){
    if (code != cudaSuccess){
        std::cout << cudaGetErrorString(code) << std::endl;
        exit(-1);
    }
}

#endif