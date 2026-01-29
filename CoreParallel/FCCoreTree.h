#ifndef FCTreeCore_H
#define FCTreeCore_H

#include "Graph/MultilayerGraph.h"

// parallel core node
struct coreNodeP{

    unsigned int k{0};
    unsigned int lmd{0};
    unsigned int length{0};
    int** degs{nullptr};
    uint *core{nullptr}; // result core
    bool *valid{nullptr}; // valid vertex
    // std::vector<bool> valid; 
    bool *t_valid{nullptr}; // valid vertex
    vector<uint> valid_vertex;
    vector<uint> diff_res;


    std::vector<bool> validsync;
    std::vector<int> res_vtx;
    int total{0}; // the number of total peeled vertex;


    uint *valid2{nullptr}; // this is just for test
    int invalidIndex{0}; // this is just for test

    uint* invalid{nullptr};
    uint* cnts{nullptr};

    bool* res{nullptr};

    uint e{0}; // e is the start point

    coreNodeP* left{nullptr};
    coreNodeP* right{nullptr};

    coreNodeP(){};

    coreNodeP(unsigned int _k, unsigned int _lmd, unsigned int _length){
        k = _k;
        lmd = _lmd;
        length = _length;
    };
     
};



class FCCoreTree
{
private:
    coreNodeP *node;
public:
    FCCoreTree(/* args */);
    ~FCCoreTree();


    FCCoreTree(unsigned int k, unsigned int lmd, unsigned int length){
       node = new coreNodeP(k, lmd, length);
    }

    // get the root Node
    inline coreNodeP* getNode() const {
        return node;
    }

    void traversal(coreNodeP *node, uint num_vertex){
        if(node != nullptr && node->k != 0 && node->lmd != 0){
            
            uint cnt = 0;
            for(uint v = 0; v < num_vertex; v ++){
                if(node->res[v] == true){
                    cnt ++;
                }
            }

            cout << node->k << " " << node->lmd << " cnt = " <<  cnt << endl;

            traversal(node->left, num_vertex);
            traversal(node->right, num_vertex);
        }
    }
};






#endif

