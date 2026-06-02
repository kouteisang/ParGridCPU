#ifndef GRAPH_H
#define GRAPH_H

#include "../header.h"

struct edge{
    uint s{0};
    uint t{0};

    edge() = default;

    edge(uint _s, uint _t){
        s = _s;
        t = _t;
    };

    bool operator<(const edge &e) const {
        return s < e.s || (s == e.s && t < e.t);
    }

};

class Graph{
protected:
    uint n{0};
    uint m{0};
    uint max_deg{0};

    uint *adj_lst_buf{nullptr};
    uint **adj_lst{nullptr};

public:
    Graph() = default;
    ~Graph();

    void BuildFromEdgeLst(edge *edge_buf, uint num_of_vtx, uint num_of_edge);
    void BuildFromGraph(uint n_, uint m_, uint **adj_lst_, uint *adj_lst_buf_);


    [[nodiscard]] uint getN() const{
        return n;
    }

    [[nodiscard]] uint getM() const{
        return m;
    }

    void SetNull(){
        adj_lst = nullptr;
        adj_lst_buf = nullptr;
    }

    [[nodiscard]] uint **GetAdjLst() {
        return adj_lst;
    }

};



#endif