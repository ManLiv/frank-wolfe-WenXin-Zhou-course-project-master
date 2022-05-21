#ifndef GRAPH_HPP_
#define GRAPH_HPP_

struct vertex_info {
    bool centroid;
};

template<typename cost_t>
struct edge_info {
    typedef cost_t cost_type;

    double weight;
    double derivative;
    cost_t cost_fun;
    double flow;
    double auxiliary_link_flow;

    edge_info() :
            weight(0.0), derivative(0.0), cost_fun(), flow(0.0), auxiliary_link_flow(0.0) {
    }

    void update(const double& flow) {
        this->flow = flow;
        this->cost_fun.update(this->flow, this->weight, this->derivative);
    }
};

#endif /*GRAPH_HPP_*/