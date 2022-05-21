#ifndef PATH_HPP_
#define PATH_HPP_

#include <boost/graph/filtered_graph.hpp>

template<typename edge_iterator>
struct CompareEdges {
    inline bool operator()(const edge_iterator& lhs, const edge_iterator& rhs) {
        return *lhs < *rhs;
    }
};

template<typename graph_type>
struct path {
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<graph_type>::edge_iterator edge_it;
    typedef std::vector<edge_it> edge_ptr_list_type;

    double path_flow;
    edge_ptr_list_type path_edges;

    vertex_t origin;
    vertex_t destination;

    std::size_t hash;

    path() :
        path_flow(0.0), path_edges(), origin(), destination(), hash() {
    }

    path(const vertex_t& _vertex, const vertex_t& _destination) :
        path_flow(0.0), path_edges(), origin(_vertex), destination(_destination), hash() {
    }

    void sort_edges() {
        std::sort(this->path_edges.begin(), this->path_edges.end(), CompareEdges<edge_it>());
    }

    size_t n_edges() const {
        return path_edges.size();
    }

    double compute_cost(const graph_type& g) const {
        double cost(0.0);

        typename edge_ptr_list_type::const_iterator it;
        for (it = this->path_edges.begin(); it != this->path_edges.end(); ++it) {
            cost += g[**it].weight;
        }

        return cost;
    }
};

#endif /*PATH_HPP_*/