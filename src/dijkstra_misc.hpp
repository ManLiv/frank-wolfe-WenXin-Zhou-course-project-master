#ifndef DIJKSTRA_MISC_HPP_
#define DIJKSTRA_MISC_HPP_

#include <boost/graph/dijkstra_shortest_paths.hpp>

struct found_goal {
};


template<class Vertex>
class dijkstra_goal_visitor: public boost::default_dijkstra_visitor {
public:
    dijkstra_goal_visitor(const Vertex& goal) :
            m_goal(goal) {
    }

    template<class Graph>
    void examine_vertex(const Vertex& u, const Graph& g) {
        if (u == m_goal) {
            throw found_goal();
        }
    }
private:
    Vertex m_goal;
};


template<class Vertex, class MatrixType>
class dijkstra_end_tree_visitor: public boost::default_dijkstra_visitor {
public:
    dijkstra_end_tree_visitor(const Vertex& r, const MatrixType& _D, const uint& destinations) :
            D(_D), root(r), n_to_be_visited(destinations) {
    }

    template<class Graph>
    void examine_vertex(const Vertex& u, const Graph& g) {
        if (g[u].centroid && D(root, u) > 0.) {
            if (--n_to_be_visited == 0) {
                throw found_goal();
            }
        }
    }
private:
    const MatrixType& D;
    Vertex root;
    int n_to_be_visited;
};


template<typename graph_type>
struct OD_node_filter {

    typedef typename graph_type::vertex_descriptor v_type;

    const graph_type *g_m;
    v_type v_source, v_destination;

    OD_node_filter() :
            g_m(NULL) {
    }

    OD_node_filter(const graph_type &g, const v_type& v, const v_type& v_end) :
            g_m(&g), v_source(v), v_destination(v_end) {
    }

    template<typename vertex_type>
    bool operator()(const vertex_type& v) const {
        if ((*g_m)[v].centroid) {
            if (v != v_destination && v != v_source) {
                return false;
            }
        }

        return true;
    }
};


template<typename graph_type>
struct O_edge_filter {

    typedef typename graph_type::vertex_descriptor v_type;

    const graph_type *g_m;
    v_type v_source;

    O_edge_filter() :
            g_m(NULL) {
    }

    O_edge_filter(const graph_type &g, const v_type& v) :
            g_m(&g), v_source(v) {
    }

    template<typename edge_type>
    bool operator()(const edge_type& e) const {
        v_type o = boost::source(e, *g_m);

        if ((*g_m)[o].centroid && o != v_source) {
            return false;
        }

        return true;
    }
};


template<typename graph_type, class MatrixType>
struct O_edge_filter_plus {

    typedef typename graph_type::vertex_descriptor v_type;

    const graph_type *g_m;
    v_type v_source;
    const MatrixType* D;

    O_edge_filter_plus() :
        g_m(NULL), D(NULL) {
    }

    O_edge_filter_plus(const graph_type &g, const v_type& v, const MatrixType* _D) :
            g_m(&g), v_source(v), D(_D) {
    }

    template<typename edge_type>
    bool operator()(const edge_type& e) const {
        v_type o = boost::source(e, *g_m);
        v_type d = boost::target(e, *g_m);

        if ((*g_m)[o].centroid && o != v_source) {
            return false;
        }

        if ((*g_m)[d].centroid && (*D)(v_source, d) == 0.){
            return false;
        }

        return true;
    }
};

#endif /*DIJKSTRA_MISC_HPP_*/