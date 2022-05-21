#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <iostream>
#include <iomanip>
#include <malloc.h>
#include <boost/numeric/ublas/io.hpp>

#include "src/io.hpp"
#include "src/dijkstra_misc.hpp"
#include "src/graph.hpp"
#include "src/cost.hpp"
#include "src/utils.hpp"
#include "src/path.hpp"
#include "src/frank_wolfe.hpp"

int main(int argc, char** argv) {
    mallopt(M_MMAP_MAX, 0);
    mallopt(M_TRIM_THRESHOLD, -1);

    typedef bpr cost_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, vertex_info, edge_info<cost_type> > graph_type;
    typedef boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef boost::graph_traits<graph_type>::edge_iterator edge_iterator;

    typedef boost::numeric::ublas::matrix<double> matrix_type;
    typedef boost::numeric::ublas::compressed_matrix<edge_iterator> edge_matrix_type;

    typedef path<graph_type> path_type;
    typedef std::list<path_type> path_list_type;
    typedef boost::numeric::ublas::matrix<path_list_type> paths_matrix_type;
    
    typedef boost::numeric::ublas::vector<double> ublas_vector;

    int num_centroids;
    bool all_centroids;
    double total_demand;
    int num_of_edges;
    int num_OD_pairs;
    char network_filename[] = "data/ChicagoSketch_net.txt";
    char trips_filename[] = "data/ChicagoSketch_trips.txt";

    graph_type g;
    matrix_type D;
    std::vector<uint> destination_count;
    std::vector<graph_type::vertex_descriptor> centroids;
    std::vector<graph_type::vertex_descriptor> p_star;
    paths_matrix_type paths_matrix;
    edge_matrix_type edge_matrix;

    load_network(network_filename, g, centroids, num_centroids, all_centroids);
    num_OD_pairs = load_trips(trips_filename, D, destination_count, total_demand);
    p_star = std::vector<graph_type::vertex_descriptor>(boost::num_vertices(g));
    paths_matrix = paths_matrix_type(D.size1(), D.size2());

    num_of_edges = boost::num_edges(g);
    edge_matrix.resize(boost::num_vertices(g), boost::num_vertices(g), false);
    std::pair<edge_iterator, edge_iterator> edges = boost::edges(g);
    bool is_multi_graph = false;
    for (edge_iterator edge = edges.first; edge != edges.second; ++edge) {
        vertex_type src = boost::source(*edge, g);
        vertex_type dst = boost::target(*edge, g);
        if (edge_matrix.find_element(src, dst)) {
            is_multi_graph = true;
        }
        edge_matrix(boost::source(*edge, g), boost::target(*edge, g)) = edge;
    }
    if (is_multi_graph) {
        exit(0);
    }

    init_graph(g, paths_matrix, all_centroids, p_star, D, destination_count, edge_matrix);
    ublas_vector final_link_flow(num_of_edges, 0);
    convex_combination_method(g, paths_matrix, all_centroids, centroids, D, destination_count, edge_matrix, final_link_flow, num_of_edges);

    edge_iterator ei1, ee1;
    double obj = compute_objective_value(g);
    std::cout << "Objective value = " << obj << std::endl;

    return 0;
}