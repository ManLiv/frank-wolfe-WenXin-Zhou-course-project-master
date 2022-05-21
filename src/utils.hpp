#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <quadmath.h>
#include <boost/numeric/ublas/matrix.hpp>
#include "dijkstra_misc.hpp"
#include "path.hpp"
#include <boost/numeric/ublas/vector.hpp>


static std::ostream& operator<<(std::ostream &o, const __float128 &value) {
    char buf[128];

    std::stringstream ss;
    ss << "%+-#*." << FLT128_DIG << "Qe";

    assert((size_t) quadmath_snprintf(buf, sizeof buf, ss.str().c_str(), FLT128_DIG, value) < sizeof buf);

    o << buf;

    return o;
}

static std::istream& operator>>(std::istream &i, __float128 &value) {
    std::string str;
    i >> str;
    value = strtoflt128(str.c_str(), NULL);
    return i;
}

template<typename graph_type, typename matrix_type, typename edge_matrix_type>
void compute_min_tree(const graph_type &g, const typename graph_type::vertex_descriptor& r,
        std::vector<typename graph_type::vertex_descriptor> &p_star, const matrix_type& D, const uint& destinations,
        const bool& all_centroid,
        const edge_matrix_type& edge_matrix) {

    typedef typename boost::edge_bundle_type<graph_type>::type edge_info_type;

    try {
        if (all_centroid) {
            boost::dijkstra_shortest_paths(g, r, boost::predecessor_map(&p_star[0]).weight_map(boost::get(&edge_info_type::weight, g))); //.
        } else {
            boost::dijkstra_shortest_paths(boost::make_filtered_graph(g, O_edge_filter<graph_type>(g, r)), r,
//            boost::dijkstra_shortest_paths(boost::make_filtered_graph(g, O_edge_filter<graph_type>(g, r)), r,
                    boost::predecessor_map(&p_star[0]).weight_map(boost::get(&edge_info_type::weight, g))
#ifdef USE_DIJKSTRA_VISITOR
                            .visitor(dijkstra_end_tree_visitor<typename graph_type::vertex_descriptor, matrix_type>(r, D, destinations))
#endif
                            );
        }
    } catch (found_goal& dummy) {
    }
}


template<typename path_type, typename p_star_type, typename edge_matrix_type>
void build_path(path_type& path, const p_star_type& p_star, const edge_matrix_type& edge_matrix) {
    // Costruzione del cammino
    typename path_type::vertex_t target(path.destination);
    boost::hash_combine(path.hash, target);
    while (target != path.origin) {
        path.path_edges.push_back(edge_matrix(p_star[target], target));
        target = p_star[target];
        boost::hash_combine(path.hash, target);
    }
}


template<typename graph_type, typename ublas_vector>
double get_directional_derivative(const graph_type& g, const ublas_vector& direction) {
    typename boost::graph_traits<graph_type>::edge_iterator ei, ee;
    int index = 0;
    double directional_derivative = 0.0;

    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        directional_derivative += g[*ei].weight * direction(index);
        index++;
    }

    return directional_derivative;
}


template<typename graph_type, typename ublas_vector>
double get_dHd(const graph_type& g, const ublas_vector& direction) {
    typename boost::graph_traits<graph_type>::edge_iterator ei, ee;
    int index = 0;
    double dHd = 0.0;

    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        dHd += g[*ei].derivative * direction(index) * direction(index);
        index++;
    }

    return dHd;
}


template<typename graph_type, typename paths_matrix_type, typename p_star_type, typename mat_type, typename edge_matrix_type>
void init_graph(graph_type& g, paths_matrix_type& paths_matrix, const bool& all_centroid, p_star_type& p_star, const mat_type& D, const std::vector<uint>& destination_count, const edge_matrix_type& edge_matrix) {
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_desc_type;
    typedef typename boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef typename paths_matrix_type::value_type paths_list_type;
    typedef typename paths_list_type::value_type path_type;

    typename boost::graph_traits<graph_type>::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].update(0.0);
    }

    // int num_paths = 0;

    typename mat_type::const_iterator1 it1;
    std::vector<uint>::const_iterator itd;

    for (it1 = D.begin1(), itd = destination_count.begin(); it1 != D.end1(); ++it1, ++itd) {
        vertex_desc_type origin = it1.index1();

        if (*itd == 0) {
            continue;
        }

        compute_min_tree(g, origin, p_star, D, *itd, all_centroid, edge_matrix);

        for (typename mat_type::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
            vertex_desc_type destination = it2.index2();
            double demand = *it2;

            if (demand == 0) {
                continue;
            }

            paths_matrix(origin, destination).push_back(path_type(origin, destination));
            path_type& path = paths_matrix(origin, destination).back();
            build_path(path, p_star, edge_matrix);

            path.sort_edges();
            path.path_flow = demand;

            for (uint i = 0; i < path.n_edges(); i++) {
                edge_desc_type current_edge = *(path.path_edges[i]);
                g[current_edge].update(g[current_edge].flow + demand);
            }

            // ++num_paths;
        }
    }

    // return num_paths;
}


template<typename graph_type, typename paths_matrix_type, typename mat_type, typename edge_matrix_type, typename ublas_vector>
void all_or_nothing_assignment(graph_type& g, paths_matrix_type& paths_matrix, const bool& all_centroid, const mat_type& D, const std::vector<uint>& destination_count, const edge_matrix_type& edge_matrix, ublas_vector& auxiliary_link_flow) {
    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_desc_type;
    typedef typename boost::graph_traits<graph_type>::edge_descriptor edge_desc_type;
    typedef typename paths_matrix_type::value_type paths_list_type;
    typedef typename paths_list_type::value_type path_type;

    typename boost::graph_traits<graph_type>::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        g[*ei].auxiliary_link_flow = 0.0;
    }

    std::vector<vertex_desc_type> _p_star;
    typename mat_type::const_iterator1 it1;
    std::vector<uint>::const_iterator itd;

    for (it1 = D.begin1(), itd = destination_count.begin(); it1 != D.end1(); ++it1, ++itd) {
        vertex_desc_type origin = it1.index1();

        if (*itd == 0) {
            continue;
        }

        _p_star = std::vector<vertex_desc_type>(boost::num_vertices(g));
        compute_min_tree(g, origin, _p_star, D, *itd, all_centroid, edge_matrix);

        for (typename mat_type::const_iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2) {
            vertex_desc_type destination = it2.index2();
            double demand = *it2;

            if (demand == 0) {
                continue;
            }

            paths_matrix(origin, destination).push_back(path_type(origin, destination));
            path_type& path = paths_matrix(origin, destination).back();
            build_path(path, _p_star, edge_matrix);

            path.sort_edges();
            path.path_flow = demand;

            for (uint i = 0; i < path.n_edges(); i++) {
                edge_desc_type current_edge = *(path.path_edges[i]);
                g[current_edge].auxiliary_link_flow += demand;
            }
        }
    }

    typename boost::graph_traits<graph_type>::edge_iterator ei1, ee1;
    int index = 0;
    for (boost::tie(ei1, ee1) = boost::edges(g); ei1 != ee1; ++ei1) {
        auxiliary_link_flow(index) = g[*ei1].auxiliary_link_flow;
        index++;
    }
}


template<typename graph_type, typename mat_type, typename edge_matrix_type, typename paths_matrix_type, typename centroids_type>
double measurement(const graph_type& g, const mat_type& D,
        const std::vector<uint>& destination_count,
        const bool& all_centroid, const edge_matrix_type& edge_matrix, 
        paths_matrix_type& paths_matrix, const centroids_type& centroids){

    typedef typename boost::graph_traits<graph_type>::vertex_descriptor vertex_type;
    typedef typename paths_matrix_type::value_type list_type;
    typedef std::list<path<graph_type> > path_list_type;
    typedef typename path_list_type::value_type path_type;

    typename mat_type::const_iterator1 it1;
    typename mat_type::const_iterator2 it2;
    vertex_type origin, destination;
    unsigned int r;
    std::vector<vertex_type> _p_star;
    double sum_d_times_miu = 0.0;

    for (r = 0; r < centroids.size(); ++r) {
        if (destination_count[r] == 0.){
            continue;
        }

        _p_star = std::vector<vertex_type>(boost::num_vertices(g));

        compute_min_tree(g, centroids[r], _p_star, D, destination_count[r], all_centroid, edge_matrix);

        it1 = D.begin1();
        std::advance(it1, r);

        for (it2 = it1.begin() ; it2 != it1.end() ; ++it2) {
            if (*it2 == 0.){
                continue;
            }

            origin = it2.index1();
            destination = it2.index2();

            path_type p(origin, destination);
            build_path(p, _p_star, edge_matrix);
            p.sort_edges();

            double demand = *it2;
            p.path_flow = demand;

            double minimal_path_cost = p.compute_cost(g);
            sum_d_times_miu += minimal_path_cost * demand;
        }
    }

    return sum_d_times_miu;
}


#endif /*UTILS_HPP_*/