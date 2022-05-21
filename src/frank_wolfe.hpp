#ifndef FRANK_WOLFE_HPP_
#define FRANK_WOLFE_HPP_

#include "utils.hpp"
#include "linesearch.hpp"
#include <fstream>
#include <float.h>
#include <chrono>

template<typename graph_type, typename edge_matrix_type, typename ublas_vector, typename centroids_type, typename paths_matrix_type, typename mat_type>
void convex_combination_method(graph_type& g, paths_matrix_type& paths_matrix, const bool& all_centroid, const centroids_type& centroids, const mat_type& D, const std::vector<uint>& destination_count, const edge_matrix_type& edge_matrix, ublas_vector& final_link_flow, const int& num_of_edges) {
    bool solved = false;
    double accuracy = 1e-4;
    std::ofstream outFile; // storing link flow on each link
    outFile.open("result_flow.csv", std::ios::out);
    outFile << "link,link_flow" << std::endl;

    std::ofstream outFile1; // storing iteration-error, time-error
    outFile1.open("result_error.csv", std::ios::out);
    outFile1 << "iteration,time,error" << std::endl;

    int index = 0;
    ublas_vector link_flow(num_of_edges, 0);
    typename boost::graph_traits<graph_type>::edge_iterator ei, ee;
    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        link_flow(index) = g[*ei].flow;
        index++;
    }

    int it = 1;
    double alpha;
    ublas_vector auxiliary_link_flow(num_of_edges, 0);
    double err;
    std::cout << "it        err" << std::endl;
    auto begin = std::chrono::system_clock::now();

    while (!solved) {
        double sum_d_times_miu = 0.0;
        double sum_t_times_v = 0.0;
        all_or_nothing_assignment(g, paths_matrix, all_centroid, D, destination_count, edge_matrix, auxiliary_link_flow);
        ublas_vector direction = auxiliary_link_flow - link_flow;
        double initial_step = std::abs(get_directional_derivative(g, direction)) / get_dHd(g, direction);
        alpha = quadratic_linesearch(g, direction, initial_step);
        link_flow = link_flow + alpha * direction;

        // first update all info of edges
        typename boost::graph_traits<graph_type>::edge_iterator ei1, ee1;
        int index1 = 0;
        for (boost::tie(ei1, ee1) = boost::edges(g); ei1 != ee1; ++ei1) {
            g[*ei1].update(link_flow(index1));
            index1++;
        }

        // then calculate convergence conditions
        sum_d_times_miu = measurement(g, D, destination_count, all_centroid, edge_matrix, paths_matrix, centroids);

        typename boost::graph_traits<graph_type>::edge_iterator ei2, ee2;
        for (boost::tie(ei2, ee2) = boost::edges(g); ei2 != ee2; ++ei2) {
            sum_t_times_v += g[*ei2].flow * g[*ei2].weight;
        }

        err = std::abs(sum_d_times_miu - sum_t_times_v) / sum_t_times_v;
        std::cout << it << "        " << err << std::endl;
        auto this_time = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(this_time - begin);
        auto beginning_to_now = double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
        outFile1 << it << "," << beginning_to_now << "," << err << std::endl;

        if (err < accuracy) {
            solved = true;
            typename boost::graph_traits<graph_type>::edge_iterator ei3, ee3;
            for (boost::tie(ei3, ee3) = boost::edges(g); ei3 != ee3; ++ei3) {
                outFile << *ei3 << "," << g[*ei3].flow << std::endl;
            }
            final_link_flow = link_flow;
        }
        else {
            it += 1;
        }
    }

    outFile.close();
    outFile1.close();
}

#endif /*FRANK_WOLFE_HPP_*/