#ifndef LINESEARCH_HPP_
#define LINESEARCH_HPP_

#define QUADRATIC_GAMMA 1e-8
#define LINESEARCH_THETA 0.5

#include "cost.hpp"
#include "utils.hpp"

template<typename graph_type, typename ublas_vector>
double golden_section(const graph_type& g, const ublas_vector& link_flow, const ublas_vector& auxiliary_link_flow, const double& accuracy=1e-8) {
    double LB = 0.0;
    double UB = 1.0;
    double golden_point = 0.618;
    double leftX = LB + (1. - golden_point) * (UB - LB);
    double rightX = LB + golden_point * (UB - LB);

    while (1) {
        graph_type* g_left = new graph_type(g);
        graph_type* g_right = new graph_type(g);

        ublas_vector left_flow = (1. - leftX) * link_flow + leftX * auxiliary_link_flow;
        ublas_vector right_flow = (1. - rightX) * link_flow + rightX * auxiliary_link_flow;

        typename graph_type::edge_iterator ei1, ee1;
        typename graph_type::edge_iterator ei2, ee2;
        int index1 = 0;
        int index2 = 0;

        for (boost::tie(ei1, ee1) = boost::edges((*g_left)); ei1 != ee1; ++ei1) {
            (*g_left)[*ei1].update(left_flow(index1));
            index1++;
        }

        for (boost::tie(ei2, ee2) = boost::edges((*g_right)); ei2 != ee2; ++ei2) {
            (*g_right)[*ei2].update(right_flow(index2));
            index2++;
        }

        double val_left = compute_objective_value((*g_left));
        double val_right = compute_objective_value((*g_right));

        if (val_left <= val_right) {
            UB = rightX;
        }
        else {
            LB = leftX;
        }

        if (abs(LB - UB) < accuracy) {
            double opt_theta = (rightX + leftX) / 2.0;
            return opt_theta;
        }
        else {
            if (val_left <= val_right) {
                rightX = leftX;
                leftX = LB + (1 - golden_point) * (UB - LB);
            }
            else {
                leftX = rightX;
                rightX = LB + golden_point*(UB - LB);
            }
        }

        delete g_left;
        delete g_right;
    }
}


template<typename float_t>
bool robust_equal(const float_t& A, const float_t& B) {
    float_t diff = std::abs(A - B);
    float_t absA = std::abs(A);
    float_t absB = std::abs(B);

    return (diff <= std::max(absB, absA) * std::numeric_limits<float_t>::epsilon());
}


template<typename graph_type, typename ublas_vector>
double quadratic_linesearch(const graph_type& g, const ublas_vector& direction, const double& initial_step) {
    double starting_z = compute_objective_value(g);
    double alpha = initial_step;
    double new_z = compute_objective_value_with_alpha(g, alpha, direction);
    double armijoLine = starting_z - alpha * alpha * QUADRATIC_GAMMA;

    while (!robust_equal<double>(new_z, armijoLine) && new_z > armijoLine) {
        alpha *= LINESEARCH_THETA;
        new_z = compute_objective_value_with_alpha(g, alpha, direction);
        armijoLine = starting_z - alpha * alpha * QUADRATIC_GAMMA;
    }

    return alpha;
}


#endif /*LINESEARCH_HPP_*/