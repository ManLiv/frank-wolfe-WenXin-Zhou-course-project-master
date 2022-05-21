#ifndef COST_HPP_
#define COST_HPP_

#include "math.h"

struct bpr {
    double capacity;
    double fft;
    double B;
    double power;

    double powerp1;
    double powerm1;
    double costanti_integral;
    double costanti_update;
    double tmp;

    bpr() :
            capacity(0.), fft(0.), B(0.), power(0.), powerp1(0.), powerm1(0.), costanti_integral(0.), costanti_update(0.), tmp(0.){
    }

    double operator()(const double& flow) const {
        return fft * (1. + B * std::pow(flow / capacity, power));
    }

    double derivative(const double& flow) const {
        return power * fft * capacity * B * std::pow(flow / capacity, powerm1);
    }

    double integral(const double& flow) const {
        return (fft * flow) + (costanti_integral * (std::pow(flow, powerp1) / powerp1));
    }

    void initialize(const double& _capacity, const double& _fft, const double& _B, const double& _power, const double& _length, const double& _toll) {
        this->capacity = _capacity;
        this->fft = _fft;
        this->B = _B;
        this->power = _power;

        this->powerp1 = this->power + 1;
        this->powerm1 = this->power - 1;
        this->costanti_integral = (this->fft * this->B) / std::pow(this->capacity, this->power);
        this->costanti_update = (fft * B) / capacity;
    }

    inline void update(const double& flow, double& weight, double& derivative) {
        double tmp = costanti_update * std::pow(flow / capacity, powerm1);

        weight = fft + tmp * flow;
        derivative = tmp * power;
    }
};

template<typename graph_t>
double compute_objective_value(const graph_t& g) {
    double f = 0.0;
    typename graph_t::edge_iterator ei, ei_end;

    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        f += g[*ei].cost_fun.integral(g[*ei].flow);
    }

    return f;
}

template<typename graph_t, typename ublas_vector>
double compute_objective_value_with_alpha(const graph_t& g, const double& alpha, const ublas_vector& direction) {
    double f = 0.0;
    typename graph_t::edge_iterator ei, ee;
    int index = 0;

    for (boost::tie(ei, ee) = boost::edges(g); ei != ee; ++ei) {
        f += g[*ei].cost_fun.integral(g[*ei].flow + alpha * direction(index));
        index++;
    }

    return f;
}

#endif /*COST_HPP_*/