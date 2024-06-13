//TermStructureHoLee.cpp
#include "TermStructureHoLee.h"
#include "TermStructure.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"

TermStructureHoLee::TermStructureHoLee(TermStructure* fitted_term, const int& n, const int& i, const double& delta, const double& pi)
    : initial_term_(fitted_term), n_(n), i_(i), delta_(delta), pi_(pi) {}

inline double hT(const double& T, const double& delta, const double& pi) {
    return (1.0 / (pi + (1 - pi) * pow(delta, T)));
}

double TermStructureHoLee::d(const double& T) const {
    double d = initial_term_->d(T + n_) / initial_term_->d(n_);


    for (int j = 1; j < n_; ++j) {
        double hT1 = hT(T + (n_ - j), delta_, pi_);
        double hT2 = hT(n_ - j, delta_, pi_);
        d *= hT1 / hT2;

    }

    double hT_final = hT(n_, delta_, pi_);
    double pow_delta = pow(delta_, T * (n_ - i_));
    d *= hT_final * pow_delta;

    return d;
}

std::vector<std::vector<TermStructureHoLee>> buildTermStructureTree(TermStructure* initial, 
                                                                    const int& no_steps, 
                                                                    const double& delta, 
                                                                    const double& pi) {
    std::vector<std::vector<TermStructureHoLee>> hl_tree;

    for (int t = 0; t < no_steps; ++t) {
        hl_tree.push_back(std::vector<TermStructureHoLee>());
        for (int j = 0; j <= t; ++j) {
            TermStructureHoLee hl(initial, t, j, delta, pi); 
            hl_tree[t].push_back(hl);
        }
    }
    return hl_tree;
}

void TermStructureHoLee::print(const size_t& row, const size_t& node) const {
    std::cout << "1 year Discount factor at row " << row << " node " << node << " is " << d(1.0) << std::endl;
}

// Functor for Levenberg-Marquardt algorithm
struct Functor {
    const TermStructureHoLee& model;
    const std::vector<double>& market_times;
    const std::vector<double>& market_prices;

    Functor(const TermStructureHoLee& model, const std::vector<double>& market_times, const std::vector<double>& market_prices)
        : model(model), market_times(market_times), market_prices(market_prices) {}

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
        double delta = x[0];
        double pi = x[1];

        for (size_t i = 0; i < market_times.size(); ++i) {
            TermStructureHoLee temp_model(model);
            temp_model.delta_ = delta;
            temp_model.pi_ = pi;
            fvec[i] = temp_model.d(market_times[i]) - market_prices[i];
        }

        return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const {
        double delta = x[0];
        double pi = x[1];
        double epsilon = 1e-8;

        for (size_t i = 0; i < market_times.size(); ++i) {
            // Partial derivative with respect to delta
            TermStructureHoLee temp_model_delta(model);
            temp_model_delta.delta_ = delta + epsilon;
            temp_model_delta.pi_ = pi;
            double f_delta = temp_model_delta.d(market_times[i]);

            temp_model_delta.delta_ = delta - epsilon;
            double f_delta_minus = temp_model_delta.d(market_times[i]);

            fjac(i, 0) = (f_delta - f_delta_minus) / (2 * epsilon);

            // Partial derivative with respect to pi
            TermStructureHoLee temp_model_pi(model);
            temp_model_pi.delta_ = delta;
            temp_model_pi.pi_ = pi + epsilon;
            double f_pi = temp_model_pi.d(market_times[i]);

            temp_model_pi.pi_ = pi - epsilon;
            double f_pi_minus = temp_model_pi.d(market_times[i]);

            fjac(i, 1) = (f_pi - f_pi_minus) / (2 * epsilon);
        }

        return 0;
    }

    int inputs() const { return 2; }
    int values() const { return market_times.size(); }
};

void TermStructureHoLee::calibrate(const std::vector<double>& market_times, const std::vector<double>& market_prices) {
    Eigen::VectorXd x(2);
    x[0] = delta_;
    x[1] = pi_;

    Functor functor(*this, market_times, market_prices);
    Eigen::LevenbergMarquardt<Functor> lm(functor);
    lm.minimize(x);

    delta_ = x[0];
    pi_ = x[1];
}
