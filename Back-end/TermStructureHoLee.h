//TermStructureHoLee.h
#ifndef TERM_STRUCTURE_HO_LEE_H
#define TERM_STRUCTURE_HO_LEE_H

#include "TermStructure.h"
#include <vector>
#include "Eigen/Dense"

class TermStructureHoLee : public TermStructure {
public:
    TermStructure* initial_term_;
    int n_; //number of step
    int i_; //statue i
    double delta_;
    double pi_; // Implied Binomial Probabilit

    TermStructureHoLee(TermStructure* fitted_term, const int& n, const int& i, const double& delta, const double& pi);

    virtual double r(const double& t) const override { return 0.0; } // Not implemented
    virtual double d(const double& T) const override;
    virtual double f(const double& t1, const double& t2) const override { return 0.0; } // Not implemented

    void print(const size_t& row, const size_t& node) const; // Declare the print method

    void calibrate(const std::vector<double>& market_times, const std::vector<double>& market_prices); // Calibration method

};

std::vector<std::vector<TermStructureHoLee>> buildTermStructureTree(TermStructure* initial, 
                                                                    const int& no_steps, 
                                                                    const double& delta, 
                                                                    const double& pi);

#endif // TERM_STRUCTURE_HO_LEE_H