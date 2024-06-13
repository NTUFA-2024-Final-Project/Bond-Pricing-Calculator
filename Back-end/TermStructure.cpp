#include "TermStructure.h"
#include <cmath>
#include <algorithm>

double TermStructure::yieldFromDiscountFactor(const double& d_t, const double& t) {
    return (-std::log(d_t) / t);
}

double TermStructure::discountFactorFromYield(const double& r, const double& t) {
    return std::exp(-r * t);
}

double TermStructure::forwardRateFromDiscountFactors(const double& d_t1, const double& d_t2, const double& time) {
    return (std::log(d_t1 / d_t2)) / time;
}

double TermStructure::forwardRateFromYields(const double& r_t1, const double& r_t2, const double& t1, const double& t2) {
    return (r_t2 * t2 - r_t1 * t1) / (t2 - t1);
}

double TermStructure::nelsonSiegelSpotRate(double t, const NelsonSiegelParams& params) {
    double term1 = params.beta0;
    double term2 = params.beta1 * (1 - exp(-t / params.tau)) / (t / params.tau);
    double term3 = params.beta2 * ((1 - exp(-t / params.tau)) / (t / params.tau) - exp(-t / params.tau));
    return term1 + term2 + term3;
}

double TermStructure::nelsonSiegelDiscountFactor(double t, const NelsonSiegelParams& params) {
    double spotRate = nelsonSiegelSpotRate(t, params);
    return exp(-spotRate * t);
}

std::vector<double> TermStructure::bootstrapYieldCurve(const std::vector<Bond>& bonds) {
    std::vector<double> spotRates(bonds.size(), 0.0);
    for (size_t i = 0; i < bonds.size(); ++i) {
        const auto& bond = bonds[i];
        double cashFlowSum = 0.0;
        double coupon = bond.couponRate * bond.faceValue / bond.frequency;
        for (int j = 1; j <= bond.maturity * bond.frequency; ++j) {
            double time = j / static_cast<double>(bond.frequency);
            double df = j - 1 < i ? discountFactorFromYield(spotRates[j - 1], time) : 1.0;  // Use previously bootstrapped rates
            cashFlowSum += coupon * df;
        }
        // Add the face value payment
        double lastDF = discountFactorFromYield(spotRates[i], bond.maturity);
        double totalValue = cashFlowSum + bond.faceValue * lastDF;
        double targetPrice = bond.price;

        // Simple solver: Adjust the spot rate until the bond is correctly priced
        double error = totalValue - targetPrice;
        while (std::abs(error) > 0.0001) {
            spotRates[i] += (targetPrice - totalValue) / 100000.0;  // Adjust rate
            lastDF = discountFactorFromYield(spotRates[i], bond.maturity);
            totalValue = cashFlowSum + bond.faceValue * lastDF;
            error = totalValue - targetPrice;
        }
    }
    return spotRates;
}

// Implementations for TermStructureFlat

TermStructureFlat::TermStructureFlat(const double& r) : R_(r) {}

TermStructureFlat::~TermStructureFlat() {}

double TermStructureFlat::r(const double& t) const {
    if (t >= 0) return R_;
    return 0;
}

double TermStructureFlat::d(const double& t) const {
    return discountFactorFromYield(R_, t);
}

double TermStructureFlat::f(const double& t1, const double& t2) const {
    double d1 = d(t1);
    double d2 = d(t2);
    return forwardRateFromDiscountFactors(d1, d2, t2 - t1);
}

void TermStructureFlat::setIntRate(const double& r) {
    R_ = r;
}

// Implementations for TermStructureInterpolated

TermStructureInterpolated::TermStructureInterpolated() {
    clear();
}

TermStructureInterpolated::TermStructureInterpolated(const std::vector<double>& times, const std::vector<double>& yields) {
    clear();
    setInterpolatedObservations(const_cast<std::vector<double>&>(times), const_cast<std::vector<double>&>(yields));
}

TermStructureInterpolated::~TermStructureInterpolated() {
    clear();
}

TermStructureInterpolated::TermStructureInterpolated(const TermStructureInterpolated& other) {
    times_ = other.times_;
    yields_ = other.yields_;
}

TermStructureInterpolated& TermStructureInterpolated::operator=(const TermStructureInterpolated& other) {
    if (this != &other) {
        times_ = other.times_;
        yields_ = other.yields_;
    }
    return *this;
}

double TermStructureInterpolated::r(const double& T) const {
    // Implementing linear interpolation for simplicity
    if (times_.empty()) return 0.0;

    auto it = std::lower_bound(times_.begin(), times_.end(), T);
    if (it == times_.begin()) return yields_.front();
    if (it == times_.end()) return yields_.back();

    size_t index = it - times_.begin();
    double t1 = times_[index - 1];
    double t2 = times_[index];
    double y1 = yields_[index - 1];
    double y2 = yields_[index];

    double interpolated_yield = y1 + (T - t1) * (y2 - y1) / (t2 - t1);
    return interpolated_yield;
}

double TermStructureInterpolated::d(const double& t) const {
    return discountFactorFromYield(r(t), t);
}

double TermStructureInterpolated::f(const double& t1, const double& t2) const {
    double d1 = d(t1);
    double d2 = d(t2);
    return forwardRateFromDiscountFactors(d1, d2, t2 - t1);
}

void TermStructureInterpolated::setInterpolatedObservations(std::vector<double>& times, std::vector<double>& yields) {
    clear();
    if (times.size() != yields.size()) return;

    times_ = times;
    yields_ = yields;
}

void TermStructureInterpolated::clear() {
    times_.clear();
    yields_.clear();
}



std::vector<double> TermStructureInterpolated::getTimes() const {
    return times_;
}

std::vector<double> TermStructureInterpolated::getDiscountFactors() const {
    std::vector<double> discount_factors;
    for (const double& time : times_) {
        discount_factors.push_back(d(time));
    }
    return discount_factors;
}

double bonds_price(const std::vector<double>& cashflow_times,
    const std::vector<double>& cashflows,
    const TermStructure& d) {
    double p = 0;
    for (int i = 0; i < cashflow_times.size(); i++) {
        p += d.d(cashflow_times[i]) * cashflows[i];
    }
    return p;
}