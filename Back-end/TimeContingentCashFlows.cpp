//TimeContingentCashFlows.cpp
#include "TimeContingentCashFlows.h"
#include "TermStructureHoLee.h"
#include "TermStructure.h"
#include <algorithm>
#include <vector>
#include <iostream>

std::vector<TimeContingentCashFlows> build_time_series_of_bond_time_contingent_cash_flows(const std::vector<double>& initial_times, 
                                                                                          const std::vector<double>& initial_cflows) {
    std::vector<TimeContingentCashFlows> vec_cf;
    std::vector<double> times = initial_times;
    std::vector<double> cflows = initial_cflows;

    while (times.size() > 0) {
        vec_cf.push_back(TimeContingentCashFlows(times, cflows));
        std::vector<double> tmp_times;
        std::vector<double> tmp_cflows;
        for (int i = 0; i < times.size(); ++i) {
            if (times[i] - 1.0 >= 0.0) {
                tmp_times.push_back(times[i] - 1.0);
                tmp_cflows.push_back(cflows[i]);
            }
        }
        times = tmp_times;
        cflows = tmp_cflows;
    }

    return vec_cf;
}

double price_european_call_option_on_bond_using_ho_lee(TermStructure* initial,
    const double& delta,
    const double& pi,
    const std::vector<double>& underlying_bond_cflow_times,
    const std::vector<double>& underlying_bond_cflows,
    const double& K,
    const double& option_time_to_maturity,
    const std::vector<double>& market_times,
    const std::vector<double>& market_prices) {

    int T = int(option_time_to_maturity + 0.0001);

    // Calibrate the Ho-Lee model
    TermStructureHoLee ho_lee_model(initial, T + 1, 0, delta, pi);
    ho_lee_model.calibrate(market_times, market_prices);
    double calibrated_delta = ho_lee_model.delta_;
    double calibrated_pi = ho_lee_model.pi_;

    // Build the term structure tree using calibrated parameters
    auto hl_tree = buildTermStructureTree(initial, T + 1, calibrated_delta, calibrated_pi);
    auto vec_cf = build_time_series_of_bond_time_contingent_cash_flows(underlying_bond_cflow_times, underlying_bond_cflows);


    for (size_t i = 0; i < hl_tree.size(); ++i) { // time i
        const auto& row = hl_tree[i];
        for (size_t j = 0; j < row.size(); ++j) { // statue j
            const auto& node = row[j];
            // std::cout << "Discount factor at row " << i << " node " << j << " is " << node.d(1.0) << std::endl;
        }
    }

    // Print vec_cf and T for debugging
    // std::cout << "T: " << T << std::endl;
    // std::cout << "vec_cf:" << std::endl;
    for (const auto& cf : vec_cf) {
        cf.print();
    }

    // std::cout << "vec_cf.sizez():" << vec_cf.size() << std::endl;
    // std::cout << "underlying_bond_cflows.sizez():" << underlying_bond_cflows.size() << std::endl;

    std::vector<double> values(T + 1, 0.0);
    for (int i = 0; i <= T; ++i) {
        // std::cout << "bonds price: " << bonds_price(vec_cf[T].times, vec_cf[T].cash_flows, hl_tree[T][i]) << std::endl;
        values[i] = std::max(0.0, bonds_price(vec_cf[T].times, vec_cf[T].cash_flows, hl_tree[T][i]) - K); // Call payoffs at maturity
        // std::cout << "check i :" << i << std::endl;
        // std::cout << "values[i] :" << values[i] << std::endl;
    }

    for (int t = T - 1; t >= 0; --t) {
        std::vector<double> values_this(t + 1, 0.0);
        for (int i = 0; i <= t; ++i) {
            values_this[i] = (pi * values[i + 1] + (1.0 - pi) * values[i]) * hl_tree[t][i].d(1);
        }
        values = values_this;
    }

    return values[0];
}

// Function to generate cash flow times and cash flows for a bond
void generate_bond_cash_flows(double face_value, double coupon_rate, double time_to_maturity,
    std::vector<double>& underlying_bond_cflow_times,
    std::vector<double>& underlying_bond_cflows) {
    int n = static_cast<int>(time_to_maturity); // Assuming annual coupon payments
    underlying_bond_cflow_times.clear();
    underlying_bond_cflows.clear();

    for (int i = 1; i <= n; ++i) {
        underlying_bond_cflow_times.push_back(static_cast<double>(i));
        underlying_bond_cflows.push_back(face_value * coupon_rate);
    }
    // Adding the face value to the last cash flow
    underlying_bond_cflows.back() += face_value;
}