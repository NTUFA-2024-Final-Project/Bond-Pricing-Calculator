//TimeContingentCashFlows.h
#ifndef TIME_CONTINGENT_CASH_FLOWS_H
#define TIME_CONTINGENT_CASH_FLOWS_H

#include <vector>
#include "TermStructure.h"
#include <iostream>

class TimeContingentCashFlows {
public:
    std::vector<double> times;
    std::vector<double> cash_flows;

    TimeContingentCashFlows(const std::vector<double>& in_times, const std::vector<double>& in_cflows)
        : times(in_times), cash_flows(in_cflows) {}

    int no_cflows() const { return int(times.size()); }

    void print() const {
        // std::cout << "Times: ";
        for (const auto& t : times) {
            // std::cout << t << " ";
        }
        // std::cout << "\nCash Flows: ";
        for (const auto& cf : cash_flows) {
            // std::cout << cf << " ";
        }
        // std::cout << std::endl;
    }
};

std::vector<TimeContingentCashFlows> build_time_series_of_bond_time_contingent_cash_flows(const std::vector<double>& initial_times,
                                                                                          const std::vector<double>& initial_cflows);

double price_european_call_option_on_bond_using_ho_lee(TermStructure* initial, 
                                                        const double& delta, 
                                                        const double& pi, 
                                                        const std::vector<double>& underlying_bond_cflow_times, 
                                                        const std::vector<double>& underlying_bond_cflows,
                                                        const double& K, 
                                                        const double& option_time_to_maturity,
                                                        const std::vector<double>& market_times,
                                                        const std::vector<double>& market_prices);
void generate_bond_cash_flows(double face_value, double coupon_rate, double time_to_maturity,
    std::vector<double>& underlying_bond_cflow_times,
    std::vector<double>& underlying_bond_cflows);
#endif // TIME_CONTINGENT_CASH_FLOWS_H