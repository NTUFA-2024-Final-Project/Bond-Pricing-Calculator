#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "date.h"
#include "date.cpp"
#include "TermStructure.h"
#include "TermStructure.cpp"
#include "TermStructureHoLee.h"
#include "TermStructureHoLee.cpp"
#include "TimeContingentCashFlows.h"
#include "TimeContingentCashFlows.cpp"

// source: https://github.com/yhirose/cpp-httplib
#include "httplib.h"

// source: https://github.com/nlohmann/json
// reference: https://hackmd.io/@tico88612/cpp-json-file-tutorial?print-pdf#/
#include "json.hpp"

using namespace std;
using namespace httplib;
using json = nlohmann::json;

// Define a structure for storing CSV data
struct CSVRow {
    vector<string> columns;
};

// Function to parse a CSV row
CSVRow parseCSVRow(const string &line) {
    CSVRow row;
    stringstream ss(line);
    string column;
    while (getline(ss, column, ',')) {
        row.columns.push_back(column);
    }
    return row;
}

// Function to read CSV and extract data
void readCSV(const string& filename, vector<double>& times, vector<double>& yields) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Unable to open file" << endl;
        return;
    }

    string line;
    // Skip the header
    getline(file, line);
    while (getline(file, line)) {
        CSVRow row = parseCSVRow(line);
        if (row.columns.size() >= 2) {
            times.push_back(stod(row.columns[0]));
            yields.push_back(stod(row.columns[1]));
        }
    }
    file.close();
}

// CORS headers
void setup_cors_headers(Response &res) {
    res.set_header("Access-Control-Allow-Origin", "*");
    res.set_header("Access-Control-Allow-Methods", "GET, POST, PUT, DELETE, OPTIONS");
    res.set_header("Access-Control-Allow-Headers", "Accept, Content-Type, Content-Length, Accept-Encoding, X-CSRF-Token, Authorization");
}

int main() {

    Server svr;

    svr.Options("/.*", [](const Request &req, Response &res) {
        setup_cors_headers(res);
        res.status = 200; // No content
    });
    
    svr.Post("/calculate", [](const Request &req, Response &res) {

        setup_cors_headers(res);

        auto params = json::parse(req.body);

        // get json data fron front-end
        double K = params["k"];
        string maturity_date = params["maturity_date"];

        double delta = params["delta"];
        double pi = params["pi"];
        double coupon_rate = params["coupon_rate"];
        double face_value = params["face_value"];
        int dcc_case = stoi(params["day_count_convention"]["value"].get<string>());

        // split maturity date
        int eYear, eMonth, eDay;
        sscanf(maturity_date.c_str(), "%d-%d-%d", &eYear, &eMonth, &eDay);

        vector<double> times_;
        vector<double> yields_;
        
        // Read data from .csv file
        readCSV("maturity_days_and_rates.csv", times_, yields_);

        // cubic term structure
        TermStructureInterpolated *initial = new TermStructureInterpolated(times_, yields_);
        
        date startingDate = date::current_date();
        date expirationDate(eDay, eMonth, eYear);

        // Determine the day count convention based on combox selection
        DayCountConvention dcc;

        switch (dcc_case) {//should be imputed by dropdown list
            case 0: dcc = DayCountConvention::Thirty360; break;
            case 1: dcc = DayCountConvention::Thirty365; break;
            case 2: dcc = DayCountConvention::Actual360; break;
            case 3: dcc = DayCountConvention::Actual365; break;
            case 4: dcc = DayCountConvention::ActualActual; break;
            default:
                //error message "Invalid day count convention selected");
                return 1;
        }
        
        // Calculate the time to maturity using years_until method
        double timeToMaturity = startingDate.years_until(expirationDate, dcc);
        
        //callable_bond_information        
        double time_to_maturity = 10;
        vector<double> underlying_bond_cflow_times;
        vector<double> underlying_bond_cflows;
        generate_bond_cash_flows(face_value, coupon_rate, time_to_maturity, underlying_bond_cflow_times, underlying_bond_cflows);

        vector<double> getTime = initial->getTimes();
        vector<double> getDiscountFactor = initial->getDiscountFactors();

        double callable_bond_price = price_european_call_option_on_bond_using_ho_lee(initial,
                                                                                    delta, 
                                                                                    pi, 
                                                                                    underlying_bond_cflow_times,
                                                                                    underlying_bond_cflows, 
                                                                                    K, 
                                                                                    timeToMaturity,
                                                                                    getTime,
                                                                                    getDiscountFactor);

        delete initial;
        
        // test
        cout << "callable bond price: " << callable_bond_price << endl;

        json response;
        response["callable_bond_price"] = callable_bond_price;
        
        res.set_content(response.dump(), "application/json");
        
        return 0;
    });

    cout << "Server is running at http://localhost:3001" << endl;
    svr.listen("localhost", 3001);

    return 0;
}


