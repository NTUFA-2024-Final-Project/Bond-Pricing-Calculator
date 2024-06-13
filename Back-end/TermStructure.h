#ifndef TERM_STRUCTURE_H
#define TERM_STRUCTURE_H

#include <vector>

class TermStructure {
public:
    struct NelsonSiegelParams {
        double beta0;
        double beta1;
        double beta2;
        double tau;

        NelsonSiegelParams(double b0, double b1, double b2, double t)
            : beta0(b0), beta1(b1), beta2(b2), tau(t) {}
    };

    struct Bond {
        double price;
        double faceValue;
        double couponRate;
        int frequency;
        double maturity;
    };

    virtual double r(const double& t) const = 0; // yield on zero coupon bond
    virtual double d(const double& t) const = 0; // discount factor/price of zero coupon bond
    virtual double f(const double& t1, const double& t2) const = 0; // forward rate
    virtual ~TermStructure() {};

    static double yieldFromDiscountFactor(const double& d_t, const double& t);
    static double discountFactorFromYield(const double& r, const double& t);
    static double forwardRateFromDiscountFactors(const double& d_t1, const double& d_t2, const double& time);
    static double forwardRateFromYields(const double& r_t1, const double& r_t2, const double& t1, const double& t2);
    static double nelsonSiegelSpotRate(double t, const NelsonSiegelParams& params);
    static double nelsonSiegelDiscountFactor(double t, const NelsonSiegelParams& params);
    static std::vector<double> bootstrapYieldCurve(const std::vector<Bond>& bonds);
};

class TermStructureFlat : public TermStructure {
private:
    double R_; // interest rate
public:
    TermStructureFlat(const double& r);
    virtual ~TermStructureFlat();
    virtual double r(const double& t) const override;
    virtual double d(const double& t) const override;
    virtual double f(const double& t1, const double& t2) const override;
    void setIntRate(const double& r);
};

class TermStructureInterpolated : public TermStructure {
private:
    std::vector<double> times_; // use to keep a list of yields
    std::vector<double> yields_;
    void clear();
public:
    TermStructureInterpolated();
    TermStructureInterpolated(const std::vector<double>& times, const std::vector<double>& yields);
    virtual ~TermStructureInterpolated();
    TermStructureInterpolated(const TermStructureInterpolated& other);
    TermStructureInterpolated& operator=(const TermStructureInterpolated& other);
    virtual double r(const double& T) const override;
    virtual double d(const double& t) const override;
    virtual double f(const double& t1, const double& t2) const override;
    void setInterpolatedObservations(std::vector<double>& times, std::vector<double>& yields);
    std::vector<double> getTimes() const; 
    std::vector<double> getDiscountFactors() const; 
};

#endif // TERM_STRUCTURE_H


double bonds_price(const std::vector<double>& cashflow_times,
                   const std::vector<double>& cashflows,
                   const TermStructure& d);