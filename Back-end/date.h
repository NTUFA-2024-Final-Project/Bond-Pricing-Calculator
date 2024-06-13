//date.h
#ifndef DATE_H
#define DATE_H
#include <algorithm>  // Include this for std::min
#include <cmath>      // For floor and leap year calculations
#include <ctime>

enum class DayCountConvention {
    Thirty360,
    Thirty365,
    Actual360,
    Actual365,
    ActualActual
};

class date {
protected:
    int year_;
    int month_;
    int day_;
public:
    date(); //default constructor
    date(const int& d, const int& m, const int& y);
    bool valid() const;
    int day() const;
    int month() const;
    int year() const;
    void set_day(const int& day);
    void set_month(const int& month);
    void set_year(const int& year);
    double years_until(const date& other, DayCountConvention dcc) const;
    bool is_leap_year(int) const;
    static date next_date(const date& d);
    static date previous_date(const date& d);
    static date current_date();
    date operator ++(); //prefix
    date operator ++(int); //postfix
    date operator --(); //prefix
    date operator --(int); //postfix
};
bool operator == (const date&, const date&); //comparisonoperators
bool operator != (const date&, const date&);
bool operator < (const date&, const date&);
bool operator > (const date&, const date&);
bool operator <= (const date&, const date&);
bool operator >= (const date&, const date&);

#endif // DATE_H
