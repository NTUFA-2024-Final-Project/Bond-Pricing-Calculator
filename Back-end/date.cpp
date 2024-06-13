//date.cpp
#include "date.h"   
date::date() {
	day_ = 1;
	month_ = 1;
	year_ = 1900;
}
date::date(const int& day, const int& month, const int& year) {
	day_ = day;
	month_ = month;
	year_ = year;
};
int date::day() const { return day_; };
int date::month() const { return month_; };
int date::year() const { return year_; };
void date::set_day(const int& day) { date::day_ = day; };
void date::set_month(const int& month) { date::month_ = month; };
void date::set_year(const int& year) { date::year_ = year; };
bool date::valid() const {
	// Basic checks
	if (year_ < 0) return false;
	if (month_ < 1 || month_ > 12) return false;
	if (day_ < 1 || day_ > 31) return false;

	// Check for months with 30 days
	if ((month_ == 4 || month_ == 6 || month_ == 9 || month_ == 11) && day_ > 30) {
		return false;
	}

	// Check for February (month 2)
	if (month_ == 2) {
		// Check if it's a leap year
		bool isLeapYear = (year_ % 4 == 0 && (year_ % 100 != 0 || year_ % 400 == 0));
		if (day_ > 29 || (day_ == 29 && !isLeapYear)) {
			return false;
		}
	}

	return true; // Passed all checks
}


date date::next_date(const date& d) {
	if (!d.valid()) return date(); // Consider handling this case differently.
	date ndat(d.day() + 1, d.month(), d.year()); // First try adding a day
	if (ndat.valid()) return ndat;
	ndat = date(1, d.month() + 1, d.year()); // Then try adding a month
	if (ndat.valid()) return ndat;
	return date(1, 1, d.year() + 1); // Must be next year
}

date date::previous_date(const date& d) {
	if (!d.valid()) return date(); // Consider handling this case differently.
	date pdat(d.day() - 1, d.month(), d.year()); // Try same month
	if (pdat.valid()) return pdat;
	pdat = date(31, d.month() - 1, d.year()); // Try previous month
	for (int day = 31; day >= 28; day--) {
		pdat = date(day, d.month() - 1, d.year());
		if (pdat.valid()) return pdat;
	}
	return date(31, 12, d.year() - 1); // Try previous year
}



double date::years_until(const date& other, DayCountConvention dcc) const {
	// Check the validity of both the current and other date instances. Return 0.0 if either is invalid.
	if (!this->valid() || !other.valid()) return 0.0;

	// Extract the day, month, and year components of the current and other date instances.
	int d1 = day_, d2 = other.day_;
	int m1 = month_, m2 = other.month_;
	int y1 = year_, y2 = other.year_;

	// Switch statement to handle different day count conventions
	switch (dcc) {
	case DayCountConvention::Thirty360: {
		// For the 30/360 convention, days are adjusted so that all months have 30 days.
		d1 = std::min(d1, 30);  // Adjust the current day not to exceed 30.
		d2 = (d1 == 30 && d2 == 31) ? 30 : d2; // If the current day was adjusted to 30 and the other day is 31, adjust the other day to 30 as well.
		// Calculate the year fraction assuming each year has 360 days (12 months * 30 days each).
		return ((y2 - y1) * 360 + (m2 - m1) * 30 + (d2 - d1)) / 360.0;
	}
	case DayCountConvention::Thirty365: {
		// Similar to the 30/360, but each year is assumed to have 365 days.
		d1 = std::min(d1, 30);
		d2 = (d1 == 30 && d2 == 31) ? 30 : d2;
		return ((y2 - y1) * 365 + (m2 - m1) * 30 + (d2 - d1)) / 365.0;
	}
	case DayCountConvention::Actual360: {
		// The actual number of days are counted, but each year is assumed to have 360 days.
		int days = (y2 - y1) * 360 + (m2 - m1) * 30 + (d2 - d1);
		return days / 360.0;
	}
	case DayCountConvention::Actual365: {
		// The actual number of days are counted, and each year is assumed to have 365 days.
		int days = (y2 - y1) * 365 + (m2 - m1) * 30 + (d2 - d1);
		return days / 365.0;
	}
	case DayCountConvention::ActualActual: {
		// This method counts each day between the two dates and divides by the actual number of days in the end year.
		// This considers leap years individually by adjusting the denominator according to the end year.
		int days = 0;
		date tmp = *this;  // Start with the current date.
		while (tmp < other) {  // Increment the date one day at a time until it reaches or surpasses 'other'.
			tmp++;  // Increment the date by one day.
			days++; // Count each day.
		}
		return days / (is_leap_year(y2) ? 366.0 : 365.0); // Divide the total days by either 366 or 365, depending on whether the end year is a leap year.
	}
	default:
		return 0.0;  // In case of an undefined day count convention, return 0.
	}
}

bool date::is_leap_year(int year) const {
	return (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
}

date date::current_date() {
	time_t t = time(nullptr);
	tm* now = localtime(&t);
	return date(now->tm_mday, now->tm_mon + 1, now->tm_year + 1900);
}

//operators overloading
date date::operator++(int) {//postfixoperator
	date d = *this;
	*this = next_date(d);
	return d;
}
date date::operator ++() {//prefixoperator
	*this = next_date(*this);
	return*this;
}
date date::operator --(int) {//postfixoperator, returncurrent value
	date d = *this;
	*this = previous_date(*this);
	return d;
}
date date::operator --() {//prefixoperator, returnnewvalue
	*this = previous_date(*this);
	return*this;
};
bool operator == (const date& d1, const date& d2) {//check forequality
	if (!(d1.valid() && (d2.valid()))) { return false; }; /*ifdatesnotvalid,not clearwhat todo.
	alternative: throwexception*/
	return((d1.day() == d2.day()) && (d1.month() == d2.month()) && (d1.year() == d2.year()));
};
bool operator < (const date& d1, const date& d2) {
	if (!(d1.valid() && (d2.valid()))) { return false; }; //seeaboveremark
	if (d1.year() == d2.year()) {//sameyear
		if (d1.month() == d2.month()) {//samemonth
			return(d1.day() < d2.day());
		}
		else {
			return(d1.month() < d2.month());
		};
	}
	else {//differentyear
		return(d1.year() < d2.year());
	};
};
bool operator<=(const date& d1, const date& d2) {
	if (d1 == d2) { return true; }
	return(d1 < d2);
}
bool operator>=(const date& d1, const date& d2) {
	if (d1 == d2) { return true; };
	return(d1 > d2);
};
bool operator>(const date& d1, const date& d2) { return !(d1 <= d2); };
bool operator !=(const date& d1, const date& d2) { return !(d1 == d2); }
