import numpy as np
import pandas as pd
from datetime import datetime, date, time, timedelta
from dateutil.relativedelta import relativedelta
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.figsize'] = 15, 10

data = pd.read_csv('cubic.csv')
data = data[data.iloc[:, 2] != 'Same']

data.reset_index(drop=True, inplace=True)

# Eample: 2010-01-04
data = data[data['date'] == '2010-01-04']

# Ensure the date column is in datetime format
data['date'] = pd.to_datetime(data['date'])

# calculate maturity_date
data['maturity_date'] = pd.to_datetime(data['date']) + pd.to_timedelta(data['days'], unit='D')

data.reset_index(drop=True, inplace=True)

def schedule(maturity_date, date_today):
    coupon_dates = []
    i = pd.to_datetime(maturity_date)
    while i >= pd.to_datetime(date_today):
        coupon_dates.append(i)
        i = i - relativedelta(months=6)
    return np.flip(coupon_dates, 0)

def Get_Input_Data(data, date_today):
    prices_data = data['rate'].values
    coupons = np.zeros(len(data))  # Assuming no coupon for zero coupon bonds
    
    X_values = []
    Y_values = []
    coupon_schedule_all_bonds = []
    for j in data.index:
        maturity_date = data.loc[j, 'maturity_date']
        X_values.append((maturity_date - pd.to_datetime(date_today)).days)
        Y_values.append(data.loc[j, 'rate'])
        coupon_schedule_all_bonds.append(schedule(maturity_date, date_today))
    
    time_T = np.round(np.array(X_values) * 12 / 365, 3)
    yields = np.array(Y_values)
    
    return time_T, yields, coupon_schedule_all_bonds, prices_data, coupons, X_values, Y_values

date_today = '2010-01-04'
time_T, yields, coupon_schedule_all_bonds, prices_data, coupons, X_values, Y_values = Get_Input_Data(data, date_today)

cs = CubicSpline(time_T, yields)

# Generate rates for each day from 1 to 5400 days
maturity_days = np.arange(1, 5401)
rates = cs(maturity_days * 12 / 365)  # Convert days to months

# Output every day's maturity days and corresponding rate
maturity_days_and_rates = pd.DataFrame({'Maturity Days': maturity_days, 'Rates': rates})
print(maturity_days_and_rates)

#  save the data to a CSV file
maturity_days_and_rates.to_csv('maturity_days_and_rates.csv', index=False)

