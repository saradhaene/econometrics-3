'''
Filename: pre_assignment1.py

Purpose: First illustration about handling the data for the assignement
    
Date:
    18 February 2026  
'''
###########################################################
### Imports

import numpy as np
import pandas as pd

# global variables
s_asY = ["lwage"]
s_asX = ["exp", "wks", "bluecol", "ind", "south", "smsa", "married", "gender", "union", "edu", "colour"]
s_iY = "lwage"
s_ivX = s_asX

###########################################################
def get_data_of_year(mD, iYear):
    if 1976 <= iYear <= 1982:
        iv = (iYear - 1976) + np.arange(0, 4164, 7)
        return mD.iloc[iv,:]
    else:
        return mD

def ols_regression(y, X):
    cn, ck = X.shape
    # Calculate beta (b) = (X'X)^-1 X'y
    mxxinv = np.linalg.inv(X.T @ X)
    vb = mxxinv @ (X.T @ y)

    ve = y - (X @ vb)
    s2 = (ve.T @ ve) / (cn - ck)
    vse = np.sqrt(np.diag(mxxinv) * s2)
    t_values = vb.T / vse
    return vb, mxxinv, vse, t_values
###########################################################

def main():
    file_path = "data_assignment1.csv"
    data = pd.read_csv(file_path)
    print(data)
    data = data.iloc[:,2:] # as stated in the assignment, don't look at the first 2 rows, see if this goes right

    year_data_test = get_data_of_year(data, 1976)
    y = year_data_test[s_asY].to_numpy()
    X = np.column_stack((np.ones(len(y)), year_data_test[s_asX].to_numpy()))  # Add constant for intercept
    beta, inv_XtX, vse, t_values = ols_regression(y, X)

    results = []
    for year in range(1976, 1983):
        year_data = get_data_of_year(data, year)
        y = year_data[s_asY].to_numpy()
        X = np.column_stack((np.ones(len(y)), year_data[s_asX].to_numpy()))  # Add constant for intercept

        beta, inv_XtX, vse, t_values = ols_regression(y, X)
        results.append((beta, t_values,vse))

    pooled_y = data[s_asY].to_numpy()
    pooled_X = np.column_stack((np.ones(len(pooled_y)), data[s_asX].to_numpy()))

    pooled_beta, pooled_inv_XtX, pooled_vse, pooled_t_values = ols_regression(pooled_y, pooled_X)

    # Summary Table
    table = pd.DataFrame()
    for year, (beta, t_values, vse) in zip(range(1976, 1983), results):
    # Flatten beta and t_values before assigning to DataFrame
        table[f"{year}_coef"] = beta.flatten()
        table[f"{year}_t"] = t_values.flatten()
        table[f"{year}_vse"] = vse.flatten()

    table["pooled_coef"] = pooled_beta.flatten()
    table["pooled_t"] = pooled_t_values.flatten()
    table["pooled_vse"] = pooled_vse.flatten()
    print(table)
    

###########################################################
### call main
if __name__ == "__main__":
    main()