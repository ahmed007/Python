"""
This File contains functions:
    1. EulersDiscretization
    2. regression_method2_call
    3. regression_method2_Put
"""
import numpy as np
import math
from scipy.stats import norm 
import pandas as pd
from scipy.optimize import curve_fit
import scipy as sy

"""
Function name: EulersDiscretization
Description: This function calculates stock price matrix using Euler's scheme
Inp  parameters: Strike price,Stock price, risk free rate(r),dividend(delta),timestep size=delta_t
N->Number of paths to simulate, T-> Expiry date
Return Value: Stock Price matrix
"""

def EulersDiscretizationRegression(S0,K,r,delta,sigma,delta_t,N,T,seed):
    stprice_matrix = np.zeros((N,int(T/delta_t)+1),dtype=float)
    stprice_matrix[:,0] = S0    
    np.random.seed(seed)
    Brownian_matrix=np.random.standard_normal(size=(N,int(T/delta_t)+1))*np.sqrt(delta_t)
    for j in range(0,N):
        for i in range(1,int(T/delta_t)+1):
            DW=Brownian_matrix[j][i]
            #DW=np.sum(Brownian_matrix[j][(1*(i-1)+1):(1*i + 1)])
            var = stprice_matrix[j][i-1] + (r-delta)*stprice_matrix[j][i-1]*delta_t+sigma*stprice_matrix[j][i-1]*DW
            #prevdw=DW
            stprice_matrix[j][i] = var
            #St.append(var)
            
    return stprice_matrix

"""
Callback function passed to curve_fit for regression
It must take the independent variable as the first argument and the parameters to fit as separate remaining arguments.
"""
def func(x,a0,a1,a2,a3):
    return a0 + a1*x + a2*(x**2) + a3*(x**3)

"""
Function name: regression_method_call
Description: This function calculates the timee zero price of American Call using regression analysis
Regression is carried out on stock price using polyfit 
[polyfit=Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg] of degree deg to points (x, y). \
Returns a vector of coefficients p that minimises the squared error.]
Return Value: time zero price of American Call
"""

def regression_method2_call(stprice_matrix,K,T,delta_t,N,r):
    M = int(T/delta_t)
    gk = np.zeros((N,1),dtype=float)
    tau_k=np.zeros((N,1),dtype=float)
    gk[:,0] = np.maximum(stprice_matrix[:,M] - K,0) #At maturity i am calculating call payoff 
    tau_k[:,0]=M 
    C_hat_vector = np.zeros((N,1),dtype=float)
    x = np.zeros((N,1),dtype=float)
    y= np.zeros((N,1),dtype=float)
    
    for j in range(M - 1,0,-1):
        #y=gk[:,0]*np.exp(-r *(tau_k[:,0] - j) * delta_t)  #i take the PV of payoff as dependent variable
        for i in range (0,N):
            if stprice_matrix[i,j] > K: #in the money
                x[i] = stprice_matrix[i,j]
                y[i]=gk[i,0]*np.exp(-r *(tau_k[i,0] - j) * delta_t)
            else:
                continue
#        a ,cov= curve_fit(func, x[:,0],y[:,0]) #calculating regression weights
#        C_hat_vector = func(x[:,0],a[0],a[1],a[2],a[3])
        reg = np.polyfit(x[:,0],y[:,0], 7)
        C_hat_vector = np.polyval(reg, x[:,0])                
        for i in range (0,N):
            if (np.maximum(stprice_matrix[i,j]-K,0)) >= C_hat_vector[i]:
                gk[i,0] = np.maximum(stprice_matrix[i,j]-K,0)
                tau_k[i,0]=j
            else:
                continue

    C_hat_0 = 0
    for i in range (0,N):
        C_hat_0 = C_hat_0 + np.exp(-r*tau_k[i,0]*delta_t)*gk[i,0]
    C_hat_0 = C_hat_0/N
    timezeroTrueVal = stprice_matrix[0,0]-K
    V_0 = np.maximum(timezeroTrueVal,C_hat_0)
    
    return ("%.6f" % V_0)
    
"""
Function name: regression_method_put
Description: This function calculates the timee zero price of American Put using regression analysis
Regression is carried out on stock price using polyfit 
[polyfit=Fit a polynomial p(x) = p[0] * x**deg + ... + p[deg] of degree deg to points (x, y). \
Returns a vector of coefficients p that minimises the squared error.]
Return Value: time zero price of American Put
"""
def regression_method2_put(stprice_matrix,K,T,delta_t,N,r):
    M = int(T/delta_t) 
    gk = np.zeros((N,1),dtype=float)
    tau_k=np.zeros((N,1),dtype=float)
    
    gk[:,0] = np.maximum(K-stprice_matrix[:,M],0) #At maturity i am calculating call payoff 
    tau_k[:,0]=M 
    P_hat_vector = np.zeros((N,1),dtype=float)
    x = np.zeros((N,1),dtype=float)
    y= np.zeros((N,1),dtype=float)
    
    for j in range(M - 1,0,-1):
        #y=gk[:,0]*np.exp(-r *(tau_k[:,0] - j) * delta_t)  #i take the PV of payoff as dependent variable
        for i in range (0,N):
            if stprice_matrix[i,j] < K: #in the money put options
                x[i] = stprice_matrix[i,j]
                y[i]=gk[i,0]*np.exp(-r *(tau_k[i,0] - j) * delta_t)
            else:
                continue
#        a ,cov= curve_fit(func, x[:,0],y[:,0]) #calculating regression weights
#        P_hat_vector = func(x[:,0],a[0],a[1],a[2],a[3])
        reg = np.polyfit(x[:,0],y[:,0], 7)
        P_hat_vector = np.polyval(reg, x[:,0])                 
        for i in range (0,N):
            if (np.maximum(K-stprice_matrix[i,j],0)) >= P_hat_vector[i]:
                gk[i,0] = np.maximum(K-stprice_matrix[i,j],0)
                tau_k[i,0]=j
            else:
                continue

    P_hat_0 = 0
    for i in range (0,N):
        P_hat_0 = P_hat_0 + np.exp(-r*tau_k[i,0]*delta_t)*gk[i,0]
    P_hat_0 = P_hat_0/N
    timezeroTrueVal = K-stprice_matrix[0,0]
    V_0 = np.maximum(timezeroTrueVal,P_hat_0)
    
    return ("%.6f" % V_0)   
    
