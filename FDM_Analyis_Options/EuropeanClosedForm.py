# -*- coding: utf-8 -*-
"""
***********************************************************************
Problem Set Description: Using Black scholes Closed form formula derive the
European Call Option price at t=0, Similarly derive the European Put Option 
price at =0.

This program calculates European Call/European PUt option Price at t=0
This program has the main function, which act as the driver.

CDF used=built in function from package scipy.stats

**************************************************************************
"""

import numpy as np
import math
from scipy.stats import norm 
import pandas as pd



"""
Function name: F(x):
Description: This function calculates cdf of X
Inp  parameters: double value
Return Value: CDF of X
"""
 
def F(x):
    return norm.cdf(x)

"""
Function name: BS_European_Call_Price
Description: This function calculates BS European call option price at time t
Inp  parameters: Strike price,Stock price, t,T,delta,interest rates, sigma
Return Value: List containing the European put price for various value of sigma
"""
def BS_European_Call_Price(K,S0,t,T,delta,r,sigma):
    
    European_Call_price=0.0
       
    d1 = (math.log(S0/K) + (r - delta + sigma**2/2)*(T-t))/(sigma*math.sqrt(T-t))
    d2 = d1 - sigma*math.sqrt(T-t)
    European_Call_price=(S0*math.exp((-delta)*(T-t))*F(d1) - K*math.exp((-r)*(T-t))*F(d2))
        
    
    return ("%.6f" % European_Call_price)

"""
Function name: BS_European_Put_Price
Description: This function calculates BS European put option price at time t
Inp  parameters: Strike price,Stock price, t,T,delta,interest rates, sigma
Return Value: List containing the European call option price
"""
def BS_European_Put_Price(K,S0,t,T,delta,r,sigma):
    
    European_Put_price=0.0
        
    d1 = (math.log(S0/K) + (r - delta + sigma**2/2)*(T-t))/(sigma*math.sqrt(T-t))
    d2 = d1 - sigma*math.sqrt(T-t)
    European_Put_price=(-1*S0*math.exp((-delta)*(T-t))*F(-d1) + K*math.exp((-r)*(T-t))*F(-d2))
        
    
    return ("%.6f" % European_Put_price)





""" main function This function resembles C main and this is the point from where
the execution starts
"""
if __name__ == '__main__':
    K=100
    S0=100
    t=0
    T=1
    delta=0.025
    r=0.05
    count=0
    sigma=[0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    V_C_Eur=BS_European_Call_Price(K,S0,t,T,delta,r,sigma)
    print('S/N','   ', 'European Call option Price')
    print('----------------------------------------')
    for i in V_C_Eur:
        print(' ',count,'     ',i)
        print('------------------------------------')
        count=count+1
        
    count=0    
    V_P_Eur=BS_European_Put_Price(K,S0,t,T,delta,r,sigma)  
    
    print("    ")
    print('S/N','   ', 'European Put option Price')
    print('---------------------------------------')
    for i in V_P_Eur:
        print(' ',count,'     ',i)
        print('----------------------------------')
        count=count+1
    
