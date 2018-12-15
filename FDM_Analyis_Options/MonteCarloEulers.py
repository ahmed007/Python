# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 20:13:09 2018

@author: admin
"""
import numpy as np
import math
from scipy.stats import norm 
import pandas as pd


"""
Function name: EulersDiscretization
Description: This function calculates stock price matrix using Euler's scheme
Inp  parameters: Strike price,Stock price, risk free rate(r),dividend(delta),timestep size=delta_t
N->Number of paths T-> Expiry date
Return Value: Stock Price matrix
"""

def EulersDiscretization(S0,K,r,delta,sigma,delta_t,N,T,seed):
    stprice_matrix = np.zeros((N,int(T/delta_t)+1),dtype=float)
    stprice_matrix[:,0] = S0    
    np.random.seed(seed)
    Brownian_matrix=np.random.standard_normal(size=(N,int(T/delta_t) + 1))*np.sqrt(delta_t)
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
Function name: MilsTeinDiscretization
Description: This function calculates stock price matrix using Milsteing scheme
Inp  parameters: Strike price,Stock price, risk free rate(r),dividend(delta),timestep size=delta_t
N->Number of paths T-> Expiry date
Return Value: Stock Price matrix
"""
def MilsTeinDiscretization(S0,K,mu,r,delta,sigma,delta_t,N,T,seed):
    stprice_matrix = np.zeros((N,int(T/delta_t)+1),dtype=float)
    stprice_matrix[:,0] = S0    
    np.random.seed(seed)
    Brownian_matrix=np.random.standard_normal(size=(N,int(T/delta_t)))*np.sqrt(delta_t)
    for j in range(0,N):
        for i in range(1,int(T/delta_t)):
            DW=Brownian_matrix[j][i]
            #DW=np.sum(Brownian_matrix[j][(1*(i-1)+1):(1*i + 1)])
            var = stprice_matrix[j][i-1] + (r-delta)*stprice_matrix[j][i-1]*delta_t+sigma*stprice_matrix[j][i-1]*DW\
            + 0.5 * sigma**2 * delta_t * (DW**2 - 1)
            #prevdw=DW
            stprice_matrix[j][i] = var
            #St.append(var)
            
    return stprice_matrix

"""
Function name: computeVCEur
Description: Function to calculate European Call price using stock matrix 
returned from Euler's or Milstein

Return Value: time zero price of Call Option
"""
def computeVCEur(stock_matrix,K,r,N,T,delta_T):
    maturity_index=int(T/delta_T)
    VCEur_ST_k=0
#    for k in range(0,N):
#        VCEur_ST_k = VCEur_ST_k + np.maximum((stock_matrix[k][maturity_index]-K),0)
    VCEur_ST_k = np.maximum((stock_matrix[:,maturity_index] - K),0)
    
    VC_Eur = np.exp(-r*T) * np.sum(VCEur_ST_k)/N
    return ("%.6f" %VC_Eur)

"""
Function name: computeVCEur
Description: Function to calculate European Put price using stock matrix 
returned from Euler's or Milstein

Return Value: time zero price of Put Option
"""	
def computeVPEur(stock_matrix,K,r,N,T,delta_T):
    maturity_index=int(T/delta_T)
    VPEur_ST_k=0
#    for k in range(0,N):
#        VPEur_ST_k = VPEur_ST_k + np.maximum((K-stock_matrix[k][maturity_index]),0)
    VPEur_ST_k = np.maximum((K-stock_matrix[:,maturity_index]),0) 
    VP_Eur = np.exp(-r*T) * np.sum(VPEur_ST_k)/N 
    return ("%.6f" %VP_Eur)   
        
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
    
    European_Call_price=[]
    
    for sig in sigma:        
        d1 = (math.log(S0/K) + (r - delta + sig**2/2)*(T-t))/(sig*math.sqrt(T-t))
        d2 = d1 - sig*math.sqrt(T-t)
        European_Call_price.append(S0*math.exp((-delta)*(T-t))*F(d1) - K*math.exp((-r)*(T-t))*F(d2))
        
    
    return European_Call_price

"""
Function name: BS_European_Put_Price
Description: This function calculates BS European put option price at time t
Inp  parameters: Strike price,Stock price, t,T,delta,interest rates, sigma
Return Value: List containing the European call option price
"""
def BS_European_Put_Price(K,S0,t,T,delta,r,sigma):
    
    European_Put_price=[]
    
    for sig in sigma:        
        d1 = (math.log(S0/K) + (r - delta + sig**2/2)*(T-t))/(sig*math.sqrt(T-t))
        d2 = d1 - sig*math.sqrt(T-t)
        European_Put_price.append(-1*S0*math.exp((-delta)*(T-t))*F(-d1) + K*math.exp((-r)*(T-t))*F(-d2))
        
    
    return European_Put_price    