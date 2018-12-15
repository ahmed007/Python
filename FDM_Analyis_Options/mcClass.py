# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 08:26:59 2018

@author: admin
"""
from EuropeanClosedForm import BS_European_Call_Price
from EuropeanClosedForm import BS_European_Put_Price
from MonteCarloEulers import EulersDiscretization
from MonteCarloEulers import computeVCEur
from MonteCarloEulers import computeVPEur
from RegressionMethod import regression_method2_call
from RegressionMethod import regression_method2_put
from RegressionMethod import EulersDiscretizationRegression

"""
class monteCarlo this class constructor initializes parameters needed in calculating
European option prices.
"""
class monteCarlo:
    def __init__(self, S0, K, r, dividend,delta_T,sigma,T,t,N,seed):
        self.S0=S0
        self.K=K
        self.r=r
        self.delta=dividend
        self.delta_T=delta_T
        self.sigma=sigma
        self.T=T
        self.t=t
        self.N=N
        self.seed=seed
    #call to Euler scheme define in file MonteCarloEulers
    def Euler(self):
            stmatrix = EulersDiscretization(self.S0,self.K,self.r,self.delta,self.sigma,self.delta_T,self.N,self.T,self.seed)
            Ecallprice=computeVCEur(stmatrix,self.K,self.r,self.N,self.T,self.delta_T)
            Eputprice=computeVPEur(stmatrix,self.K,self.r,self.N,self.T,self.delta_T)
            
            return Ecallprice , Eputprice
        
    #call to Regression scheme define in file RegressionMethod   
    def regressionTwoInterface(self):
        stmatrix = EulersDiscretizationRegression(self.S0,self.K,self.r,self.delta,self.sigma,self.delta_T,self.N,self.T,self.seed)
        american_call=regression_method2_call(stmatrix,self.K,self.T,self.delta_T,self.N,self.r)
        stmatrix = EulersDiscretizationRegression(self.S0,self.K,self.r,self.delta,self.sigma,self.delta_T,self.N,self.T,self.seed)
        american_put=regression_method2_put(stmatrix,self.K,self.T,self.delta_T,self.N,self.r)
        
        return american_call,american_put
        
        