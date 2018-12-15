# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 09:01:02 2018
"""
from EuropeanClosedForm import BS_European_Call_Price
from EuropeanClosedForm import BS_European_Put_Price
"""
class blackScholes this class constructor initializes parameters needed in calculating
European option prices.
"""
class blackScholes:
    def __init__(self, S0, K, r, dividend,sigma,T,t):
        self.S0=S0
        self.K=K
        self.r=r
        self.delta=dividend
        self.sigma=sigma
        self.T=T
        self.t=t
    #from here i make call to closed form formula implement in file EuropeanClosedForm
    def closedFormInterface(self):
        Ecallprice=BS_European_Call_Price(self.K,self.S0,self.t,self.T,self.delta,self.r,self.sigma)
        Eputprice=BS_European_Put_Price(self.K,self.S0,self.t,self.T,self.delta,self.r,self.sigma)
        
        return Ecallprice,Eputprice
        
