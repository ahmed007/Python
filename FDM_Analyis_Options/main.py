# -*- coding: utf-8 -*-
"""

Course Name: Numerical Methods for Financial Derivatives
Course Number=6204
Assignment Description:
Part1 
I am asked to implement Finite Difference Methods to simulate European and 
American options price.The implemented Finite Difference methods include Explicit
FDM,Implicit FDM and Thomas algorithm[ European Only] Implicit FDM and 
Brennan-Schwartz algorithm (for American options) Implicit FDM and SOR 
(for European options), Implicit FDM and Projected SOR (for American 
options), Crank-Nicholson and Thomas algorithm (for European options), Crank-
Nicholson and Brennan-Schwartz algorithm (for American options),Crank-Nicholson
FDM and SOR (for European options), Crank-Nicholson FDM and Projected SOR 
(for American options).

Part2 : Monte Carlo integration of risk-neutral expectation (European options),
 Regression method II (for American options).
 
Part3: Closed-form solution formulas (for European options BS )
This is the main file from where i instantiate classes , initialize algorithmic
 parameters.
 
Total Number of files in project = 11

Dependent Files: thomasAlgo.py , SORsolver.py, RegressionMethod.py, projectedSOR.py
MonteCarloEulers.py, mcClass.py, FDM.py, EuropeanClosedForm.py, closedFormClass.py,
brennanschwartz.py

Important Note: To change any algorothmic parameter only chnage the values in
main function.
"""

import sys
import time
import numpy as np
import pandas as pd
from mcClass import monteCarlo
from closedFormClass import blackScholes
from FDM import FDM
import time
from threading import Thread

"""
Function Name:absoluteError
argu: truevalue from BS closed form a solution,estimated value
returns the absolute value of the error

"""
def absoluteError(trueValue,estimatedValue):
    retVal=[]
    indx=0
    
    for v in trueValue:
        retVal.append(np.absolute(float(v) - float(estimatedValue[indx])))
        indx = indx+1
        
    return retVal
    

"""
function name : MonteCarloExplicitEulers
argu: object to class monteCarlo
interface function between mcClass and driver file gives the call to member 
function of class mcClass to calculate European option price.
"""
def MonteCarloExplicitEulers(monteCarlo,priceList):
    
    print(" Please wait: Calculating Monte Carlo Explicit Eulers European:")

    Ecallprice,Eputprice=monteCarlo.Euler()
    Value=[Ecallprice,Eputprice]
    abosluteError=absoluteError(priceList,Value)
    print(" ")
    print(" ===============Monte Carlo explicit Eulers O/P===================================== ")
    df=pd.DataFrame([[Value[0],Value[1],abosluteError[0],abosluteError[1]]],columns=['European Call','European Put','Call Abs error','Put Abs error'],\
                    index = ['  1. Option Price']) 
    print(df)
  
    print(" ==================================================================================== ")
    
"""
function name : EuropeanClosed
argu: object to class blackScholes
interface function between class and driver file
interface function between closedFormClass and driver file gives the call to member 
function of class closedFormClass to calculate European option price.
"""
def EuropeanClosed(blackScholes):
    
    Ecallprice,Eputprice=blackScholes.closedFormInterface()
    Value=[Ecallprice,Eputprice]
    print(" ")
    print(" ===============Black Scholes O/P================================================== ")
    df=pd.DataFrame([Value],columns=['European Call','European Put'],\
                    index = ['  1. Option Price']) 
    print(df)
  
    print(" ==================================================================================== ")   
    
    return Value
    
"""
function name : RegressionMethod2
argu: object to class monteCarlo
interface function between mcClass and driver file gives the call to member 
function of class mcClass to calculate American option price.
"""    
def RegressionMethod2(monteCarlo):
   
    american_call,american_put=monteCarlo.regressionTwoInterface()
    Value=[american_call,american_put]
    print(" ")
    print(" ===============Regression 2 O/P=================================================== ")
    df=pd.DataFrame([[Value[0],Value[1]]],columns=['American Call','American Put'],\
                    index = ['  1. Option Price']) 
    print(df)
  
    print(" ==================================================================================== ")    
   
"""
function name : FiniteDifferenceMethods
arg1: object to class FDM
arg2: val indicating which method to simulate first
interface function between FDM class and main file, uses numerical methods 
explicit,implict, and Cranknilcolson to value European and AMerican option
"""
def FiniteDifferenceMethods(FDM,val,priceList) :
    Ecallprice = 0.0
    Eputprice=0.0
    Amcallprice = 0.0
    Amputprice = 0.0
    varStr=""
    
    if val == "11" :
        varStr ="Value of European Option using Explicit Method "
        Ecallprice = FDM.calculateOptionPrice(method="explicit",option="call",algo="none")
        Eputprice = FDM.calculateOptionPrice(method="explicit",option="put",algo="none")
    elif val == "12" :
        varStr ="Value of American Option using Explicit Method "
        Amcallprice = FDM.calculateOptionPrice(method="explicit",option="amc",algo="none")
        Amputprice = FDM.calculateOptionPrice(method="explicit",option="amp",algo="none")
    elif val == "13" :
        varStr ="Value of European Options using Implicit Method algo=thomas "
        Ecallprice = FDM.calculateOptionPrice(method="implicit",option="call",algo="thomas")
        Eputprice = FDM.calculateOptionPrice(method="implicit",option="put",algo="thomas")        
    elif val == "14" :
        varStr ="Value of American Options using Implicit Method algo=Brennan-Schwartz"
        Amcallprice = FDM.calculateOptionPrice(method="implicit",option="amc",algo="BSCh")
        Amputprice = FDM.calculateOptionPrice(method="implicit",option="amp",algo="BSCh")           
    elif val == "15" :
        varStr ="Value of European Options using Implicit Method algo=SOR"
        Ecallprice = FDM.calculateOptionPrice(method="implicit",option="call",algo="SOR")
        Eputprice = FDM.calculateOptionPrice(method="implicit",option="put",algo="SOR")             
    elif val == "16" :
        varStr ="Value of American Options using Implicit Method algo=PSOR"
        Amcallprice = FDM.calculateOptionPrice(method="implicit",option="amc",algo="PSOR")
        Amputprice = FDM.calculateOptionPrice(method="implicit",option="amp",algo="PSOR")          
    elif val == "17" :
        varStr ="Value of European Options using Implicit Method algo=thomas"
        Ecallprice = FDM.calculateOptionPrice(method="crank",option="call",algo="thomas")
        Eputprice = FDM.calculateOptionPrice(method="crank",option="put",algo="thomas")           
    elif val == "18" :
        varStr ="Value of American Options using Crank Method algo=thomas"
        Amcallprice = FDM.calculateOptionPrice(method="crank",option="amc",algo="BSCh")
        Amputprice = FDM.calculateOptionPrice(method="crank",option="amp",algo="BSCh")      
    elif val == "19" :
        varStr ="Value of European Options using Crank Method algo=SOR"
        Ecallprice = FDM.calculateOptionPrice(method="crank",option="call",algo="SOR")
        Eputprice = FDM.calculateOptionPrice(method="crank",option="put",algo="SOR") 
    elif val == "20" :
        varStr ="Value of American Options using Crank Method algo=PSOR"
        Amcallprice = FDM.calculateOptionPrice(method="crank",option="amc",algo="PSOR")
        Amputprice = FDM.calculateOptionPrice(method="crank",option="amp",algo="PSOR") 
    
    if Amputprice == 0.0 and Amcallprice == 0.0:
        Value=[Ecallprice,Eputprice]
        abosluteError=absoluteError(priceList,Value)
        print(" ")
        print(" ==============="+varStr+" O/P========================================== ")
        df=pd.DataFrame([[Value[0],Value[1],abosluteError[0],abosluteError[1]]],columns=['European Call','European Put','abs_call_err','abs_put_err'],\
                    index = ['  1. Option Price']) 
        print(df)
  
        print(" ==================================================================================== ")  
    else:
        Value=[Amcallprice,Amputprice]
        print(" ")
        print(" ==============="+varStr+" O/P========================================== ")
        df=pd.DataFrame([Value],columns=['American Call','American Put'],\
                    index = ['  1. Option Price']) 
        print(df)
  
        print(" ==================================================================================== ")
        

"""
main function from where the execution starts
    S0=100 #stock price
    K=100 #Strike Price
    mu=0.05
    r=0.02 #risk free int rate
    delta=0.01 #dividend
    delta_T = 0.00125 # for MC simulation
    sigma=0.6 #volatility
    T=1
    t=0    
    N=500 #paths
    seed=10   #for GBM
    delta_x = 0.05 #the space step in FDM
    delta_tau = 0.00125 #for explicit implicit simulation
    xmin = -2.5 #negative infinite x domain value
    xmax = 2.5 #positive infinite x domain value
    tol=10**-6 #absolute error tolerance
    omega = 1.10 #relaxxation parameter
"""       
        
if __name__=='__main__':
    "Algorithmic Parameters "
    S0=100 #stock price
    K=100 #Strike Price
    mu=0.05
    r=0.02 #risk free int rate
    delta=0.01 #dividend
    delta_T = 0.00125 # for MC simulation
    sigma=0.6 #volatility
    T=1
    t=0    
    N=500 #paths
    seed=123   #for GBM
    delta_x = 0.05 #the space step in FDM
    delta_tau = 0.00125 #for explicit implicit simulation
    xmin = -2.5 #negative infinite x domain value
    xmax = 2.5 #positive infinite x domain value
    tol=10**-6 #absolute error tolerance
    omega = 1.10 #relaxxation parameter
    
    print("initializing algorithmic parameters for MonteCarlo,BlackScholes and FDM....")
    monteCarlo  = monteCarlo(S0,K,r,delta,delta_T,sigma,T,t,N,seed)
    blackScholes = blackScholes(S0,K,r,delta,sigma,T,t)
    #S0, K, r, dividend,sigma,T,t,delta_x,delta_tau,xmin,xmax)
    FDM = FDM(S0,K,r,delta,sigma,T,t,delta_x,delta_tau,xmin,xmax,omega,tol)    
    print("initialization done..")

    priceList = EuropeanClosed(blackScholes)
    MonteCarloExplicitEulers(monteCarlo,priceList)
    RegressionMethod2(monteCarlo)
    
    # i am using below list to iterate through explicit and implicit combinations               
    methodsArray = ["11","12","13","14","15","16","17","18","19","20"]
    for val in methodsArray :
        FiniteDifferenceMethods(FDM,val,priceList)
    
    
    