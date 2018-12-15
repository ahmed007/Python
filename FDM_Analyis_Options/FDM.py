# -*- coding: utf-8 -*-
"""
This file contains code for FDM implicit and explicit method.
This file handles all transformation and setting up boundary conditions
for both European and American options.
Reference Book: Tools for Computational finance 
"""
import sys
import numpy as np
from scipy.sparse import diags
from thomasAlgo import applyThomas
from brennanschwartz import applyBrennanScwartz
from SORsolver import SORsolver
from projectedSOR import projectedSOR

class FDM:
    def __init__(self, S0, K, r, dividend,sigma,T,t,delta_x,delta_tau,xmin,xmax,omega,tol):
        self.S0=S0
        self.K=K
        self.r=r
        self.delta = dividend
        self.T = T #maturity time
        self.t = t
        self.sigma = sigma #vol
        self.delta_x = delta_x #on grid delta_x represents mesh x axis spacing
        self.delta_tau = delta_tau #on grid delta_tau represents mesh time interval
        self.a=xmin #grid xmin this is to fix domain -inifity transform to xmin
        self.b=xmax #grid xmax this is to fix domain +inifity transform to xmax
        #black scholes eqation is transformed into equivalent heat equation below are transformations
        self.q_delta=2*(self.r-self.delta)/self.sigma**2
        self.q = 2*(self.r)/self.sigma**2
        self.m = int((self.b - self.a)/self.delta_x) #total number of points on x axix
        self.tau_max = 0.5 * (self.sigma**2) * self.T
        self.v_max = int(self.tau_max/self.delta_tau) #total number of points on Tau axix
        self.tau_v_mat=np.zeros((1,self.v_max+1),dtype=float)
        self.xmat = np.zeros((1,self.m+1),dtype=float)
        self.w_matrix = self.__utilgetShape(self.xmat,self.tau_v_mat)
        self.v_matrix = self.__utilgetShape(self.xmat,self.tau_v_mat)
        self.V_vector = np.zeros((1,self.m+1),dtype=float)
        self.g_matrix = self.__utilgetShape(self.xmat,self.tau_v_mat) #for american options
        self.sor_omega = omega
        self.sor_tol = tol #tolerance level for SOR PSOR
        self.sor_iter = 200
        self.theta = 1 #default value would be overridden in code below
        
    """
    Funnction: __utilgetShape 
    Util function to get shape
    
    """    
    def __utilgetShape(self,tau,xaxis):
        row,col = tau.shape
        row1,col1 = xaxis.shape
        return np.zeros((col,col1),dtype=float)
     
    """
    Funnction: Transformation 
    This defines a two-dimensional uniform grid,
    equidistant grid in the book is defined in terms of x and tau
    Reference: page 185 Seydel
    """        
    def __transformation(self,optionType):

        for v in range(0,self.v_max+1):
            self.tau_v_mat[0,v] = v*self.delta_tau
        for i in range(0,self.m+1):
            self.xmat[0,i] = self.a + i*self.delta_x
        
        
    """
    Funnction: __initialConditionBoundry private function
    This function handles only the intial boundary conditions for American and 
    plain vanilla options.
    Reference: page 182 Seydel equation 4.5 and 4.6 for european. W_matrix is the
    main grid, g_matrix contains value for fixed tau.
    American initial boundary condition :w_matrix[x][0] = g[x][0] and w_matrix[0][tau]=g[0][tau]
    Reference:equation 4.41 Seydel
    """    
   #self.call_matrix.shape = (145,101)=(xmat,tau)
    def __initialConditionBoundry(self,optionType):
        j=0
        if optionType == "call" :
            for i in self.xmat[0,:]:
                self.w_matrix[j,0] = max(np.exp(i/2*(self.q_delta+1))-np.exp(i/2*(self.q_delta-1)),0)
                j=j+1
#            self.put_matrix[j,0] =  max(np.exp(i/2*(self.q_delta-1))-np.exp(i/2*(self.q_delta+1)),0)
            for v in range(0,self.v_max+1):
                self.w_matrix[self.m,v] = np.exp(0.5*(self.q_delta+1)*self.b + 0.25*(self.q_delta+1)**2 * self.tau_v_mat[0,v])
        elif optionType == "put" :
            for i in self.xmat[0,:]:
                self.w_matrix[j,0] = max(np.exp(i/2*(self.q_delta-1))-np.exp(i/2*(self.q_delta+1)),0)
                j=j+1            
            for v in range(0,self.v_max+1):
                self.w_matrix[0,v] = np.exp(0.5*(self.q_delta-1)*self.a + 0.25*(self.q_delta-1)**2 * self.tau_v_mat[0,v])
        elif optionType == "amp" :
            for v in range(0,self.v_max+1): #i is x, v is tau
                for i in range(0,self.m+1):
                    self.g_matrix[i,v] = np.exp(self.tau_v_mat[0,v]/4*((self.q_delta-1)**2 + 4*self.q)) * \
                    max(np.exp(0.5*self.xmat[0,i]*(self.q_delta-1))-np.exp(0.5*self.xmat[0,i]*(self.q_delta+1)),0)            
#            for v in self.tau_v_mat[0,:]: #i is x v is tau
#                for i in self.xmat[0,:]:
#                    self.g_vector[v,i] = np.exp(v/4*((self.q_delta-1)**2 + 4*self.q)) * \
#                    max(np.exp(0.5*i*(self.q_delta-1))-np.exp(0.5*i*(self.q_delta+1)),0)
            for i in range (1,self.m):
                self.w_matrix[i,0]=self.g_matrix[i,0]
            for v in range (1,self.v_max+1):
                self.w_matrix[0,v]=self.g_matrix[0,v]  
            for v in range (1,self.v_max+1):
                self.w_matrix[self.m,v]=self.g_matrix[self.m,v] 
        elif optionType == "amc" :
            for v in range(0,self.v_max+1): #i is x v is tau
                for i in range(0,self.m+1):
                    self.g_matrix[i,v] = np.exp(self.tau_v_mat[0,v]/4*((self.q_delta-1)**2 + 4*self.q)) * \
                    max(np.exp(0.5*self.xmat[0,i]*(self.q_delta+1))-np.exp(0.5*self.xmat[0,i]*(self.q_delta-1)),0)            
#            for v in self.tau_v_mat[0,:]: #i is x v is tau
#                for i in self.xmat[0,:]:
#                    self.g_vector[v,i] = np.exp(v/4*((self.q_delta-1)**2 + 4*self.q)) * \
#                    max(np.exp(0.5*i*(self.q_delta-1))-np.exp(0.5*i*(self.q_delta+1)),0)
            for i in range (1,self.m):
                self.w_matrix[i,0]=self.g_matrix[i,0]
            for v in range (1,self.v_max+1):
                self.w_matrix[0,v]=self.g_matrix[0,v]  
            for v in range (1,self.v_max+1):
                self.w_matrix[self.m,v]=self.g_matrix[self.m,v]                 
                
    """
    function name=createTriDiagonalMatrix.
    descrip: this function creates tridiagonal matrix needed in Crank Nilcolsol,
    Explicit and implicit Method.
    parameters:  alpha_0 main diagonal element.
    beta_0 upper diagonal element.
    gamma_1 lower diagonal elements.
    N square matrix size.
    returns G a tridiagonal matrix.
    """

    def __createTriDiagonalMatrix(self,alpha_0,beta_0,gamma_1,N):
    
    
        k = np.array([gamma_1*np.ones(N-1,dtype=float),alpha_0*np.ones(N,dtype=float)\
                      ,beta_0*np.ones(N-1,dtype=float)])        
        offset = [-1,0,1]
        G = diags(k,offset).toarray()
    
        return G    
               
    """
    function name:getARMatrix,Explicit method need AR matrix
    parameters:  alpha_0 main diagonal element
    beta_0 upper diagonal element
    gamma_1 lower diagonal elements
    N square matrix size    
    """
       
    def getARMatrix(self,theta):
        
        lamda = self.delta_tau/self.delta_x**2
        alpha_0 = 1-2*(1-theta)*lamda
        beta_0 = (1-theta)*lamda
        gamma_0 = (1-theta)*lamda
        N = self.m -1
        AR_mat=self.__createTriDiagonalMatrix(alpha_0,beta_0,gamma_0,N)
       
        return AR_mat

    """
    function name:getALMatrix.
    Implicit method need AL matrix on LHS.
    parameters:  alpha_0 main diagonal element.
    beta_0 upper diagonal element.
    gamma_1 lower diagonal elements.
    N square matrix size    
    """
        
    def getALMatrix(self,theta):
        
        lamda = self.delta_tau/self.delta_x**2
        alpha_0 = 1+2*theta*lamda
        beta_0 = -theta*lamda
        gamma_0 = -theta*lamda
        N = self.m - 1
        AL_mat=self.__createTriDiagonalMatrix(alpha_0,beta_0,gamma_0,N)
       
        return AL_mat    
            
    """
    function name:explicitMethod  
    parameters:  option="call or put"
    returns option price to callee
    """        
    def explicitMethod(self,option="call or put"):
        
        
        self.__transformation(option)
        self.__initialConditionBoundry(option)
        
        lamda = self.delta_tau/self.delta_x**2
        
#        for i in range(0,self.m+1):
#            self.w_matrix[i,0] = self.call_matrix[i,0]
        
        for v in range(0,self.v_max):
            for i in range(1,self.m):
                self.w_matrix[i,v+1] = lamda*self.w_matrix[i-1,v] + (1-2*lamda)*self.w_matrix[i,v]\
                +lamda*self.w_matrix[i+1,v]
                
        for v in range(0,self.v_max+1):
            for i in range(0,self.m+1):
                self.v_matrix[i,v] = self.K * np.exp(-0.5*(self.q_delta-1)*self.xmat[0,i]\
                - (0.25*(self.q_delta-1)**2+self.q)*self.tau_v_mat[0,v]) * self.w_matrix[i,v]
                
        return self.__getoptionPrice()

    """
    function name:__retriveSubMatrix_w  
    parameters:  indx
    Decription: This function retrieve sub matrix w from the main grid needed by
    the implicit method function.This submatrix excludes boundary values i.e in
    my case top most row, bottom most row and leftmost column values are excluded
    and i take values in between.
    This is helper/utility function
    """   
    def __retriveSubMatrix_w(self,indx):
        sub_w_matrix = np.zeros((self.m-1,1),dtype=float)
        
        for i in range(1,self.m):
            sub_w_matrix[i-1,0] = self.w_matrix[i,indx]
            
        return sub_w_matrix

    """
    function name:__retriveSubMatrix_g  
    parameters:  indx
    Decription: This function retrieve sub matrix g from the g grid needed by
    the implicit method function,gmatrix contains initial boundary condition
    and it is initialized in function __initialConditionBoundry
    This is helper/utility function
    """    
    def __retriveSubMatrix_g(self,indx):
        sub_g_matrix = np.zeros((self.m-1,1),dtype=float)
        
        for i in range(1,self.m):
            sub_g_matrix[i-1,0] = self.g_matrix[i,indx]
            
        return sub_g_matrix    


    """
    function name:__updateSubMatrix_w  
    parameters:  solvector,indx
    Decription: This function updates matrix w with the solution vector obtained at
    each iteration. the implicit method function needs this
    This is helper/utility function
    """   
    def __updateSubMatrix_w(self,solvector,indx):
        for i in range(1,self.m):
            self.w_matrix[i,indx+1]=solvector[i-1] 
     
    """
    function name:__crankNicholsonMethod 
    Crank and Nicolson suggested to average the forward- and the backward 
    difference method.
    parameters:  optiontype, algo type
    returns option price
    """  
    def __crankNicholsonMethod(self,option="call or put",algo="thomas or BSCh"):  
        price=self.__implicitMethod(option,algo)
        return price

    """
    function name:promptError  
    parameters:  logs error if any
    """
    def promptError(self,errorType):
        print("Error::"+errorType)
        sys.exit()
        

    """
    function name:calculateOptionPrice  
    parameters:  method type, option type, algo type supported by this programme
    return option price
    """    
    def calculateOptionPrice(self,method="explicit or implicit or crank",option="call/put/amc/amp" \
                             ,algo="thomas/BSCh/SOR/PSOR") :
        retVal = 0.0;
        
        if method == "explicit" :
            self.theta=0
            algo ="none"
            retVal=self.explicitMethod(option)
        elif method == "implicit" :
            if algo != "thomas" and algo != "BSCh" and algo != "SOR" and algo != "PSOR" :
                self.promptError("calculateOptionPrice, Unsupported algo type!")
            self.theta = 1
            retVal = self.__implicitMethod(option,algo)
        elif method == "crank" :
            if algo != "thomas" and algo != "BSCh" and algo != "SOR" and algo != "PSOR" :
                self.promptError("Unsupported algo type!")  
            self.theta = 0.5
            retVal=self.__crankNicholsonMethod(option,algo)
        else:
            print("Wrong arguments specified")
            sys.exit(1)
            
        return retVal

    
    

    """
    function name:__implicitMethod.
    we try the backward difference in implicit method
    parameters:  option type and algo type
    return option price
    """    
    #private function using double underscore
    def __implicitMethod(self,option="call or put",algo="thomas/BSCh/SOR/PSOR"):

        
        self.__transformation(option)
        AL=self.getALMatrix(self.theta)
        AR = self.getARMatrix(self.theta)
        self.__initialConditionBoundry(option)
        
        rows,col = self.w_matrix.shape
        f_matrix = np.zeros((rows-2,2),dtype=float)
        lamda = self.delta_tau/self.delta_x**2
        
        for i in range(0,self.v_max):
             f_matrix[0,0] = lamda*self.w_matrix[0,i] #implicit boundary condition
             f_matrix[self.m-2,0] = lamda*self.w_matrix[self.m,i] #implicit boundary condition 
             
             f_matrix[0,1] = -lamda*self.w_matrix[0,i+1] #implicit boundary condition
             f_matrix[self.m-2,1] = -lamda*self.w_matrix[self.m,i+1] #implicit boundary condition
             
             sub_w_matrix = self.__retriveSubMatrix_w(i)
             b_vector = np.matmul(AR,sub_w_matrix[:,0]) - f_matrix[:,0] -f_matrix[:,1]
             #A.x=b thomas direct solver method x-> unknowns and b-> are RHS vector
             #generic Thomas algorithm solver below
             if algo == "thomas" : #european
                 solvector=applyThomas(AL,b_vector) 
                 
             elif algo == "BSCh" : #american
                 sub_g_matrix = self.__retriveSubMatrix_g(i)
                 solvector=applyBrennanScwartz(AL,b_vector,sub_g_matrix)
             elif algo == "SOR" : #european
                 #Any initial guess for X i am assuming all ones         
                 X0 = np.ones((rows-2,1),dtype=float)                 
                 solvector=SORsolver(AL,X0,b_vector,self.sor_iter,self.sor_omega,self.sor_tol)
             elif algo == "PSOR" :      #american  
                 #Any initial guess for X i am assuming all ones         
                 X0 = np.ones((rows-2,1),dtype=float)                  
                 sub_g_matrix = self.__retriveSubMatrix_g(i)
                 solvector=projectedSOR(AL,X0,b_vector,self.sor_iter,self.sor_omega,sub_g_matrix,self.sor_tol)
                 
             self.__updateSubMatrix_w(solvector,i)
            
            
            
        for v in range(0,self.v_max+1):
            for i in range(0,self.m+1):
                self.v_matrix[i,v] = self.K * np.exp(-0.5*(self.q_delta-1)*self.xmat[0,i]\
                - (0.25*(self.q_delta-1)**2+self.q)*self.tau_v_mat[0,v]) * self.w_matrix[i,v]       
        
        return self.__getoptionPrice()
                 
                

    """
    function name:__getoptionPrice  
    parameters:  none
    return interpolate and returns option price
    """             
                
    def __getoptionPrice(self):
        
        for i in range(0,self.m+1):
            self.V_vector[0,i] = self.K * np.exp(self.xmat[0,i])
            
        price = np.interp(self.S0, self.V_vector[0,:],self.v_matrix[:,self.v_max])
        
        return ("%.6f" %price)
        
        
                
            
            
   
        