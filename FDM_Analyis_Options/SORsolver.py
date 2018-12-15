# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 11:54:23 2018

"""


import sys
import numpy as np
from scipy.sparse import diags
from numpy.linalg import eigh
from numpy.linalg import inv

"""
function Name:createTriDiagonalMatrix
input parameters: alpha_0 main diag element,beta_0 upper diag element,
gamma_1 lower diag element, N size of square matrix diagonally.
retruns Tridiagonal matrix
"""
def createTriDiagonalMatrix(alpha_0,beta_0,gamma_1,N):
    
    
    k = np.array([gamma_1*np.ones(N-1,dtype=float),alpha_0*np.ones(N,dtype=float),beta_0*np.ones(N-1,dtype=float)])        
    offset = [-1,0,1]
    G = diags(k,offset).toarray()
    
    return G 

"""
SOR method was devised simultaneously by David M. Young, Jr. and by Stanley P. Frankel 
in 1950 for the purpose of automatically solving linear systems on digital computers. 
Over-relaxation methods had been used before the work of Young and Frankel.(wiki).
Given a square system of n linear equations with unknown x: A x = b 
Iterative Methods for Linear Systems of Equations
Input Parameters: A tridiagonal matrix, X0 initial guess vector, itera number of
iteration , omega is overrelaxation parameter, tolerance is the tolerance level.

"""
def SORsolver(A,X0,bvector,itera,omega,tolerance):
    
    Nrows = len(A)
    NCol = len(A[0])
    
    
    if Nrows != NCol:
        print("The A vector is not a square matrix")
        print("Exiting program")
        sys.exit()
   #A=D+L+U
    D=np.diag(np.diag(A))
    L=np.tril(A) - D #tri_lower_diag = np.tril(a, k=0)
    U=np.triu(A)-D #tri_upper_diag = np.tril(a, k=0)
    mat=np.zeros((Nrows,NCol),dtype=float)
    mat=inv(D+omega*L) * (D*(1-omega)-omega*U)
    vals, vecs=eigh(mat)
    if np.abs(max(vals)) >=1:
        print("Since the modulus of largest eigen value of iterative matrix is not less than 1")
        print("This process is not convergent. Please try some other process.")
        sys.exit()
    
    X=np.zeros((Nrows,itera),dtype=float)
    X[:,0] = X0[:,0]
    
    for k in range(0,itera-1):
        temp=0
        temp1=0
        for i in range (0,Nrows) :
            
            for j in range (0,i):
                temp=temp+A[i,j]*X[j,k+1]
            for j in range(i+1,Nrows):
                temp1=temp1+A[i,j]*X[j,k]
            X[i,k+1] = (1-omega)*X[i,k] + omega/A[i,i] *( bvector[i] - temp - temp1 )
            temp=0
            temp1=0
        if k > 0:
            if np.abs(np.sum(X[:,k+1] - X[:,k])) < tolerance :
                break;   
            
            
    return X[:,k]         
        