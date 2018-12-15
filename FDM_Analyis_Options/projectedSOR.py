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
Function Name :  Projected SOR [Cryer (1971)]:
parameters: Matrix A, initial guess X0, bvector holding the knowns,
omega overrelaxation parameter,gvector initial boundary condition and
tolerance level to stop execution/further processing.
Problem which we solve in PSOR w>=g, Aw-b>=0,(Aw-b)'(w-g)=0
"""

def projectedSOR(A,X0,bvector,itera,omega,gvector,tolerance):
    
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
    X[:,0] = np.maximum(gvector[:,0],X0[:,0])
    
    for k in range(0,itera-1):
        temp=0
        temp1=0
        for i in range (0,Nrows) :
            
            for j in range (0,i):
                temp=temp+A[i,j]*X[j,k+1]
            for j in range(i+1,Nrows):
                temp1=temp1+A[i,j]*X[j,k]
            X[i,k+1] = max(gvector[i,0],(1-omega)*X[i,k] + omega/A[i,i] *( bvector[i] - temp - temp1 ))
            temp=0
            temp1=0
        if k > 0:
            if np.abs(np.sum(X[:,k+1] - X[:,k])) < tolerance :
                break;   
            
            
    return X[:,k]         
        