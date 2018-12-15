# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 10:47:34 2018

"""
import sys
import numpy as np
from scipy.sparse import diags
from numba import jit    

"""
function name=createTriDiagonalMatrix
descrip: this function creates tridiagonal matrix.
parameters:  alpha_0 main diagonal element
beta_0 upper diagonal element
gamma_1 lower diagonal elements
N square matrix size
returns G a tridiagonal matrix
"""

def createTriDiagonalMatrix(alpha_0,beta_0,gamma_1,N):
    
    
    k = np.array([gamma_1*np.ones(N-1,dtype=float),alpha_0*np.ones(N,dtype=float),beta_0*np.ones(N-1,dtype=float)])        
    offset = [-1,0,1]
    G = diags(k,offset).toarray()
    
    return G    

"""
function name=applyThomas
descrip: interface function which is visible to end user

"""
def applyThomas(matrix,bvector):
    return algoThomas(matrix,bvector)
    
"""
In numerical linear algebra, the tridiagonal matrix algorithm, also known as the 
Thomas algorithm (named after Llewellyn Thomas), is a simplified form of 
Gaussian elimination that can be used to solve tridiagonal systems of equations.
function name=algoThomas
descrip: this function implements Thomas algorithm.
parameters:  matrix -> Tridiagonal matrix, bvector->RHS value of equation \
A.x = b
"""   
@jit
def algoThomas(matrix,bvector):
    
    rows=len(matrix)
    col=len(matrix[0])
    lowd_gm = np.diag(matrix,-1) #lower diag
    maind_alph = np.diag(matrix,0) #main diag
    uppd_beta = np.diag(matrix,1) #upper diag
     
    
    if rows != col:
        print("Not a sqaure matrix")
        sys.exit()
    alpha0_hat = np.zeros(rows,float)
    b_hat =np.zeros(rows,float)
    solvect_x =np.zeros(rows,float)
    

    
    alpha0_hat[0] = matrix[0][0]
    b_hat[0] = bvector[0]
    
    #step1 forward loop
    for i in range(1,rows):
        alpha0_hat[i] = maind_alph[i] - uppd_beta[i-1] *(lowd_gm[i-1]/alpha0_hat[i-1])
        b_hat[i] = bvector[i] - b_hat[i-1] *(lowd_gm[i-1]/alpha0_hat[i-1])
    
    #step 2 bachward loop  
    solvect_x[rows-1] = b_hat[rows-1]/alpha0_hat[rows-1]
    for j in range (rows-2,-1,-1):
        solvect_x[j] = (b_hat[j] - uppd_beta[j] * solvect_x[j+1])/alpha0_hat[j]
        
    return solvect_x

"""
unused function kept it for reference for future. source stackoverflow
function name=algoThomas
descrip: this function implements Thomas algorithm.
parameters:  matrix -> Tridiagonal matrix, bvector->RHS value of equation \
A.x = b
""" 
def TDMA(a,b,c,d):
    n = len(d)
    w= np.zeros(n-1,float)
    g= np.zeros(n, float)
    p = np.zeros(n,float)

    w[0] = c[0]/b[0]
    g[0] = d[0]/b[0]

    for i in range(1,n-1):
        w[i] = c[i]/(b[i] - a[i-1]*w[i-1])
    for i in range(1,n):
        g[i] = (d[i] - a[i-1]*g[i-1])/(b[i] - a[i-1]*w[i-1])
    p[n-1] = g[n-1]
    for i in range(n-1,0,-1):
        p[i-1] = g[i-1] - w[i-1]*p[i]
    return p  
    