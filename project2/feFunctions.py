import numpy as np
import matplotlib.pyplot as plt

#Function definitions

def initialConditions(x):
    '''
    This function defines the initial conditions of the system at u(x,0).
    '''
    return np.sin(np.pi*x)

def forcingFunction(x,t):
    '''
    This function defines the forcing function of the ODE.
    '''
    return (np.pi**2 - 1) * np.exp(-1 * t) * np.sin(np.pi*x)

def mappedForcingFunction(x1,x2,x,t):
    '''
    This function returns the forcing function when mapped to range (-1, 1).
    '''
    xMapped = ((x2-x1)/2)*(x + 1) + x1
    return forcingFunction(xMapped,t)

def lagrange1(x):
    return (1-x)/2

def lagrange2(x):
    return (1+x)/2

def analyticalSolution(x,t):
    '''
    This function defines the analytical solution of the ODE.
    '''
    return np.exp(-1*t) * np.sin(np.pi*x)

def quadrature(x1,x2,t):
    '''
    This function returns f(x) * v(x) integrated over (-1, 1) using 2nd-order Gaussian quadrature.

    Args:
        x1 (float): position of node 1
        x2 (float): position of node 2
        t (float): timestep
    Returns:
        fLocal (array): 2x1 matrix containing f(x) * v(x) integrated at both nodes
    '''
    #Define 2nd-order quadrature points
    q = [(1/np.sqrt(3)), (-1/np.sqrt(3))]
    w = [1,1]
    
    dxdxi = (x2-x1)/2

    int1 = 0
    for i in range(len(q)):
        int1 += w[i] * mappedForcingFunction(x1,x2,q[i],t) * lagrange1(q[i]) * dxdxi
    int2 = 0
    for i in range(len(q)):
        int2 += w[i] * mappedForcingFunction(x1,x2,q[i],t) * lagrange2(q[i]) * dxdxi

    return [int1, int2]

def globalAssemble(G,l,i):
    '''
    This function adds local matrices to the global M or K matrix.

    Args:
        G (2D array): global matrix
        l (2D array): local matrix
        i (int): index of element
    Returns:
        G (2D array): global matrix with local matrix i added in
    '''
    G[i][i] += l[0][0]
    G[i][i+1] += l[0][1]
    G[i+1][i] += l[1][0]
    G[i+1][i+1] += l[1][1]

    return G