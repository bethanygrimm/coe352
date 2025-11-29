import numpy as np
import matplotlib.pyplot as plt
from feFunctions import initialConditions, forcingFunction, mappedForcingFunction, lagrange1, lagrange2, analyticalSolution, quadrature, globalAssemble

#Variable definitions

invalidInput = True
while invalidInput:
    print("Enter number of nodes (positive integer): ")
    try:
        Nn = int(input())
        if Nn < 1:
            pass
        else:
            invalidInput = False
    except ValueError:
        continue

invalidInput = True
while invalidInput:
    print("Enter number of timesteps (positive integer): ")
    try:
        Nt = int(input())
        if Nt < 1:
            pass
        else:
            invalidInput = False
    except ValueError:
        continue

invalidInput = True
while invalidInput:
    print("Forward or backward Euler (\"f\"/\"b\"): ")
    try:
        fb = str(input())
        if fb == "f":
            forward = True
            invalidInput = False
        elif fb == "b":
            forward = False
            invalidInput = False
    except ValueError:
        continue

Ne = Nn - 1

t0 = 0
tf = 1
x0 = 0
xf = 1

#Assemble F (for t=0), M, and K
#Define the global matrices first, then integrate using quadrature or mapping

F = np.zeros((Nn, 1))
M = np.zeros((Nn, Nn))
K = np.zeros((Nn, Nn))
dx = (xf-x0)/Ne
dt = (tf-t0)/Nt
t = t0

for i in range(Ne):
    fLocal = quadrature(i*dx,(i+1)*dx,t)
    #mLocal and kLocal are constant can be directly computed
    mLocal = [[(dx/3), (1*dx/6)], [(1*dx/6), (dx/3)]]
    kLocal = [[(1/dx), (-1/dx)], [(-1/dx), (1/dx)]]

    #now assemble the global matrices
    F[i][0] += fLocal[0]
    F[i+1][0] += fLocal[1]
    M = globalAssemble(M, mLocal, i)
    K = globalAssemble(K, kLocal, i)

#Apply Dirichlet boundary conditions
#Generalize this later (?) but because u(0,t) = u(1,t) = 0, this is easy
for i in range(Nn):
    for j in range(Nn):
        if((i == 0 and j == 0) or (i==(Nn-1) and j == (Nn-1))):
            M[i][j] = 1
            K[i][j] = 1
        elif(i == 0 or j == 0 or i == (Nn-1) or j == (Nn-1)):
            M[i][j] = 0
            K[i][j] = 0
    if(i==0 or i==(Nn-1)):
        F[i][0] = 0

# print(F)
# print(M)
# print(K)

#Find u(x,0)
xRange = np.linspace(x0,xf,Nn)
xRange = xRange.reshape(Nn,1)
u = initialConditions(xRange)
uNext = initialConditions(xRange)

#Then iterate through using Euler (two methods)

if forward:
    #Forward Euler here
    #uNext = u + dt*MInv*(F-Ku)
    MInv = np.linalg.inv(M)

    for i in range(Nt - 1):
        #Recalculate F according to new t
        t += dt
        F = np.zeros((Nn, 1))
        for i in range(Ne):
            fLocal = quadrature(i*dx,(i+1)*dx,t)
            F[i][0] += fLocal[0]
            F[i+1][0] += fLocal[1]
        F[0][0] = 0
        F[Nn-1][0] = 0

        uNext = u + dt * np.matmul(MInv, (F - np.matmul(K,u)))
        u = uNext

else:
    #Backward Euler here
    #uNext = ((1/dt)*M+K)Inv * ((1/dt)*M*u + F)
    B = np.linalg.inv(M/dt + K)

    for i in range(Nt - 1):
        #Recalculate F according to new t
        t += dt
        F = np.zeros((Nn, 1))
        for i in range(Ne):
            fLocal = quadrature(i*dx,(i+1)*dx,t)
            F[i][0] += fLocal[0]
            F[i+1][0] += fLocal[1]
        F[0][0] = 0
        F[Nn-1][0] = 0
        
        uNext = np.matmul(B, ((np.matmul(M,u))/dt + F))
        u = uNext

#Plot result along with analytical solution

uResult = uNext.reshape(Nn)
if forward:
    title = "Heat Transfer Finite Element Solution with Foward Euler"
else:
    title = "Heat Transfer Finite Element Solution with Backward Euler"
xRangeA = np.linspace(x0,xf,100)
soln = analyticalSolution(xRangeA,t)
plt.plot(xRangeA,soln,label="Analytical Solution")
plt.plot(xRange,uResult,label="Finite Element Result")
plt.xlabel("Space")
plt.ylabel("Time")
plt.annotate("Spatial nodes: " + str(Nn) + "\nTimesteps: " + str(Nt), xy=(0,0.9*max(uResult)))
plt.title(title)
plt.legend()
plt.show()