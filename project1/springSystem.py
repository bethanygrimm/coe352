from svd import InverseException, svd
import numpy as np

#Define gravitational constant
g = 9.81

invalidInput = True
while invalidInput:
    print("Enter number of masses (positive integer): ")
    try:
        nMass = int(input())
        if nMass < 1:
            pass
        else:
            invalidInput = False
    except ValueError:
        continue

invalidInput = True
while invalidInput:
    print("Enter number of fixed ends (1 or 2): ")
    try:
        fixedEnds = int(input())
        if (fixedEnds == 1) or (fixedEnds == 2):
            invalidInput = False
        else:
            pass
    except ValueError:
        continue

nSp = nMass + fixedEnds - 1

masses = []
for i in range(nMass):
    invalidInput = True
    while invalidInput:
        print("Enter mass " + str(i+1) + " (positive number): ")
        try:
            m_i = float(input())
            if m_i <= 0:
                pass
            else:
                masses.append(m_i)
                invalidInput = False
        except ValueError:
            continue

spCs = []
for i in range(nSp):
    invalidInput = True
    while invalidInput:
        print("Enter spring constant " + str(i+1) + " (positive number): ")
        try:
            k_i = float(input())
            if k_i <= 0:
                pass
            else:
                spCs.append(k_i)
                invalidInput = False
        except ValueError:
            continue

if (len(masses) != nMass) or (len(spCs) != nSp):
    print("Error: number of masses or spring constants do not match")

#Assign A, C, and AT
A1 = np.diag(np.ones(nMass))
A2 = np.diag((np.ones(nMass-1) * -1), -1)
A = A1+A2
#Start with the fixed/free matrix, and add/subtract accordingly
#If fixed/fixed, there must be an extra row, with last value -1
if fixedEnds == 2:
    A = np.append(A, np.zeros([1,nMass]), axis=0)
    A[nSp-1,nMass-1] = -1
#If free/free, we can remove the first row
if fixedEnds == 0:
    A = np.delete(A, (0), axis=0)
# print(A)

C = np.diag(spCs)
# print(C)
AT = np.transpose(A)

#Assign f = mg
f = np.array(masses) * g

#K = ATCA
K = np.matmul(np.matmul(AT, C), A)
# print(K)

#Ku = f; u = KIf
try:
    (SVD, condN, KI) = svd(K)

    u = np.matmul(KI, f)
    print("Displacement vector of masses: " + str(u))
    e = np.matmul(A, u)
    print("Elongation vector of springs: " + str(e))
    w = np.matmul(C, e)
    print("Internal stresses in springs: " + str(w))
    S = SVD[1]
    sVals = np.diagonal(S)
    eVals = sVals ** 2
    print("Singular values of K: " + str(sVals) + ", and eigenvalues of K: " + str(eVals))
    print("L2 condition number of K: " + str(condN))
except InverseException:
    print("Error: stiffness matrix K is not invertible. Try defining one or two fixed ends instead.")
