import numpy as np
import numpy.typing as npt

class InverseException(Exception):
    pass

def add_evs(A: npt.ArrayLike, e_c: npt.ArrayLike, evs: npt.ArrayLike) -> npt.ArrayLike:
    '''
    This function is a quick way to add the remaining eigenvectors that do not correspond to an eigenvalue to either U or V.)

    Args:
        A (np matrix): matrix that is to have eigenvectors added to it
        e_c (np array): array of associated eigenvalues, with any "used" eigenvalues replaced with 0
        evs (np matrx): matrix of associated eigenvectors
    Returns:
        A (np matrix): matrix A with new eigenvectors added to it
    '''
    #For every nonzero value in e_c, add that vector to A
    try:
        A_rows = len(A)
        for i in range(len(e_c)):
            v = e_c[i]
            if v != -1:
                A[A_rows-(i)] = np.ndarray.flatten(evs[:,i:i+1])
    except IndexError:
        #May occur if matrices are the wrong size
        pass
    return A

def L2_norm(A: npt.ArrayLike) -> float:
    '''
    This function finds the L2-norm of a matrix.

    Args:
        A (np matrix): matrix for which the L2-norm is to be found
    Returns:
        norm (float): L2-norm of A
    '''
    #Calculate the L2-norm of A
    val = 0
    for i in A:
        for j in i:
            val += float(j)**2
    return(np.sqrt(val))

def svd(A: npt.ArrayLike) -> tuple:
    '''
    This function finds the SVD decomposition of a matrix A, its condition number, and its inverse.

    Args:
        A (np matrix): matrix for which SVD decomposition, condition number, and inverse are to be found
    Returns:
        result (tuple): result[0] is itself a tuple, containing matrices U, S, and V of the SVD decomposition for A, for which U*S*V_T = A; 
            result[1] is the condition number of A, computed using A's infinity-norm;
            result[2] is the inverse of A, calculated by V_T*S_inv*U
    '''
    #Identify dimensions m, n
    try:
        m = np.shape(A)[0]
        n = np.shape(A)[1]
    except IndexError:
        raise Exception("Input is not a matrix")

    #Find eigenvalues of AAT, throw an error if l=0
    AAT = np.matmul(A, np.transpose(A))
    ATA = np.matmul(np.transpose(A), A)

    eigvals = np.linalg.eigvals(AAT)
    has_inverse = True
    try:
        for i in eigvals:
            if np.isclose(i,0):
                has_inverse = False
    except TypeError:
        raise Exception("Matrix has no real eigenvalues")
    
    #Create diagonal matrix (or vector)
    svs = np.sqrt(eigvals)

    #Find eigenvectors of AAT and ATA
    e_aat, ev_aat = np.linalg.eig(AAT)
    e_ata, ev_ata = np.linalg.eig(ATA)
    # print(np.linalg.eig(AAT))
    # print(np.linalg.eig(ATA))

    #Be sure the eigenvectors and values match when assigning U and V
    e_aat_c = e_aat.copy()
    e_ata_c = e_ata.copy()
    UT = np.zeros((m,m), float)
    VT = np.zeros((n,n), float)
    #Defined as such, as it is easier to build the tranpose matrix

    for i in range(len(eigvals)):
        #Identify what index it occurs at for both AAT and ATA, and assign the first col of U/V appropriately
        #No repeats allowed, so replace value with -1 when done
        #This is fine because no eigenvalues are allowed to be negative anyway
        #Copies have been made so that original eigenvalues are not modified
        v = eigvals[i]
        aat_i = int(np.where(np.isclose(e_aat_c, v))[0])
        ata_i = int(np.where(np.isclose(e_ata_c, v))[0])
        UT[i] = (np.ndarray.flatten(ev_aat[:,aat_i:aat_i+1]))
        VT[i] = (np.ndarray.flatten(ev_ata[:,ata_i:ata_i+1]))
        e_aat_c[aat_i] = -1
        e_ata_c[ata_i] = -1

    #We must account for the other eigenvectors of the bigger matrix
    if m > n:
        UT = add_evs(UT, e_aat_c, ev_aat)
    elif n > m:
        VT = add_evs(VT, e_ata_c, ev_ata)
    else:
        pass

    #Build S and SI from svs
    S = np.zeros((m,n), float)
    SI = np.zeros((n,m), float)
    for i in range(len(svs)):
        S[i][i] = svs[i]
        SI[i][i] = 1/svs[i]

    #Return the decomp USVT
    U = np.transpose(UT)
    V = np.transpose(VT)
    # print(U)
    # print(S)
    # print(V)
    SVD_decomp = (U,S,V)

    # print(np.matmul((np.matmul(U, S)),VT))

    #Return the inverse VS-1UT
    if has_inverse:
        AI = np.matmul((np.matmul(V, SI)),UT)
    else:
        raise InverseException(Exception)
    # print(np.matmul((np.matmul(V, SI)),UT))

    #Calculate condition number by calculating inf-norm of A and A-1
    if has_inverse:
        condition = float(L2_norm(A) * L2_norm(AI))
    else:
        raise Exception("Matrix has no finite condition numbers")
    # print(condition)

    return(SVD_decomp,condition,AI)
