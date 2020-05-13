import numpy as np
import scipy
import scipy.linalg

def mgs(A):

    V = np.zeros(A.shape, dtype='float64')
    R = np.zeros(A.shape, dtype='float64')
    Q = np.zeros(A.shape, dtype='float64')
    dim = A.shape[1]
    for i in range(dim):
        V[:,i] = A[:,i]
    for i in range(0, dim):
        R[i,i] = scipy.linalg.norm(V[:,i])
        Q[:, i] = V[:,i]/R[i,i]
        for j in range(i+1, dim):
            R[i,j] = np.dot(np.transpose(Q[:,i]), V[:,j])
            V[:,j] = V[:,j] - np.dot(R[i,j],Q[:,i])
    return Q,R


A = scipy.array([[1,3,2], [4,5,6], [9,8,7]])
b = scipy.array([[13],[32],[46]])
print ("Matrix A =")
print(A)
print ("\nVector B =")
print(b)
Q,R=mgs(A)
print("\nQ Matrix=\n",Q)
print("\nR Matrix=\n",R)

print("\nDot of QR=\n",np.dot(Q,R))
bT = np.dot(np.transpose(Q), b)
X = scipy.linalg.solve(R, bT)

print("\nSolution of X=\n",X)

