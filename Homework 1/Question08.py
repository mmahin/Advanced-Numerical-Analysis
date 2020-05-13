import numpy as np
from numpy import array
import scipy
from scipy.linalg import svd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle
import seaborn as sns
import math


def my_SVD(X):
    # SVD
    U, s, VT = scipy.linalg.svd(X)
    print("U= \n",U)
    print("V Transpose= \n", VT)
    m=len(X)
    n=len(X[0])
    sigma = np.zeros((m, n))
    for i in range(min(m, n)):
        sigma[i, i] = s[i]
    print(np.diagonal(sigma))
    Usigma=np.dot(U,sigma)
    print("Usigma = \n", Usigma)
    V=np.transpose(VT)
    v1=V[:,0]
    v2=V[:,1]
    print("\nRight Singular v1 = ")
    print(v1)
    print("\nRight Singular v2 = ")
    print(v2)

    Usigma1=Usigma[:,0]
    Usigma2=Usigma[:,1]
    print("\nLeft Singular u1 = ")
    print(Usigma1)
    print("\nLeft Singular u2 = ")
    print(Usigma2)

    x = np.linspace(-1, 1, 100000)
    y = np.sqrt(1 - (x ** 2))
    #plt.plot(x, y, sns.color_palette().as_hex()[0])
    #plt.plot(x, -y, sns.color_palette().as_hex()[0])
    fig, ax = plt.subplots(1,2)


    circle1 = plt.Circle((0, 0), 1)
    ax[0].add_artist(circle1)
    ax[0].arrow(0, 0, v1[0], v1[1], head_width=0.05, head_length=0.01, fc='k', ec='k')
    ax[0].arrow(0, 0, v2[0], v2[1], head_width=0.05, head_length=0.01, fc='k', ec='k')
    leng=np.linalg.norm(v1)
    ax[0].set_xlim(-1*leng, leng)
    ax[0].set_ylim(-1*leng, leng)
    ax[0].text(v1[0], v1[1], r'$\vec{v1}$', color='black', size=18)
    ax[0].text(v2[0], v2[1], r'$\vec{v2}$', color='black', size=18)
    import matplotlib as mpl
    mean = [0, 0]
    leng = np.linalg.norm(Usigma1)
    if np.linalg.norm(Usigma2)< np.linalg.norm(Usigma1):
        width = np.linalg.norm(Usigma2)*2
        height = np.linalg.norm(Usigma1)*2
    else:
        width = np.linalg.norm(Usigma1) * 2
        height = np.linalg.norm(Usigma2) * 2
        leng = np.linalg.norm(Usigma2)

    def dotproduct(v1, v2):
        return sum((a * b) for a, b in zip(v1, v2))

    def length(v):
        return math.sqrt(dotproduct(v, v))

    def angle(v1, v2):
        return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

    ang= angle(Usigma1,[1,0])
    ang=np.degrees(ang)
    ell = mpl.patches.Ellipse(xy=mean, width=width, height=height, angle=90+ang)
    ax[1].arrow(0, 0, Usigma1[0], Usigma1[1], head_width=0.05, head_length=0.01, fc='k', ec='k')
    ax[1].arrow(0, 0, Usigma2[0], Usigma2[1], head_width=0.05, head_length=0.01, fc='k', ec='k')
    ax[1].set_xlim(-1*leng, leng)
    ax[1].set_ylim(-1*leng, leng)

    ax[1].add_patch(ell)
    ax[1].text(Usigma1[0], Usigma1[1], r'$\vec{sigma.U1}$', color='black', size=18)
    ax[1].text(Usigma2[0], Usigma2[1], r'$\vec{sigma.U2}$' ,color='black', size=18)

    # plt.show()
    plt.show()

A = scipy.array([[1, 2], [0, 2]])
B = scipy.array([[3, 0], [0, -2]])
C = scipy.array([[1, 1], [0, 0]])
print("\nSVD of A: \n",A)
Result_A=my_SVD(A)
print("\nSVD of B: \n",B)
Result_A=my_SVD(B)
print("\nSVD of C: \n",C)
Result_A=my_SVD(C)



