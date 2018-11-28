#Projet MAP

import numpy as np
import scipy.integrate as integrate
import numpy.linalg as LA
from scipy import signal

#Parameters
lbda=10.
a=2.5
esp=0.2
eps_0=1.
eps_1=1.
J1=3.
J2=3.
s1=0.25
x0=0.5

a_0=0.28
a_1=0.063*2

#Utils
def inte(f):
    """Compute the integrale of f between -pi/2 and pi/2"""
    inte=integrate.quad(f,-np.pi/2,np.pi/2)
    return (1/np.pi)*inte[0]

def cos_ax(x):
    """return the function x -> cos(ax)"""
    return np.cos(a*x)

def sin_ax(x):
    """return the function x -> sin(ax)"""
    return np.sin(a*x)

def cos_ax2(x):
    """return the function x -> cos(ax)**2"""
    return np.cos(a*x)*np.cos(a*x)

def sin_ax2(x):
    """return the function x -> sin(ax)**2"""
    return np.sin(a*x)*np.sin(a*x)

def rank(M):
    """return the rank of the matrix M"""
    return LA.matrix_rank(M)

# Controllability matrix

def Kalman(lbda=10.,
           a=2.,
           esp=0.2,
           eps_0=-1.,
           eps_1=1.,
           esp_2=1.,
           J1=-7.3,
           J2=3.,
           s1=0.25,
           x0=0.5):
    K= np.array([[eps_0 , eps_0*np.sqrt(abs(J1))*inte(cos_ax) , 0],
            [eps_1*np.sqrt(abs(J1))*inte(cos_ax) , J1*inte(cos_ax2) , 0],
            [0 , 0 , esp_2*J2*inte(sin_ax2)]])

    A_lin = -1.*np.identity(n=3,dtype=np.float32) + lbda*s1*K

    B = esp * np.array([[1. , 0. , 0.],
                        [0. , 1/np.sqrt(abs(J1)) , 1/np.sqrt(abs(J1))],
                        [0 , 1/np.sqrt(abs(J2)) , -1/np.sqrt(abs(J2))]])

    B2=esp*np.array([[a_0,a_1*np.cos(a+x0),-a_1*np.cos(a*x0)]])

    C=np.concatenate((B,np.dot(A_lin,B),np.dot(A_lin**2,B)),axis=0)

    print(C)
    return rank(C)


print(Kalman())


#K=signal.place_poles(A, B_2D, P)

#print("Matrix K is :")
#   print(K.gain_matrix)


