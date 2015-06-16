import scipy.integrate
import numpy as np
from numpy import linalg
from numpy import random
import matplotlib.pyplot as plt
import rebound

#------------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------------

def alpha(a1,a2):
    if a1<a2:
        return a1/a2
    else: 
        return a2/a1

def alphabar(a1,a2):
    if a1<a2:
        return a1/a2
    else: 
        return 1.

def b_integrand(x,j,aa):
    return np.cos(j*x)/(1-2.*aa*np.cos(x)+aa**2)**(3./2.)

def b(j,aa):
    I = scipy.integrate.quad(b_integrand,0.,2.*np.pi,args=(j,aa))
    return I[0]/np.pi

# look at summing it in the function
def A(Planets,i,j,diag):
    aa = alpha(Planets[i,1],Planets[j,1])
    aab = alphabar(Planets[i,1],Planets[j,1])
    if diag == True:
        return Planets[i,2]/4*(Planets[j,0]/(1+Planets[i,0])*aa*aab*b(1.,aa))
    if diag == False:
        return -Planets[i,2]/4*(Planets[j,0]/(1+Planets[i,0])*aa*aab*b(2.,aa))

def EigFrequencies(Planets):

    NumPlanets = len(Planets)
    AMat = np.zeros([NumPlanets,NumPlanets])

    for i in range(NumPlanets):
        for j in range(NumPlanets):
            if i == j:
                AMat[i,i] = sum(A(Planets,i,j,True) for j in range(NumPlanets) if not j==i)
            else:
                AMat[i,j] = A(Planets,i,j,False)

    Eig = linalg.eig(AMat)
    EigFreq = Eig[0]

    return EigFreq



#------------------------------------------------------------------------------
# Initialize variables for the planets
#------------------------------------------------------------------------------

# Planets : m_planet/m_star, a (AU), n (degrees/year)

K11 = np.array([[[1e-3,5.,5.76],
                 [1e-3,5.+i,360./(i**(3/2))]] for i in np.linspace(1.6,1.7,11)])

print np.linspace(1.6,1.7,11)

pool = rebound.InterruptiblePool()
Frequencies = pool.map(EigFrequencies,K11)

Frequencies = np.array(Frequencies)

print Frequencies
