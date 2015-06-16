import scipy.integrate
import numpy as np
from numpy import linalg
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

# Planet = [ Mp/Ms, a(AU), n(deg/y), ecosw, esinw,I(deg), W(deg) ]

PlanetNames = ['K11b', 'K11c', 'K11d', 'K11e', 'K11f', 'K11g','Jup1','Jup2']
Planets = np.array([[[ 1.3274e-5, 0.091089, 13140, 0.010, -0.011, 88.40, 0. ],
                    [ 2.9076e-5, 0.106522, 10086, 0.005, -0.004, 91.17, 3.2 ],
                    [ 2.8113e-5, 0.154241, 5788.5, -0.013, -0.009, 89.18, -33. ],
                    [ 3.3817e-5, 0.193937, 4105.6, 0.020, -0.016, 88.743, -31. ],
                    [ 1.1278e-5, 0.249489, 2813.7, 0.006, -0.017, 89.30, -32. ],
                    [ 5.6889e-5, 0.463918, 1109.7, -0.26, 0.008, 90.23, -65. ],
                    [ 1.0000e-3, 5.000000, 123.12, 0.05, 0.000, 89.00, -32. ],
                    [ 1.0000e-3, 7.600000+i, 360./(7.6+i)**(3/2), 0.00, 0.035, 91.00, -33. ]] for i in np.linspace(0.01,0.3,10)])


pool = rebound.InterruptiblePool()
Freq = pool.map(EigFrequencies, Planets)

Freq = np.array(Freq)

print (Freq[0,:]-Freq[9,:]

plt.figure()
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,0],label='Freq1')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,1],label='Freq2')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,2],label='Freq3')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,3],label='Freq4')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,4],label='Freq5')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,5],label='Freq6')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,6],label='Freq7')
plt.plot(7.6+np.linspace(0.01,0.3,10),Freq[:,7],label='Freq8')
plt.legend()
plt.show()
