import scipy.integrate
import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

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
def A(i,j,diag):
    aa = alpha(Planets[i,1],Planets[j,1])
    aab = alphabar(Planets[i,1],Planets[j,1])
    if diag == True:
        return Planets[i,2]/4*(Planets[j,0]/(1+Planets[i,0])*aa*aab*b(1.,aa))
    if diag == False:
        return -Planets[i,2]/4*(Planets[j,0]/(1+Planets[i,0])*aa*aab*b(2.,aa))


#------------------------------------------------------------------------------
# Initialize variables for the planets
#------------------------------------------------------------------------------

# Planet = [ Mp/Ms, a(AU), n(deg/y), ecosw, esinw,I(deg), W(deg) ]

PlanetNames = ['K11b', 'K11c', 'K11d', 'K11e', 'K11f', 'K11g']
Planets = np.array([[ 1.3274e-5, 0.091089, 13140, 0.010, -0.011, 88.40, 0. ],
           [ 2.9076e-5, 0.106522, 10086, 0.005, -0.004, 91.17, 3.2 ],
           [ 2.8113e-5, 0.154241, 5788.5, -0.013, -0.009, 89.18, -33. ],
           [ 3.3817e-5, 0.193937, 4105.6, 0.020, -0.016, 88.743, -31. ],
           [ 1.1278e-5, 0.249489, 2813.7, 0.006, -0.017, 89.30, -32. ],
           [ 5.6889e-5, 0.463918, 1109.7, -0.26, 0.008, 90.23, -65. ]])


NumPlanets = len(Planets)
AMat = np.zeros([NumPlanets,NumPlanets])

for i in range(NumPlanets):
    for j in range(NumPlanets):
        if i == j:
            AMat[i,i] = sum(A(i,j,True) for j in range(NumPlanets) if not j==i)
        else:
            AMat[i,j] = A(i,j,False)

Eig = linalg.eig(AMat)
EigFreq = Eig[0]
EigVect = Eig[1]


SsinB = linalg.solve(EigVect,Planets[:,3])
ScosB = linalg.solve(EigVect,Planets[:,4])

B = np.arctan2(SsinB,ScosB)
S = ScosB/np.cos(B)

eig = S*EigVect

time = np.linspace(-100000.,100000.,10000)
h_i = np.array([[sum([eig[i,j]*np.sin(EigFreq[j]*t*np.pi/180.+B[j]) for j in range(NumPlanets)])  for i in range(NumPlanets)] for t in time])
k_i = np.array([[sum([eig[i,j]*np.cos(EigFreq[j]*t*np.pi/180.+B[j]) for j in range(NumPlanets)])  for i in range(NumPlanets)] for t in time])



ecc_j = (h_i**2+k_i**2)**(1./2.)

TotEcc = sum(np.transpose(ecc_j))

# These are important for getting labels in latex
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif','serif':['Helvetica']})

plt.figure()
plt.xlabel('Time (y)')
plt.ylabel('Eccentricity')
for i in range(NumPlanets):
    plt.plot(time, ecc_j[:,i], label = PlanetNames[i])
plt.legend()
plt.savefig('plots/6Planets_IndivEcc2.pdf')
'''
plt.figure()
plt.xlabel('Time (y)')
plt.ylabel('Eccentricity')
for i in range(NumPlanets-1):
    plt.plot(time, ecc_j[:,i], label = PlanetNames[i])
plt.legend()
plt.savefig('plots/5Planets_IndivEcc.pdf')
'''
plt.figure()
plt.xlabel('Time (y)')
plt.ylabel('Total Eccentricity of the System')
plt.plot(time, TotEcc)
#plt.savefig('plots/TotalEcc_Kepler11.pdf')

plt.show()
