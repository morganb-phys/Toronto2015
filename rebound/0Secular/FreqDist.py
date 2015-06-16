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


def pdf(x,xmin,xmax,mean):

    pi = np.pi
    norm = 6/np.sqrt(2*pi)/(1+(xmax-mean)**2/(mean-xmin))

    if x<mean:
        pdf = norm/(mean-xmin)*np.exp(-9./2.*(x-mean)**2/(mean-xmin)**2)
    else:
        pdf = norm/(mean-xmin)*np.exp(-9./2.*(x-mean)**2/(xmax-mean)**2)

    return pdf
 
def randomvariate(pdf,n,xmin,xmax,mean):  
 
    # Calculates the minimal and maximum values of the PDF in the desired  
    # interval. The rejection method needs these values in order to work  
    # properly.  
    x=np.linspace(xmin,xmax,1000)
    y=np.zeros(len(x))
    for i,xi in list(enumerate(x)):
        y[i]=pdf(xi,xmin,xmax,mean)  
    pmin=y.min()  
    pmax=y.max() 

    # Counters  
    naccept=0  
    ntrial=0  

    # Keeps generating numbers until we achieve the desired n  
    ran=np.array([]) # output list of random numbers  
 
    while naccept<n:  
        x=np.random.uniform(xmin,xmax) # x'  
        y=np.random.uniform(pmin,pmax) # y'  
      
        if y<pdf(x,xmin,xmax,mean):  
            ran = np.append(ran,x)
            naccept=naccept+1  
        ntrial=ntrial+1  
    
    ran=np.asarray(ran)  
    
    return ran,ntrial  




#------------------------------------------------------------------------------
# Initialize variables for the planets
#------------------------------------------------------------------------------

mass = np.array([ 4.2,9.2,8.9,10.7,3.6,18. ])/320071.
emass = np.array([[3.0,2.4],[6.9,3.8],[2.7,3.5],[2.1,2.4],[2.6,5.4],[15.,24.]])/320071.

a = [ 0.091089,0.106522,0.154241,0.193937,0.249489,0.463918 ]
ea = [[1.1e-5,1.3e-5],[1.2e-5,7e-6],[1e-5,1.9e-5],[2.1e-5,1.5e-5],[2.5e-5,3.9e-5],[2.8e-5,5.9e-5]]

n = [  12728.,10065.,5776.6,4097.2,2807.9,1107.4]
en = [[3.,2.],[2.,1.5],[1.1,0.6],[0.8,0.5],[0.7,0.4],[0.2,0.1]]

NumFreq = 10000

K11 = np.zeros([NumFreq,6,3])

for i in range(6):

    K11[:,i,0] = randomvariate(pdf,n=NumFreq,mean=mass[i],xmin=mass[i]-emass[i][0],xmax=mass[i]+emass[i][1])[0]
    K11[:,i,1] = randomvariate(pdf,n=NumFreq,mean=a[i],xmin=a[i]-ea[i][0],xmax=a[i]+ea[i][1])[0]
    K11[:,i,2] = randomvariate(pdf,n=NumFreq,mean=n[i],xmin=n[i]-en[i][0],xmax=n[i]+en[i][1])[0]

'''
K11 = np.zeros([NumFreq,6,3])
for k in range(NumFreq):
    for i in range(6):
        K11[k,i,0] = random.triangular(mass[i]-emass[i][0],mass[i],mass[i]+emass[i][1])
        K11[k,i,1] = random.triangular(a[i]-ea[i][0],a[i],a[i]+ea[i][1])
        K11[k,i,2] = random.triangular(n[i]-en[i][0],n[i],n[i]+en[i][1])
'''
'''
K11 = np.zeros([NumFreq,6,3])
for k in range(NumFreq):
    for i in range(6):

        K11[k,i,0] = random.uniform(mass[i]-emass[i][0],mass[i]+emass[i][1])
        K11[k,i,1] = random.uniform(a[i]-ea[i][0],a[i]+ea[i][1])
        K11[k,i,2] = random.uniform(n[i]-en[i][0],n[i]+en[i][1])
'''

pool = rebound.InterruptiblePool()
Frequencies = pool.map(EigFrequencies,K11)

Frequencies = np.array(Frequencies)

#files = open('data/Freq_triangular'+str(NumFreq)+'.txt','a')
#files = open('data/Freq_uniform'+str(NumFreq)+'.txt','a')
files = open('data/Freq_Gaussian'+str(NumFreq)+'.txt','a')

for freq in Frequencies:
    files.write(str(freq[0])+' '+str(freq[1])+' '+str(freq[2])+' '+str(freq[3])+' '+str(freq[4])+' '+str(freq[5])+'\n')
files.write('\n')
files.close


plt.figure()
plt.hist(Frequencies[:,0])

plt.figure()
plt.hist(Frequencies[:,1])

plt.figure()
plt.hist(Frequencies[:,2])

plt.figure()
plt.hist(Frequencies[:,3])

plt.figure()
plt.hist(Frequencies[:,4])

plt.figure()
plt.hist(Frequencies[:,5])

plt.show()

