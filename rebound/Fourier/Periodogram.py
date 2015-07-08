import FourierCoord as fc
import string
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

Nout = 1001
tmax = 5000.
NumPlanet = 6
times,[e,O,o] = fc.VarfromFile(Nout,NumPlanet,FileName='data/dt5_eOo.txt')

ecospulmega = e*np.cos(O+o)

Npts = 3000
logPmin = np.log10(10.)
logPmax = np.log10(1.e5)
Ps = np.logspace(logPmin,logPmax,Npts)
ws = np.asarray([2*np.pi/P for P in Ps])

periodogram = np.zeros([6,3000])

for i in range(NumPlanet):
    periodogram[i] = signal.lombscargle(times,ecospulmega[i],ws)

    fig = plt.figure(figsize=(12,5))
    ax = plt.subplot(111)
    ax.plot(Ps,np.sqrt(4*periodogram[i]/Nout))
    ax.set_xlim([1,100])
    ax.set_ylim([0,0.005])
    ax.set_xlabel("Period (yrs)", fontsize=20)
    ax.set_ylabel("Power", fontsize=20)
    ax.tick_params(labelsize=20)

plt.show()
