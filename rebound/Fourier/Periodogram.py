import FourierCoord as fc
import string
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

NumPlanet = 6

# Choose the range you want to look for periods
PeriodRange = [10.,1.e3]

# Read variables from file
[e,O,o],times,Nout = fc.VarfromFile(NumPlanet,'data/dt5_eOo.txt')

# Calculate variable of which you want to find the periods 
ecospomega = e*np.cos(O+o)

# Initialize number of points in periodogram
Npts = Nout*3

# Find power of period range
logPmin = np.log10(PeriodRange[0])
logPmax = np.log10(PeriodRange[1])

# Create logscale of periods in PeriodRange
Ps = np.logspace(logPmin,logPmax,Npts)

# Find correspunding angular frequencies
ws = np.asarray([2*np.pi/P for P in Ps])

# Initialize array to hold periodograms
periodogram = np.zeros([NumPlanet,Npts])

# Calculate periodogram for each planet
for i in range(NumPlanet):

    # Note: power given by periodogram is not actually the power
    # however the relative size of peaks to each other is accurate
    periodogram[i] = signal.lombscargle(times,ecospomega[i],ws)

    # Plot periodogram for the ith planet
    fig = plt.figure(figsize=(12,5))
    ax = plt.subplot(111)
    ax.plot(Ps,np.sqrt(4*periodogram[i]/Nout))
    ax.set_xlim([1,100])
    ax.set_ylim([0,0.005])
    ax.set_xlabel("Period (yrs)", fontsize=20)
    ax.set_ylabel("Power", fontsize=20)
    ax.tick_params(labelsize=20)

plt.show()
