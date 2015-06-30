import FourierX as fx
import string
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

Nout = 10001
tmax = 1.e5
x = fx.VarfromFile(10001,tmax,6)
times = np.linspace(0,tmax,Nout)

Npts = 3000
logPmin = np.log10(10.)
logPmax = np.log10(1.e5)
Ps = np.logspace(logPmin,logPmax,Npts)
ws = np.asarray([2*np.pi/P for P in Ps])

periodogram = signal.lombscargle(times,x[0][3],ws)

fig = plt.figure(figsize=(12,5))
ax = plt.subplot(111)
ax.plot(Ps,np.sqrt(4*periodogram/Nout))
#ax.set_xscale('log')
#ax.set_xlim([10**logPmin,10**logPmax])
ax.set_xlim([1,100])
ax.set_xlabel("Period (yrs)", fontsize=20)
ax.set_ylabel("Power", fontsize=20)
ax.tick_params(labelsize=20)

plt.show()
