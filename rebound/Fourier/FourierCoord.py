import numpy as np
import scipy
import matplotlib.pyplot as plt
import string
from scipy.fftpack import fft, fftfreq, fftshift
from scipy import signal
import integration.integrator as itg
import sys
import scipy.signal
from scipy.stats import norm
import math


def FindPeak(fourier,freq,dt,span=None):
    # Based on code to find Giant Pulses in Matthew Quenneville's
    # pulsar-analysis repository at https://github.com/MatthewQuenneville
    

    locMin = np.argmin(np.absolute(fourier))
    
    threshold = np.mean(np.abs(fourier))
    #threshold=np.mean(np.abs(fourier[locMin-100:locMin+100])*40.)
    
    #print 'Threshold:', threshold
    freqSpace = freq[1]-freq[0]
    
    if span==None:
        span = [freq.min(),freq.max()]
    
    # Define a 'mask' such that previously found pulses are ignored
    mask=np.ones(len(freq))
    pulseDict={}

    i = 0
    while True:
        i+=1
        #print i
        # Find largest pulse in masked time series
        currentMax=np.amax(np.abs(fourier)*mask)

        # If larger than lowerbound, count as new pulse. Otherwise,
        # break and notify user that not enough giant pulses exist
        if currentMax>threshold:

            # Find location of current max
            currentMaxLoc=np.argmax(np.abs(fourier)*mask)
            pulseDict[currentMaxLoc]=currentMax

            
            # Mask pulse from future searching
            width = currentMax/5.e4
            #print np.arange(max(0,math.floor(currentMaxLoc-width/freqSpace)),
                           # math.ceil(currentMaxLoc+width/freqSpace),1)
            mask[max(0,math.floor(currentMaxLoc-width/freqSpace)):
                     math.ceil(currentMaxLoc+width/freqSpace)] = 0.     
        
        else:
            plt.figure()
            plt.plot(freq,np.abs(fourier)*mask)
            plt.plot(freq,np.abs(fourier))
            plt.show()
            break

    pulseList=sorted(pulseDict.keys(),key=lambda x: -pulseDict[x])

    return [(i,pulseDict[i],freq[i]) for i in pulseList]    

    

def PeakWidth():
    pass
    
def FindFreq():
    pass

def VarfromFile(Nout,NumP,tmax = None,FileName = None):
    # Read the data from the x-axis
    
    if FileName == None:
        FileName = 'data/'+"{:.0E}".format(tmax)+'Yanqin_X.txt'

    File = open(FileName,'r')
    NumVar = len(string.split(File.readline()))/NumP

    x = np.zeros([NumVar,NumP,Nout])

    if len(string.split(File.readline()))%NumP == 0:

        for i,line in list(enumerate(File)):
            Line = string.split(line)
            for k in range(NumVar):
                Vars=Line[k*NumP:(k+1)*NumP]
                for j,word in list(enumerate(Line)):
                    x[k][j][i] = float(word)

        File.close()
        return x


    else:
        IncTime = raw_input('Did your file include a time index (y/n)?\t')

        if IncTime=='y':

            time = np.zeros(Nout)
                
            for i,line in list(enumerate(File)):
                    
                Line = string.split(line)
                if len(Line) != 0: 

                    time[i] = Line[0]
                
                    for k in range(int(NumVar)):
                        Vars=Line[k*NumP+1:(k+1)*NumP+1]
                        
                        for j,word in list(enumerate(Vars)):       
                            x[k][j][i] = float(word)
     
            File.close()
            return time,x
                            
        else:
            print 'Unable to read the data file.'
            sys.exit()
            


if __name__=="__main__":

    Nout = 10001
    tmax=1.e2
    x = VarfromFile(Nout,6,tmax=tmax)[0]
    time = np.linspace(0,tmax,Nout)
    dt = tmax/(Nout-1.)
    Planets = ['k11b','k11c','k11d','k11e','k11f','k11g']
    Resonances = {}


    ft= fft(x)
    freq = np.linspace(0.,1./(2.*dt), Nout/2)


    for i,planet in enumerate(Planets):
        print i
        Resonances[planet] = FindPeak(ft[i][0:Nout/2],freq,dt)
        
    
    for planet in Planets:
        print planet
        for frequency in Resonances[planet]:
            print frequency[2]
         
    plt.figure()
    plt.xlabel('Frequency $(yr^{-1})$')
    for i,planet in enumerate(Planets):
        plt.plot(freq,np.absolute(ft[i][0:Nout/2]),label=planet)
    plt.legend()
    plt.show()
    '''
    plt.figure()
    plt.ylabel('Fourier Amplitude')
    plt.xlabel('Frequency $(^\circ yr^{-1})$')
    for i in range(6):
        plt.plot(180./np.pi*freq, np.abs(ft[i][0:Nout/2]))
        
    plt.figure()
    plt.ylabel('Fourier Amplitude')
    plt.xlabel('Period (yr)')
    for i in range(6):
        plt.semilogx(1./freq, np.abs(ft[i][0:Nout/2]))

    plt.figure()
    plt.ylabel('Fourier Amplitude')
    plt.xlabel('Period (yr)')
    plt.xlim([20,2000])
    for i in range(6):
        plt.plot(1./freq, np.abs(ft[i][0:Nout/2]))

    plt.figure()
    plt.xlim([5*10**2,10**5])
    plt.ylim([0,30])
    plt.ylabel('Fourier Amplitude')
    plt.xlabel('Period (yr)')
    for i in range(6):
        plt.semilogx(1./freq, np.abs(ft[i][0:Nout/2]))

    plt.show()


'''
