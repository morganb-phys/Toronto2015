import numpy as np

import scipy
from scipy.stats import norm
from scipy.fftpack import fft, fftfreq, fftshift
from scipy import signal
from scipy import optimize

import matplotlib.pyplot as plt

import math
import string
import sys

import integration.integrator as itg



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
        
        # Find largest peak in masked time series
        currentMax=np.amax(np.abs(fourier)*mask)

        #Make sure we're not just finding peaks in the noise
        if currentMax>threshold:

            # Find location of current max
            currentMaxLoc=np.argmax(np.abs(fourier)*mask)

            pulseDict[currentMaxLoc]=currentMax

            # Mask pulse from future searching
            width = PeakWidth(fourier,freq,currentMax,
                              currentMaxLoc,freqSpace) 

            # Zero out that location so it won't find same pulse
            mask[width[0]:width[1]] = 0.  
 
        
        else:
            # Plots the fourier transform for a single planet and the 
            # masked array to make sure all have been removed
            plt.figure()
            plt.plot(freq,np.abs(fourier)*mask)
            plt.plot(freq,np.abs(fourier))
            plt.show()
            break

    # Sorts the peaks from largest to smallest
    pulseList=sorted(pulseDict.keys(),key=lambda x: -pulseDict[x])

    # return the number of the pulse, height of pulse and 
    # frequency location of pulse
    return [(i,pulseDict[i],freq[i]) for i in pulseList]    

#Fitting functions to find width of peak
def f(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2.*c**2))

def laplace(x,a,b,c):
    return a*np.exp(-np.abs(x-b)/c)

# Find peak width
def PeakWidth(fourier,freq,currentMax,currentMaxLoc,freqSpace,method=None):

    # General rule is that the width of the peak is
    # height/5e2
    if method == 'FracHeight' or method == None:
        width = currentMax/5.e2
        
        # returns an array of the approximate location of peak
        return [int(max(0,math.floor(currentMaxLoc-width/freqSpace))),
                    int(math.ceil(currentMaxLoc+width/freqSpace))]

    #fits to the peak and uses std dev to define width
    elif method == 'fit':
 
        bins = PeakWidth(fourier,freq,currentMax,currentMaxLoc,freqSpace)
       
        popt,pcov = optimize.curve_fit(laplace,freq[bins[0]:bins[1]],
                                       np.abs(fourier[bins[0]:bins[1]]))
        freqmin = popt[1]-3.*popt[2]
        freqmax = popt[1]+3.*popt[2]
        
        # returns an array of the approximate location of peak
        return [max(0,math.floor(currentMaxLoc-freqmin/freqSpace)),
                    math.ceil(currentMaxLoc+freqmax/freqSpace)]
        

def VarfromFile(Nout,NumP,tmax = None,FileName = None):
    
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
            return x,time
                            
        else:
            print 'Unable to read the data file.'
            sys.exit()
            


if __name__=="__main__":

    Nout = 1001
    tmax=5000.
    
    x = VarfromFile(Nout,6,tmax=tmax,FileName='data/dt5_eOo.txt')[0]
    time = np.linspace(0,tmax,Nout)
    dt = tmax/(Nout-1.)
    Planets = ['k11b','k11c','k11d','k11e','k11f','k11g']
    Resonances = {}

    ecosPomega = x[0]*np.cos(x[1]+x[2])


    ft= fft(ecosPomega)
    freq = np.linspace(0.,1./(2.*dt), Nout/2)

    '''
    for i,planet in enumerate(Planets):
        print i
        Resonances[planet] = FindPeak(ft[i][0:Nout/2],freq,dt)
        
    
    for planet in Planets:
        print planet
        for frequency in Resonances[planet]:
            print frequency[2]
    '''
     
    itg.PlotLatex()

    plt.figure()
    plt.xlabel('Frequency $(yr^{-1})$')
    plt.xlim([0.01,0.1])
    plt.ylim([0,0.6])
    for i,planet in enumerate(Planets):
        plt.plot(freq,np.absolute(ft[i][0:Nout/2]),label=planet)
    plt.legend()
    plt.show()
