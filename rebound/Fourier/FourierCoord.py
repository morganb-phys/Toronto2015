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
import integration.cfg as cfg


#-----------------------------------------------------------------------------
# FUNCTIONS TO FIND PEAKS IN THE FOURIER TRANSFORM (FREQUENCIES OF THE SYSTEM)
#-----------------------------------------------------------------------------

def FindPeak(fourier,freq,dt,span=None,Plot=False):
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

    while True:

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
            if Plot == True:
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


# Find peak width
def PeakWidth(fourier,freq,currentMax,currentMaxLoc,freqSpace,method=None):

    # General rule is that the width of the peak is
    # height/5e2
    if method == 'FracHeight' or method == None:
        width = currentMax/5.e2
        
        # returns an array of the approximate location of peak
        return [int(max(0,math.floor(currentMaxLoc-width/freqSpace))),
                    int(math.ceil(currentMaxLoc+width/freqSpace))]

    #fits to the peak and uses std dev to define width - NOT WORKING
    elif method == 'fit':
 
        bins = PeakWidth(fourier,freq,currentMax,currentMaxLoc,freqSpace)
       
        popt,pcov = optimize.curve_fit(laplace,freq[bins[0]:bins[1]],
                                       np.abs(fourier[bins[0]:bins[1]]))
        freqmin = popt[1]-3.*popt[2]
        freqmax = popt[1]+3.*popt[2]
        
        # returns an array of the approximate location of peak
        return [max(0,math.floor(currentMaxLoc-freqmin/freqSpace)),
                    math.ceil(currentMaxLoc+freqmax/freqSpace)]


#Fitting functions to find width of peak
def f(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2.*c**2))

def laplace(x,a,b,c):
    return a*np.exp(-np.abs(x-b)/c)



#-----------------------------------------------------------------------------
# USEFUL FOURIER TRANSFORM OUTPUTS
#-----------------------------------------------------------------------------

# Print frequencies found for each planet to screen
# Note for later (add in print to file)
def printFreq(Planets,PrintTo='screen',FName=None,Title=None):

    if PrintTo == 'screen':

        # Run through each planet
        for planet in Planets:
            
            # Print planet name
            print planet
            
            # Print frequencies of that planet
            for frequency in Resonances[planet]:
                print frequency[2]

    elif PrintTo == 'file':
        
        # If no file name given, use a generic one
        if FName == None:
            FName = 'FourierTransformData.txt'
            
        # Open file to write to it
        File = open('data/'+FName,'a')

        # Title for the data set - name of file read in from
        if Title != None:
            File.write('\n'+Title+'\n\n')

        for planet in Planets:

            File.write(planet+'\n')
            for frequency in Resonances[planet]:
                File.write(str(frequency[2])+'\n')
            File.write('\n')
            
        File.close()


def PlotAllFreq(freq,fourier,xspan=None,save=False,PlotSave=None):

    plt.figure()

    plt.xlabel('Frequency $(yr^{-1})$')

    # If there is a specific range over which to plot
    if xspan != None:

        # Set the x-axis range to values specified 
        plt.xlim(xspan)
        
        # Find the index number of the frequency closest to the min
        # and max that was specified
        xminInd = (abs(freq-xspan[0])).argmin()
        xmaxInd = (abs(freq-xspan[1])).argmin()

        # Find the maximum y value in that range and take 110% of it
        ymax = 1.1*(np.abs(fourier[:,xminInd:xmaxInd])).max()

        # Scale the y-axis to the ymax
        plt.ylim([0,ymax])
    
    # For each planet, plot the fourier transform
    for i,planet in enumerate(Planets):
        plt.plot(freq,np.absolute(ft[i][0:Nout/2]),
                 color = cfg.Colours[i],label=planet)

    # plot the legend
    plt.legend()

    if save == True:
        plt.savefig('plots/'+PlotSave+'.png')


#-----------------------------------------------------------------------------
# READ VARIABLES FROM FILE
#-----------------------------------------------------------------------------


# Read variables from file with format:
# x1 x2 x3 y1 y2 y3 z1 z2 z3 (Grouped by coordinate, not planet)
def VarfromFile(NumP,FileName,tmax = None):

    # Open the file
    File = open(FileName,'r')
    
    # find the number of variables
    # number of variabless in line/number of planets
    NumVar = len(string.split(File.readline()))/NumP
       
    # Find the number of outputs be summing the number of non-empty lines
    Nout = 0
    for line in File:
        if line.strip():
            Nout += 1
    
    # Places the read cursor at the beginning of the file
    File.seek(0)

    # Initialize variable to hold read values
    x = np.zeros([NumVar,NumP,Nout])

    # If number of values in a line is a multiple of the number of planets
    if len(string.split(File.readline()))%NumP == 0:

        # Read each line
        for i,line in list(enumerate(File)):
            
            # Split line into variables
            Line = string.split(line)

            # For each variable
            for k in range(NumVar):
                
                # Select the chunk of line that contains that variable
                Vars=Line[k*NumP:(k+1)*NumP]

                # Store it in x 
                for j,word in list(enumerate(Line)):
                    x[k][j][i] = float(word)

        File.close()
        return x, Nout

    # If the number of values is not a multpile of the number of planets, 
    # there is most likely a time column included
    else:

        # Check if a time column was included
        IncTime = raw_input('Did your file include a time index (y/n)?\t')
        
        # If it was, read in the variables
        if IncTime=='y':

            # initialize time variable
            time = np.zeros(Nout)
                
            # Run through each line
            for i,line in list(enumerate(File)):
                    
                # Split line into variables
                Line = string.split(line)
                
                # Sometimes there is a blank line at beginning 
                # or end of the file
                if line.strip(): 

                    # Set line variable
                    time[i] = Line[0]
                
                    # For each variable
                    for k in range(int(NumVar)):
                        
                        # Select the range of the line with that variable
                        Vars=Line[k*NumP+1:(k+1)*NumP+1]
                        
                        # Save the variables to x
                        for j,word in list(enumerate(Vars)):       
                            x[k][j][i] = float(word)
     
            File.close()
            return x,time,Nout
                            
        # If it wasn't, the file is not in a recognizeable format
        else:
            print 'Unable to read the data file.'
            sys.exit()
            


if __name__=="__main__":

    NumPlanet=6
    Planets = cfg.Names[:NumPlanet]
    varToFour = 'pomega'

    # Check to see if a file was included
    if len(sys.argv) == 1:
        
        # No file. Are they okay with the standard format
        print 'The standard format is: data/tmaxYanqin_X.txt'
        cont = raw_input('Do you want to continue with this format (y/n)?\t')
        
        # If they are, read in the file and fourier transform x
        if cont=='y':

            # Read in tmax to read the correct file
            tmax = float(raw_input('What was the maximum time?\t'))
            File = 'data/'+"{:.0E}".format(tmax)+'Yanqin_X.txt'
            varToFour = 'x'
        
        # They did not specify and don't want to use stadard format. Quit.
        else:
            print 'Please run the code as:'
            print '\tpython FourierCoord.py path/to/file.txt'
            sys.exit()

    # If they included a file, use that one
    else:
        File = sys.argv[1]
   
    # Read variables from file
    ReadIn = VarfromFile(NumPlanet,tmax=tmax,FileName=File)
    
    # If there was a time column
    if len(ReadIn) == 2:
        
        x = ReadIn[0]
        Nout = ReadIn[1]
        time = np.linspace(0,tmax,Nout)

    # If there wasn't a time column
    else:
        
        x = ReadIn[0]
        time = ReadIn[1]
        Nout = ReadIn[2]

    # Assumes evenly spaced time steps
    dt = time[1]-time[0]

    if varToFour == 'pomega':
        # Calculate variable to be fourier transformed
        ecosPomega = x[0]*np.cos(x[1]+x[2])

        # Fourier transform the variable and calculate frequencies
        ft= fft(ecosPomega)

    if varToFour == 'x':
        x = x[0]
        ft = fft(x)

    # Upper limit on frequency is the Nyquist Frequency
    freq = np.linspace(0.,1./(2.*dt), Nout/2)

    # Initialize dictionary
    Resonances = {}

    # Runs through each planet
    for i,planet in enumerate(Planets):

        # Index of dictionary is name of planet
        # Value is a list of frequency peak locations
        Resonances[planet] = FindPeak(ft[i][0:Nout/2],freq,dt)

    # Prints the frequencies for each planet
    printFreq(Planets,PrintTo='file',Title=File)
     
    # Plot nicely in latex
    itg.PlotLatex()

    # Set range for x-axis in plot
    xspan = [0.01,0.1]

    # If a span is specified, include it in your plots save name
    if len(xspan)==2:
        PlotName = varToFour+'_'+str(xspan[0])+'to'+str(xspan[1])
    else:
        PlotName = varToFour+'_all'

    # Plot all of the fourier transforms of each planet together
    PlotAllFreq(freq,ft,xspan=[0.01,0.1],save=True,PlotSave=PlotName)

    plt.show()
