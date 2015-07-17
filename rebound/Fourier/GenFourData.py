import rebound
import reboundxf

import numpy as np
from numpy import ndarray 
from numpy import cos, sin

import integration.integrator as itg
import integration.cfg as cfg
import integration.destabilization as dtb
import FourierCoord as fc

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.fftpack import fft

import sys
import time


def RotateSystem():

    # Find the angles corresponding to the angular momentum vector
    phi, theta = CalcRotAngles()

    # For the planets and the star rotate the system
    for i in range(cfg.NumPlanet+1):
        
        # Find the current position and velocity vector for
        # the planet or star
        x = np.array([[cfg.ps[i].x],[cfg.ps[i].y],[cfg.ps[i].z]])
        vx = np.array([[cfg.ps[i].vx],[cfg.ps[i].vy],[cfg.ps[i].vz]])

        # Calculate the matrix to rotate around the Y-axis by theta
        RotY = np.matrix([[cos(theta),0,-sin(theta)],
                          [0,1,0],
                          [sin(theta),0,cos(theta)]])

        # Calculate the matrix to roatte around Z-axis by phi
        RotZ = np.matrix([[cos(phi),sin(phi),0],
                         [-sin(phi),cos(phi),0],
                         [0,0,1]])

        # Rotate the position and velocity to new coordinates
        rotated = np.dot(RotY,np.dot(RotZ,x)) 
        rotatedvel = np.dot(RotY,np.dot(RotZ,vx))  

        # Update the particles in rebound
        cfg.ps[i].x = rotated[0]
        cfg.ps[i].y = rotated[1]
        cfg.ps[i].z = rotated[2]

        cfg.ps[i].vx = rotatedvel[0]
        cfg.ps[i].vy = rotatedvel[1]
        cfg.ps[i].vz = rotatedvel[2]



def CalcRotAngles():

    # Calculate angular momentum vector
    L = CalcAngMomVect()

    # Calculate the angles by which to rotate the system 
    phi = np.arctan2(L[1],L[0])
    theta = np.arccos(L[2])

    return phi, theta


def CalcAngMomVect():
    
    # Calculate the angular momentum of the star using
    # L = m(r x v)
    L = cfg.ps[0].m*np.cross(np.array([cfg.ps[0].x,cfg.ps[0].y,cfg.ps[0].z]),
                             np.array([cfg.ps[0].vx,cfg.ps[0].vy,cfg.ps[0].vz]))
    
    # Read the angular momentum of each planet from rebound
    # and add it to the total angular momentum
    for i in range(cfg.NumPlanet): 
        L = L+np.array(cfg.ps[i+1].calculate_orbit().hvec)*cfg.ps[i+1].m

    # Normalize the angular momentum vector to a unit vector
    L = L/np.sqrt(L[0]**2+L[1]**2+L[2]**2)
    
    return L

def FourierTransform(array,tmax,Nout, Plot=None):
    
    # fourier transform the array
    ft = fft(array)

    # Calculate the spacing between timesteps 
    dt = tmax/(Nout-1.)
    
    # Use time step to calculate frequencies
    # Max frequency is the Nyquist frequency
    freq = np.linspace(0.,1./(2.*dt),Nout/2)

    # If you want the frequency plotted
    if Plot == 'freq':

        # Plot all fourier transforms on one plot
        fc.PlotAllFreq(freq,ft)

    # Plot the fourier transform vs. period
    elif Plot == 'period':

        plt.figure()
        plt.xlabel('Period (yr)')
        plt.xlim([0,100])
        for i in range(cfg.NumPlanet):
            # Don't include first point (0)
            plt.plot(1./freq[1:],np.absolute(ft[i][1:Nout/2]),
                     label=cfg.Names[i])
        plt.legend()
        plt.show()

if __name__=="__main__":
    
    model='MigaIVa'
    integrator='ias15'

    VarOut = [['Omega','e','omega','a'],['x','y','z']]

    tmax = 50. #years
    Nout = 1001

    time1 = time.time()

    # Set up and perform the simulation in new coordinate frame
    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    rebound.add( m=1e-3, a=8.2, e=0.05, inc=89.*np.pi/180. )  #Jup1
    
    itg.ArrangeSys()

    cfg.NumPlanet = rebound.N-1

    RotateSystem()

    itg.ArrangeSys()
    dtb.InitOrbElem2(VarOut,Nout,tmax)
    
    dtb.Destab(tmax,VarOut,Nout)

    print 'Computational Time:', time.time()-time1

    # Calculate h and k
    ecospomega = np.array(cfg.e)*cos(np.array(cfg.Omega)+np.array(cfg.omega))
    esinpomega = np.array(cfg.e)*sin(np.array(cfg.Omega)+np.array(cfg.omega))

    # Print variables used to calculate h and k to file
    FileName = 'dt'+str(tmax/(Nout-1.))+'_Nout'+str(Nout)+'_eOo.txt'
    times = np.reshape(cfg.times,[1,Nout])
    eOo = np.concatenate((times,cfg.e,cfg.Omega,cfg.omega))
    itg.Print2File(eOo,FileName)

    # Plot h and k vs. time
    itg.PlotvTime(ecospomega,'$e\cdot cos(\~{\omega})$')
    itg.PlotvTime(esinpomega,'$e\cdot sin(\~{\omega})$')

    # Plot orbit and semi-major axis as reassurances
    itg.Plot_Orbit()
    itg.Plot_a()

    # Fourier transform h and plot fourier transform vs. period
    FourierTransform(ecospomega,tmax,Nout,Plot='period')

    plt.show()
