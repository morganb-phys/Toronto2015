import rebound
import reboundxf

import numpy as np
from numpy import ndarray 
from numpy import cos, sin

import integration.integrator as itg
import integration.cfg as cfg
import integration.destabilization as dtb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.fftpack import fft

import sys
import time


def RotateSystem():

    phi, theta = CalcRotAngles()

    for i in range(cfg.NumPlanet+1):
        x = np.array([[cfg.ps[i].x],[cfg.ps[i].y],[cfg.ps[i].z]])
        vx = np.array([[cfg.ps[i].vx],[cfg.ps[i].vy],[cfg.ps[i].vz]])

        RotY = np.matrix([[cos(theta),0,-sin(theta)],
                          [0,1,0],
                          [sin(theta),0,cos(theta)]])
        RotZ = np.matrix([[cos(phi),sin(phi),0],
                         [-sin(phi),cos(phi),0],
                         [0,0,1]])
        rotated = np.dot(RotY,np.dot(RotZ,x)) 
        rotatedvel = np.dot(RotY,np.dot(RotZ,vx))  

        cfg.ps[i].x = rotated[0]
        cfg.ps[i].y = rotated[1]
        cfg.ps[i].z = rotated[2]

        cfg.ps[i].vx = rotatedvel[0]
        cfg.ps[i].vy = rotatedvel[1]
        cfg.ps[i].vz = rotatedvel[2]



def CalcRotAngles():
    L = CalcAngMomVect()

    phi = np.arctan2(L[1],L[0])
    theta = np.arccos(L[2])

    return phi, theta


def CalcAngMomVect():
    
    L = cfg.ps[0].m*np.cross(np.array([cfg.ps[0].x,cfg.ps[0].y,cfg.ps[0].z]),
                             np.array([cfg.ps[0].vx,cfg.ps[0].vy,cfg.ps[0].vz]))
    for i in range(cfg.NumPlanet): 
        L = L+np.array(cfg.ps[i+1].calculate_orbit().hvec)*cfg.ps[i+1].m
    L = L/np.sqrt(L[0]**2+L[1]**2+L[2]**2)
    
    return L

def FourierTransform(array,tmax,Nout, Plot=None):
    
    ft = fft(array)
    dt = tmax/(Nout-1.)
    freq = np.linspace(0.,1./(2.*dt),Nout/2)

    if Plot == 'freq':
        plt.figure()
        plt.xlabel('Frequency $(yr^{-1})$')
        for i in range(cfg.NumPlanet):
            plt.plot(freq,np.absolute(ft[i][0:Nout/2]),label=cfg.Names[i])
        plt.legend()
        plt.show()

    elif Plot == 'period':
        plt.figure()
        plt.xlabel('Period (yr)')
        plt.xlim([0,100])
        for i in range(cfg.NumPlanet):
            plt.plot(1./freq[1:],np.absolute(ft[i][1:Nout/2]),
                     label=cfg.Names[i])
        plt.legend()
        plt.show()

if __name__=="__main__":
    
    model='MigaIVa'
    integrator='ias15'

    VarOut = [['Omega','e','omega','a'],['x','y','z']]

    tmax = 5000. #years
    Nout = 1001

    time1 = time.time()

    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    #rebound.add( m=1e-3, a=8.2, e=0.05, inc=89.*np.pi/180. )  #Jup1
    
    itg.ArrangeSys()

    cfg.NumPlanet = rebound.N-1

    RotateSystem()

    itg.ArrangeSys()
    dtb.InitOrbElem2(VarOut,Nout,tmax)
    
    dtb.Destab(tmax,VarOut,Nout)

    print 'Computational Time:', time.time()-time1

    ecospulmega = np.array(cfg.e)*cos(np.array(cfg.Omega)+np.array(cfg.omega))
    esinpulmega = np.array(cfg.e)*sin(np.array(cfg.Omega)+np.array(cfg.omega))

    print np.array(cfg.times).shape
    FileName = 'dt'+str(tmax/(Nout-1.))+'_Nout'+str(Nout)+'_eOo.txt'
    times = np.reshape(cfg.times,[1,Nout])
    eOo = np.concatenate((times,cfg.e,cfg.Omega,cfg.omega))
    itg.Print2File(eOo,FileName)

    itg.PlotvTime(ecospulmega,'$e\cdot cos(\~{\omega})$')
    itg.PlotvTime(esinpulmega,'$e\cdot sin(\~{\omega})$')

    itg.Plot_Orbit()
    itg.Plot_a()
    
    plt.show()

    FourierTransform(ecospulmega,tmax,Nout,Plot='period')

    
