import numpy as np
import rebound
from numpy import ndarray 
import sys
import integrator as itg
import cfg as cfg
import Secular.migration as mig
import reboundxf
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from numpy import cos, sin

def InitOrbElem2(VarOut,Nout,tmax):

    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    #Number of Planets
    cfg.NumPlanet = rebound.N-1
    cfg.Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    cfg.Colours = ['r','y','g','b','c','k','m','#fa8a29']
    
    # Initialize Orbital Elements
    for OrbVar in OrbOut:
        setattr(cfg,OrbVar,[[] for x in range(cfg.NumPlanet)])

    for CoordVar in CoordOut:
        setattr(cfg,CoordVar,[[] for x in range(cfg.NumPlanet)])

    # The times at which it will output the orbital elements and coordinates
    cfg.times = []

def Destab(tmax,VarOut,Nout):

    dt = tmax/(Nout-1.)

    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    step = 0
    
    a0 = np.zeros(cfg.NumPlanet)
    for i in range(cfg.NumPlanet):
        a0[i] = cfg.ps[i+1].calculate_orbit().a
    a = a0
    aHigh = 1.2*a0
    aLow = 0.8*a0
    Time = 0
    
    print '\nTime ellapsed:'

    while CheckSMA(a,aHigh,aLow)==True and Time<tmax:
 
        if step%(Nout/10) == 1:
            print '\t'+str(Time)+' yr'

        Time = rebound.t
        
        cfg.times.append(Time)
        
        for j in range(cfg.NumPlanet):
            a[j] = cfg.ps[j+1].calculate_orbit().a
            for OrbVar in OrbOut:
                getattr(cfg,OrbVar)[j].append(getattr(cfg.ps[j+1].calculate_orbit(),OrbVar))
          
            for CoordVar in CoordOut:
                getattr(cfg,CoordVar)[j].append(getattr(cfg.ps[j+1],CoordVar))
                
        rebound.integrate(Time+dt)

        step += 1
        
    print Time

def CheckSMA(a,aHigh,aLow):

    for i in range(cfg.NumPlanet):
        if (a[i]<aLow[i] or a[i]>aHigh[i]) and i != len(a)-1:
            print cfg.Names[i], a[i]
            return False
        elif a[i]>aLow[i] and a[i]<aHigh[i] and i == len(a)-1:
            return True


if __name__=="__main__":
    
    model='MigaIVa'
    integrator='ias15'

    VarOut = [['Omega','e','omega','a'],['x','y','z']]

    tmax = 1000. #years
    Nout = 11
    dt = tmax/(Nout-1.)

    time1 = time.time()

    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    rebound.add( m=1e-3, a=8.2, e=0.05, inc=89.*np.pi/180. )  #Jup1
    
    itg.ArrangeSys()

    Destab(tmax,VarOut,Nout)

    print 'Computational Time:', time.time()-time1

    plt.figure()
    plt.plot(cfg.times,np.array(cfg.e[0])*cos(np.array(cfg.Omega[0])+np.array(cfg.omega[0])))

    itg.Plot_Orbit()
    mig.Plot_innerOrbit()

    itg.Plot_a()
    
    plt.show()
    
