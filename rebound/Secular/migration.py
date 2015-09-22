import rebound
import reboundxf

import numpy as np
from numpy import ndarray
 
import sys
import time

import integration.integrator as itg
import Evection.Phi as phi
import integration.cfg as cfg
import integration.destab as dtb

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D




def migration(planet,tau):
    rebound.post_timestep_modifications = reboundxf.modify_elements()

    #rebound.additional_forces = reboundxf.forces()

    tauas = np.zeros(cfg.NumPlanet+1)
    tauas[planet] = tau
    reboundxf.set_migration(tauas)


def Plot_innerOrbit():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel('x')
    plt.ylabel('y')
    Colours = ['r','y','g','b','c','k']
    for Planet in range(6):
        ax.scatter(cfg.x[Planet],cfg.y[Planet],cfg.z[Planet],'o',s=5,
                   c=Colours[Planet],edgecolors='none', label=cfg.Names[Planet])
    plt.legend()
    

def CloseEnc(VarOut,Noutputs):

    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    # Initializes a counter for the number of steps
    step = 0

    # Runs through the times 
    try:
        for i,time in enumerate(cfg.times):
        
            #Updates the counter
            step += 1

            # Prints an update to screen on how far along integration is.
            if step%(Noutputs/10) == 1:
                print time

            # Integrate from rebound.t (previous time) to the new time
            rebound.integrate(time,minD=0.001)
        
            # Writes the orbital elements and coordinates to cfg
            for j in range(cfg.NumPlanet):
                com = rebound.calculate_com(j)
                for OrbVar in OrbOut:
                    getattr(cfg,
                            OrbVar)[j][i] = getattr(cfg.ps[j+1].calculate_orbit(),OrbVar)
                for CoordVar in CoordOut:
                    getattr(cfg,
                            CoordVar)[j][i] = getattr(cfg.ps[j+1],CoordVar)
            if cfg.ps[1].x == 'nan':
                print 'Not a number'
                sys.exit()
                
    except rebound.CloseEncounter as e:
        print "Close encounter detected at t=%f, between particles %d and %d." % (rebound.t, e.id1, e.id2)

def migr(para):

    isma = para[0]
    fsma = para[1]
    tmax = para[2]

    model='MigaIVa'
    integrator='ias15'

    VarOut = [['a','e','Omega','omega','l'],['x','y','z']]

    tau = tmax/np.log(isma/fsma)
    Nout = 1001

    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    rebound.add( m=1e-3, a=isma, e=0.05, inc=89.*np.pi/180. )  #Jup1

    itg.ArrangeSys()

    dtb.InitOrbElem2(VarOut,Nout,tmax)

    migration(7,tau)

    dtb.Destab(tmax,VarOut,Nout)

    Nout = len(cfg.times)

    itg.PlotLatex()
    Phi = phi.calcPhi(Nout)
    phi.Plot_phi(isma,Phi,'Phi')
    
    times = np.reshape(cfg.times,[1,Nout])
    OrbElem = np.concatenate((times,cfg.a,cfg.e))
    Coord = np.concatenate((times,cfg.x,cfg.y,cfg.z))
    Phi = np.concatenate((times,cfg.Omega,cfg.omega,cfg.l))

    itg.Print2File(OrbElem,
                     'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
                     +str(int(isma))+'_OrbElem')
    itg.Print2File(Coord,
                     'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
                     +str(int(isma))+'_Coord')
    itg.Print2File(Coord,
                     'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
                     +str(int(isma))+'_Phi')

    
    strSMA = str(isma).replace('.','_')

    itg.PlotvTime(cfg.a,'Semi-Major Axis (AU)',Title=str(isma),
                  save=True,File='plots/'+strSMA+'/sma.png')
    itg.PlotvTime(cfg.e,'Eccentricity',Title=str(isma),save=True,
                  File='plots/'+strSMA+'/ecc.png')
    
    plt.figure()
    plt.xlabel('Time (yr)')
    plt.ylabel('Semi-Major Axis (AU)')
    plt.title(str(isma))
    
    for Planet in range(cfg.NumPlanet-1):
        plt.plot(cfg.times,cfg.a[Planet],
                 color=cfg.Colours[Planet],label=cfg.Names[Planet])
    plt.legend()
    plt.savefig('plots/'+strSMA+'/sma_zoom.png')

    plt.show()
    

if __name__=="__main__":

    #ismas = [5.10,5.20,5.55,6.10,5.10,6.37,6.73,7.17,7.53,7.97,8.62,
    #         10.22,12.18,12.84,14.31,17.77,18.31]
    #fsmas = [5.04,5.16,5.45,6.00,5.00,6.33,6.69,7.13,7.47,7.93,8.58,
    #         10.18,12.12,12.78,14.28,17.73,18.27]

    tmax = 10000.

    ismas = [5.2]
    fsmas = [5.16]

    parameters = [(ismas[i],fsmas[i],tmax) for i in range(len(ismas))]

    pool = rebound.InterruptiblePool()
    pool.map(migr,parameters)
    

