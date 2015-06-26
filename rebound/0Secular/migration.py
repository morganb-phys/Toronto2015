import numpy as np
import rebound
from numpy import ndarray 
import sys
import integration.integrator as itg
import integration.cfg as cfg
import reboundxf
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D


def migration(planet,tau):
    rebound.post_timestep_modifications = reboundxf.modify_elements()

    tauas = np.zeros(cfg.NumPlanet+1)
    tauas[planet] = tau
    reboundxf.set_migration(tauas)

def Plot_innerOrbit():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel('x')
    plt.ylabel('y')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g']
    Colours = ['r','y','g','b','c','k']
    for Planet in range(cfg.NumPlanet-2):
        ax.scatter(cfg.x[Planet],cfg.y[Planet],cfg.z[Planet],'o',s=5,
                   c=Colours[Planet],edgecolors='none', label=Names[Planet])
    plt.legend()
    

if __name__=="__main__":
    
    model='Yanqin'
    integrator='ias15'

    tmax = 100000. #years
    tau = 1.70e5
    Noutputs = 1001

    time1 = time.time()

    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    rebound.add( m=1e-3, a=9.0, e=0.05, inc=89.*np.pi/180. )  #Jup1
    #rebound.add( m=1e-3, a=7.5, e=0.03, inc=88.5*np.pi/180. )  #Jup2

    itg.ArrangeSys()

    itg.InitOrbitalElem(Noutputs,tmax)

    migration(7,tau)

    itg.OutputOrbit(Noutputs)
    
    print time.time() - time1

    print cfg.a[:,Noutputs-1]


    times = np.reshape(cfg.times,[1,Noutputs])
    OrbElem = np.concatenate((times,cfg.a,cfg.e,cfg.inc))
    Coord = np.concatenate((times,cfg.x,cfg.y,cfg.z))

    itg.Print2File(OrbElem,
                     'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
                     +model+'_'+integrator+'_OrbElem')
    itg.Print2File(Coord,
                     'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
                     +model+'_'+integrator+'_Coord')

    

    itg.PlotLatex()
    itg.Plot_a()
    itg.Plot_e()
    itg.Plot_inc()
    itg.Plot_Orbit()
    Plot_innerOrbit()

    plt.show()

    
