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
            print step

            # Prints an update to screen on how far along integration is.
            if step%(Noutputs/10) == 0:
                print time

            # Integrate from rebound.t (previous time) to the new time
            rebound.integrate(time,minD=0.001)
        
            # Writes the orbital elements and coordinates to cfg
            for j in range(cfg.NumPlanet):
                com = rebound.calculate_com(j)
                for OrbVar in OrbOut:
                    getattr(cfg,
                            OrbVar)[j][i] = getattr(cfg.ps[j+1].calculate_orbit(primary=com),OrbVar)
                for CoordVar in CoordOut:
                    getattr(cfg,
                            CoordVar)[j][i] = getattr(cfg.ps[j+1],CoordVar)
            if cfg.ps[1].x == 'nan':
                print 'Not a number'
                sys.exit()
                
    except rebound.CloseEncounter as e:
        print "Close encounter detected at t=%f, between particles %d and %d." % (rebound.t, e.id1, e.id2)


if __name__=="__main__":
    
    model='MigaIVa'
    integrator='ias15'

    VarOut = [['a','e','inc'],['x','y','z']]

    tmax = 100000. #years
    tau = 1.316e6
    Noutputs = 1001

    time1 = time.time()

    itg.InitRebound(integrator)

    itg.AddPlanets(model)
    rebound.add( m=1e-3, a=8.2, e=0.05, inc=89.*np.pi/180. )  #Jup1

    itg.ArrangeSys()

    itg.InitOrbitalElem(VarOut,Noutputs,tmax)

    migration(7,tau)

    CloseEnc(VarOut,Noutputs)
    
    print time.time() - time1

    #times = np.reshape(cfg.times,[1,Noutputs])
    #OrbElem = np.concatenate((times,cfg.a,cfg.e,cfg.inc))
    #Coord = np.concatenate((times,cfg.x,cfg.y,cfg.z))

    #itg.Print2File(OrbElem,
    #                 'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
    #                 +model+'_'+integrator+'_OrbElem')
    #itg.Print2File(Coord,
    #                 'tmax'+str(int(tmax))+'_tau'+str(int(tau))+'_'
    #                 +model+'_'+integrator+'_Coord')

    itg.PlotLatex()
    itg.Plot_a()
    itg.Plot_e()
    itg.Plot_inc()
    itg.Plot_Orbit()
    Plot_innerOrbit()

    plt.show()

    
