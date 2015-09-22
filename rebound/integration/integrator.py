import numpy as np
import rebound
from numpy import ndarray 
import sys
import cfg

import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

def Simple(sim):

    K11 = rebound.Particle( m=0.950 )
    sim.add( K11 )
    sim.add( primary=K11, m=6.00492e-6, a=0.09100, 
                 e=0.045 ) #K11b
    sim.add( primary=K11, m=9.16540e-6, a=0.10700,
                 e=0.026 ) #K11c
    sim.add( primary=K11, m=2.30715e-5, a=0.15500,
                 e=0.004 ) #K11d
    sim.add( primary=K11, m=2.52839e-5, a=0.19500,
                 e=0.012 ) #K11e
    sim.add( primary=K11, m=6.32097e-6, a=0.25000,
                 e=0.013 ) #K11f
    sim.add( primary=K11, m=7.90121e-5, a=0.46600, 
                 e=0. ) #K11g

# The next model is taken from: 
# Migaszewski, C., Slonina, M., Gozdziewski, K., 2012, MNRAS, 427, 770

# Model IVa 
def MigaIVa(sim):

    K11 = rebound.Particle( m=0.950 )
    sim.add( K11 )
    sim.add( primary=K11, m=7.456e-6, a=0.091113, MEAN=3.123995,
                 e=0.04423, omega=0.3604, inc=1.5558 )  #K11b
    sim.add( primary=K11, m=1.070e-5, a=0.106505, MEAN=3.658224,
                 e=0.01719, omega=0.9726, inc=1.592 )  #K11c
    sim.add( primary=K11, m=1.779e-5, a=0.154243, MEAN=0.711965,
                 e=0.00633, omega=2.4566, inc=1.5591 )  #K11d
    sim.add( primary=K11, m=3.426e-5, a=0.193940, MEAN=5.559193,
                 e=0.00258, omega=4.1323, inc=1.5505 )  #K11e
    sim.add( primary=K11, m=2.378e-5, a=0.249511, MEAN=1.598297,
                 e=0.01073, omega=6.2107, inc=1.5602 )  #K11f
    sim.add( primary=K11, m=7.952e-5, a=0.463991, MEAN=5.868932,
                 e=0., inc=1.5668 )  #K11g


# Data taken from Mahajan N., Wu Y., 2014. ApJ, 795, 32
def Yanqin(sim):

    K11 = rebound.Particle( m=0.950 )
    sim.add( K11 )
    sim.add( primary=K11, m=1.2e-5, a=0.091110, MEAN=0.4555,
                 e=0.0031, inc=1.54, omega=-1.861 )  #K11b       
    sim.add( primary=K11, m=3.6e-5, a=0.106497, anom=2.9514,
                 e=0.0028, inc=1.55, omega=-1.170 )  #K11c      
    sim.add( primary=K11, m=2.3e-5, a=0.154141, anom=6.1244,
                 e=0.0195, inc=1.59, omega=-0.709 )  #K11d      
    sim.add( primary=K11, m=2.8e-5, a=0.193972, anom=1.3474,
                 e=0.0160, inc=1.55, omega=-1.570 )  #K11e   
    sim.add( primary=K11, m=1.0e-5, a=0.249498, anom=0.4014,
                 e=0.0125, inc=1.56, omega=-1.872 )  #K11f   
    sim.add( primary=K11, m=2.8e-6, a=0.463981, anom=0.6021,
                 e=0. , inc=1.57 )  #K11g


# Initializes all of the important variables that are passed around functions
# To use any of these variables in another script call cfg.variable_name
def InitOrbitalElem(VarOut,Noutputs,tmax,sim):
    
    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    #Number of Planets
    cfg.NumPlanet = sim.N-1
    print 'Number of planets',cfg.NumPlanet

    # Initialize Orbital Elements
    for OrbVar in OrbOut:
        setattr(cfg,OrbVar,np.zeros((cfg.NumPlanet,Noutputs)))

    for CoordVar in CoordOut:
        setattr(cfg,CoordVar,np.zeros((cfg.NumPlanet,Noutputs)))

    # The times at which it will output the orbital elements and coordinates
    cfg.times = np.linspace(0.,tmax,Noutputs)
    

# This will set up the Kepler system so that it can be integrated
def InitRebound(integrator,sim):

    # Chooses an integrator
    sim.integrator=integrator

    # Sets G so units are AU, years, and solar masses
    sim.G=4.*np.pi**2

def AddPlanets(model,sim):
    # Chooses which initial conditions to use
    if model=='simple':
        Simple(sim)
    elif model=='NASA':
        NASAPlanets(sim)
    elif model=='MigaI':
        MigaI(sim)
    elif model=='MigaII':
        MigaII(sim)
    elif model=='MigaIII':
        MigaIII(sim)
    elif model=='MigaIVa':
        MigaIVa(sim)
    elif model=='Yanqin':
        Yanqin(sim)
    else:
        print '\nNot a valid model. Please choose from:\nsimple, MigaI, MigaII, MigaIII, or Yanqin\n'
        sys.exit()
    
def ArrangeSys(sim):
    # Moves the system to the center of momentum frame
    sim.move_to_com()

    # Pointer for the planets where ps[0]=star and ps[1] is the first planet
    cfg.ps = sim.particles

    # Sets the time step to tenth of Inner Period
    sim.dt=0.01*np.sqrt(sim.calculate_orbits()[0].a**3)

def InitializeInteg(model,integrator,sim):
    InitRebound(integrator,sim)
    AddPlanets(model,sim)
    ArrangeSys(sim)

# Integrates from 0 to tmax and outputs orbital elements 
# and coordinates at each time step.
def OutputOrbit(VarOut,Noutputs,sim,TimeUpdate=True):

    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    # Initializes a counter for the number of steps
    step = 0

    # Runs through the times 
    for i,time in enumerate(cfg.times):

        if TimeUpdate==True:
            #Updates the counter
            step += 1

            # Prints an update to screen on how far along integration is.
            if step%(Noutputs/10) == 1:
                print time

        # Integrate from rebound.t (previous time) to the new time
        sim.integrate(time)
        
        '''
        print 'Time: ', time
        for i in range(cfg.NumPlanet):
            print cfg.ps[i+1].calculate_orbit().e
        os = rebound.calculate_orbits()
        '''

        # Writes the orbital elements and coordinates to cfg
        for j in range(cfg.NumPlanet):
            
            Orbits = sim.calculate_orbits()

            for OrbVar in OrbOut:
                getattr(cfg,OrbVar)[j][i] = getattr(Orbits[j], OrbVar)

            for CoordVar in CoordOut:
                getattr(cfg,
                        CoordVar)[j][i] = getattr(cfg.ps[j+1],CoordVar)
        

# prints entered array to text file with FileName
def Print2File(array,FileName):

    Data = open('data/'+FileName+'.txt', 'a')
    
    # Finds the number of Planets and Times
    [Planet,Time] = np.ma.shape(array)
    
    for t in range(Time):
        output = ''
    
        # Adds the next element in teh array to the t-th line
        for p in range(Planet):
            output = output+str(array[p,t])+' '

            # Writes the entire line to file
        Data.write(output+'\n')
    Data.write('\n')
    Data.close()


# Plots an array vs. cfg.times
def PlotvTime(array,ylabel=None,Title=None,save=False,File=None):

    plt.figure()
    plt.xlabel('Time (y)')
    plt.ylabel(ylabel)
    plt.title(Title)

    for Planet in range(cfg.NumPlanet):

        plt.plot(cfg.times,array[Planet],
                 color=cfg.Colours[Planet],label=cfg.Names[Planet])
    plt.legend()
    
    if save==True:
        plt.savefig(File)


# Plots Semi-Major axis vs. Time
def Plot_a(Title=None):
    PlotvTime(cfg.a,'Semi-Major Axis (AU)',Title=Title)
    
def Plot_inc(Title=None):
    PlotvTime(cfg.inc*180/np.pi,'Inclination ( $^{\circ}$ )',Title=Title)

# Plots Eccentricity vs. Time
def Plot_e(Title=None):
    PlotvTime(cfg.e,'Eccentricity',Title=Title)

# Plot in xy-Plane
def Plot_xy(Title=None):

    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(Title)

    for Planet in range(cfg.NumPlanet):
        plt.plot(cfg.x[Planet],cfg.y[Planet],'.',
                 color=cfg.Colours[Planet],label=cfg.Names[Planet])
    plt.legend()


# 3D plot of the orbits
def Plot_Orbit(Title=None):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(Title)

    # edge colour dominates at this size of marker
    for Planet in range(cfg.NumPlanet):
        ax.scatter(cfg.x[Planet],cfg.y[Planet],cfg.z[Planet],'o',s=5,
                   c=cfg.Colours[Planet],edgecolors='none',
                   label=cfg.Names[Planet])
    plt.legend()


# Writes plot text in LateX
def PlotLatex():
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Helvetica']})

if __name__=="__main__":
    
    model='MigaIVa'
    integrator = 'ias15'

    tmax = 1000. # Final time in years
    Noutputs = 1001 # Number of outputs
    
    OutputVar = [['a','e'],['x','y','z']]

    # Begins timing the integration
    time1 = time.time()

    InitializeInteg(model,integrator,sim)
    InitOrbitalElem(OutputVar,Noutputs,tmax,sim)
    OutputOrbit(OutputVar,Noutputs,sim)

    # Prints total computational time to screen
    print time.time() - time1

    times = np.reshape(cfg.times,[1,Noutputs])
    OrbElem = np.concatenate((times,cfg.a,cfg.e))
    Coord = np.concatenate((times,cfg.x,cfg.y,cfg.z))

    Print2File(OrbElem,str(int(tmax))+model+integrator+'_OrbElem')
    Print2File(Coord,str(int(tmax))+model+integrator+'_Coord')

    PlotLatex()
    Plot_a()
    Plot_e()

    Plot_Orbit()
    Plot_xy()

    # Show the plots
    plt.show()

    
