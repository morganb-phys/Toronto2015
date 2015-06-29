import numpy as np
import rebound
from numpy import ndarray 
import sys
import cfg

import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

def Simple():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=6.00492e-6, a=0.09100, 
                 e=0.045 ) #K11b
    rebound.add( primary=K11, m=9.16540e-6, a=0.10700,
                 e=0.026 ) #K11c
    rebound.add( primary=K11, m=2.30715e-5, a=0.15500,
                 e=0.004 ) #K11d
    rebound.add( primary=K11, m=2.52839e-5, a=0.19500,
                 e=0.012 ) #K11e
    rebound.add( primary=K11, m=6.32097e-6, a=0.25000,
                 e=0.013 ) #K11f
    rebound.add( primary=K11, m=7.90121e-5, a=0.46600, 
                 e=0. ) #K11g

# Data taken from http://kepler.nasa.gov/Mission/discoveries/
def NASAPlanets():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=6.00492e-6, a=0.09100, 
                 e=0.045, inc=1.5645 ) #K11b
    rebound.add( primary=K11, m=9.16540e-6, a=0.10700,
                 e=0.026, inc=1.5636 ) #K11c
    rebound.add( primary=K11, m=2.30715e-5, a=0.15500,
                 e=0.004, inc=1.5650 ) #K11d
    rebound.add( primary=K11, m=2.52839e-5, a=0.19500,
                 e=0.012, inc=1.5514 ) #K11e
    rebound.add( primary=K11, m=6.32097e-6, a=0.25000,
                 e=0.013, inc=1.5615 ) #K11f
    rebound.add( primary=K11, m=7.90121e-5, a=0.46600, 
                 e=0., inc=1.5685 ) #K11g



# The next three models are taken from: 
# Migaszewski, C., Slonina, M., Gozdziewski, K., 2012, MNRAS, 427, 770
# Anomaly is given by 


#Model I - Only the mass of the star is fixed
def MigaI():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=1.3122e-5, a=0.091089, anom=0.4528,
                 e=0.0149, omega=-0.83, inc=1.54, Omega=0. )  #K11b
    rebound.add( primary=K11, m=2.8744e-5, a=0.106522, anom=2.9540,
                 e=0.0064, omega=-0.67, inc=1.59, Omega=0.056 )  #K11c
    rebound.add( primary=K11, m=2.7806e-5, a=0.154241, anom=-0.1530,
                 e=0.0158, omega=-2.54, inc=1.56, Omega=-0.576 )  #K11d
    rebound.add( primary=K11, m=3.3430e-5, a=0.193937, anom=1.3164,
                 e=0.0256, omega=-2.47, inc=1.55, Omega=-0.541 )  #K11e
    rebound.add( primary=K11, m=1.1248e-5, a=0.249489, anom=0.3916,
                 e=0.0180, omega=-1.91, inc=1.56, Omega=-0.559 )  #K11f
    rebound.add( primary=K11, m=5.6238e-5, a=0.463918, anom=0.6021,
                 e=0.2601, omega=-3.11, inc=1.57, Omega=-0.6 )  #K11g

# Model II - Eccentricity of Kepler-11g is fixed at 0
def MigaII():
 
    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=1.3122e-5, a=0.091087, anom=0.4528,
                 e=0.0209, omega=1.278, inc=1.54, Omega=0. )  #K11b
    rebound.add( primary=K11, m=2.8744e-5, a=0.106521, anom=2.9540,
                 e=0.0240, omega=1.528, inc=1.59, Omega=0.059 )  #K11c
    rebound.add( primary=K11, m=2.8744e-5, a=0.154233, anom=-0.1530,
                 e=0.0153, omega=2.588, inc=1.56, Omega=-0.401 )  #K11d
    rebound.add( primary=K11, m=3.2805e-5, a=0.193926, anom=1.3164,
                 e=0.0201, omega=-3.041, inc=1.55, Omega=-0.384 )  #K11e
    rebound.add( primary=K11, m=1.3747e-5, a=0.249511, anom=0.3916,
                 e=0.0081, omega=-2.090, inc=1.56, Omega=-0.436 )  #K11f
    rebound.add( primary=K11, m=9.3729e-6, a=0.463924, anom=0.6021,
                 e=0., inc=1.58, Omega=0.646 )  #K11g

# Model III - Eccentricity of Kepler-11g fixed and Longitude of ascending 
#             Node fixed for Kepler-11b and Kepler-11g
# PRETTY GOOD!!!
def MigaIII():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=1.2497e-5, a=0.091088, anom=0.4528,
                 e=0.0166, omega=0.436, inc=1.54, Omega=0. )  #K11b
    rebound.add( primary=K11, m=2.8431e-5, a=0.106519, anom=2.9540,
                 e=0.0144, omega=0.983, inc=1.55, Omega=-0.061 )  #K11c
    rebound.add( primary=K11, m=2.8431e-5, a=0.154234, anom=-0.1530,
                 e=0.0080, omega=3.142, inc=1.59, Omega=0.367 )  #K11d
    rebound.add( primary=K11, m=3.3118e-5, a=0.193924, anom=1.3164,
                 e=0.0180, omega=-2.553, inc=1.55, Omega=0.349 )  #K11e
    rebound.add( primary=K11, m=1.3434e-5, a=0.249511, anom=0.3916,
                 e=0.0112, omega=-1.751, inc=1.56, Omega=0.484 )  #K11f
    rebound.add( primary=K11, m=3.1243e-6, a=0.463917, anom=0.6021,
                 e=0., inc=1.57, Omega=0. )  #K11g


# Model IVa 
def MigaIVa():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=7.456e-6, a=0.091113, MEAN=3.122075,
                 e=0.04423, omega=0.3604, inc=1.54 )  #K11b
    rebound.add( primary=K11, m=1.070e-5, a=0.106505, MEAN=3.658224,
                 e=0.01719, omega=0.9726, inc=1.55 )  #K11c
    rebound.add( primary=K11, m=1.779e-5, a=0.154243, MEAN=0.711965,
                 e=0.00633, omega=2.4566, inc=1.59 )  #K11d
    rebound.add( primary=K11, m=3.426e-5, a=0.193940, MEAN=5.559193,
                 e=0.00258, omega=4.1323, inc=1.55 )  #K11e
    rebound.add( primary=K11, m=2.378e-5, a=0.249511, MEAN=1.598297,
                 e=0.01073, omega=6.2107, inc=1.56 )  #K11f
    rebound.add( primary=K11, m=7.952e-5, a=0.463991, MEAN=5.868932,
                 e=0., inc=1.57, Omega=0. )  #K11g


# Data taken from Mahajan N., Wu Y., 2014. ApJ, 795, 32
def Yanqin():

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=1.2e-5, a=0.091110, anom=0.4528,
                 e=0.0031, inc=1.54, omega=-1.861 )  #K11b       
    rebound.add( primary=K11, m=3.6e-5, a=0.106497, anom=2.9540,
                 e=0.0028, inc=1.55, omega=-1.170 )  #K11c      
    rebound.add( primary=K11, m=2.3e-5, a=0.154141, anom=-0.1530,
                 e=0.0195, inc=1.59, omega=-0.709 )  #K11d      
    rebound.add( primary=K11, m=2.8e-5, a=0.193972, anom=1.3164,
                 e=0.0160, inc=1.55, omega=-1.570 )  #K11e   
    rebound.add( primary=K11, m=1.0e-5, a=0.249498, anom=0.3916,
                 e=0.0125, inc=1.56, omega=-1.872 )  #K11f   
    rebound.add( primary=K11, m=2.8e-6, a=0.463981, anom=0.6021,
                 e=0. , inc=1.57 )  #K11g


# Initializes all of the important variables that are passed around functions
# To use any of these variables in another script call cfg.variable_name
def InitOrbitalElem(VarOut,Noutputs,tmax):
    
    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    #Number of Planets
    cfg.NumPlanet = rebound.N-1
    
    # Initialize Orbital Elements
    for OrbVar in OrbOut:
        setattr(cfg,OrbVar,np.zeros((cfg.NumPlanet,Noutputs)))

    for CoordVar in CoordOut:
        setattr(cfg,CoordVar,np.zeros((cfg.NumPlanet,Noutputs)))

    # The times at which it will output the orbital elements and coordinates
    cfg.times = np.linspace(0.,tmax,Noutputs)
    

# This will set up the Kepler system so that it can be integrated
def InitRebound(integrator):

    # Resets any values stored by rebound
    rebound.reset()

    # Chooses whfast as an integrator
    rebound.integrator=integrator

    # Sets G so units are AU, years, and solar masses
    rebound.G=4.*np.pi**2

def AddPlanets(model):
    # Chooses which initial conditions to use
    if model=='simple':
        Simple()
    elif model=='NASA':
        NASAPlanets()
    elif model=='MigaI':
        MigaI()
    elif model=='MigaII':
        MigaII()
    elif model=='MigaIII':
        MigaIII()
    elif model=='MigaIVa':
        MigaIVa()
    elif model=='Yanqin':
        Yanqin()
    else:
        print '\nNot a valid model. Please choose from:\nsimple, MigaI, MigaII, MigaIII, or Yanqin\n'
        sys.exit()
    
def ArrangeSys():
    # Moves the system to the center of momentum frame
    rebound.move_to_com()

    # Pointer for the planets where ps[0]=star and ps[1] is the first planet
    cfg.ps = rebound.particles
    
    # Sets the time step to tenth of Inner Period
    rebound.dt=0.01*np.sqrt(cfg.ps[1].calculate_orbit().a**3)

def InitializeInteg(model,integrator):
    InitRebound(integrator)
    AddPlanets(model)
    ArrangeSys()

# Integrates from 0 to tmax and outputs orbital elements 
# and coordinates at each time step.
def OutputOrbit(VarOut,Noutputs):

    OrbOut = VarOut[0]
    CoordOut = VarOut[1]

    # Initializes a counter for the number of steps
    step = 0

    # Runs through the times 
    for i,time in enumerate(cfg.times):
        
        #Updates the counter
        step += 1

        # Prints an update to screen on how far along integration is.
        if step%(Noutputs/10) == 0:
            print time

        # Integrate from rebound.t (previous time) to the new time
        rebound.integrate(time)
        
        # Writes the orbital elements and coordinates to cfg
        for j in range(cfg.NumPlanet):
            for OrbVar in OrbOut:
                print OrbVar
                getattr(cfg,OrbVar)[j][i] = getattr(cfg.ps[j+1].calculate_orbit(),
                                                    OrbVar)


def Print2File(array,FileName):
    Data = open('data/'+FileName+'.txt', 'a')
    [Planet,Time] = np.ma.shape(array)
    for t in range(Time):
        output = ''
        for p in range(Planet):
            output = output+str(array[p,t])+' '
        Data.write(output+'\n')
    Data.write('\n')
    Data.close()

def PlotvTime(array,ylabel=None):
    plt.figure()
    plt.xlabel('Time (y)')
    plt.ylabel(ylabel)
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    Colours = ['r','y','g','b','c','k','m','#fa8a29']
    for Planet in range(cfg.NumPlanet):
        plt.plot(cfg.times,array[Planet],
                 color=Colours[Planet],label=Names[Planet])
    plt.legend()

# Plots Semi-Major axis vs. Time
def Plot_a():
    PlotvTime(cfg.a,'Semi-Major Axis (AU)')
    
def Plot_inc():
    PlotvTime(cfg.inc*180/np.pi,'Inclination ( $^{\circ}$ )')

# Plots Eccentricity vs. Time
def Plot_e():
    PlotvTime(cfg.e,'Eccentricity')

# Plot in xy-Plane
def Plot_xy():
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    Colours = ['r','y','g','b','c','k','m','#fa8a29']
    for Planet in range(cfg.NumPlanet):
        plt.plot(cfg.x[Planet],cfg.y[Planet],'.',
                 color=Colours[Planet],label=Names[Planet])
    plt.legend()

# 3D plot of the orbits
def Plot_Orbit():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel('x')
    plt.ylabel('y')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    Colours = ['r','y','g','b','c','k','m','#fa8a29']
    print cfg.NumPlanet
    for Planet in range(cfg.NumPlanet):
        ax.scatter(cfg.x[Planet],cfg.y[Planet],cfg.z[Planet],'o',s=5,
                   c=Colours[Planet],edgecolors='none', label=Names[Planet])
    plt.legend()


# Writes plot text in LateX
def PlotLatex():
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Helvetica']})

if __name__=="__main__":
    
    model='MigaIVa'
    integrator = 'ias15'

    tmax = 100. # Final time in years
    Noutputs = 1001 # Number of outputs
    
    OutputVar = [['a','e'],['x','y','z']]

    # Begins timing the integration
    time1 = time.time()

    InitializeInteg(model,integrator)
    InitOrbitalElem(OutputVar,Noutputs,tmax)
    OutputOrbit(OutputVar,Noutputs)

    # Prints total computational time to screen
    print time.time() - time1

    times = np.reshape(cfg.times,[1,Noutputs])
    OrbElem = np.concatenate((times,cfg.a,cfg.e))
    Coord = np.concatenate((times,cfg.x,cfg.y,cfg.z))

    #Print2File(OrbElem,str(int(tmax))+model+integrator+'_OrbElem')
    #Print2File(OrbElem,str(int(tmax))+model+integrator+'_Coord')

    PlotLatex()
    Plot_a()
    Plot_e()
    names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    Plot_Orbit()
    Plot_xy()

    # Show the plots
    plt.show()

    
