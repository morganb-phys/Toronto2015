import numpy as np
import rebound
from numpy import ndarray 
import sys
import cfg

import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D


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

def MigaI():
    # Model I
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


def MigaII():
    # Model II
    rebound.add( primary=K11, m=1.3122e-5, a=0.091087, anom=0.4528,
                 e=0.0209, omega=1.278, inc=1.54)#, Omega=0. )  #K11b
    rebound.add( primary=K11, m=2.8744e-5, a=0.106521, anom=2.9540,
                 e=0.0240, omega=1.528, inc=1.59)#, Omega=0.059 )  #K11c
    rebound.add( primary=K11, m=2.8744e-5, a=0.154233, anom=-0.1530,
                 e=0.0153, omega=2.588, inc=1.56)#, Omega=-0.401 )  #K11d
    rebound.add( primary=K11, m=3.2805e-5, a=0.193926, anom=1.3164,
                 e=0.0201, omega=-3.041, inc=1.55)#, Omega=-0.384 )  #K11e
    rebound.add( primary=K11, m=1.3747e-5, a=0.249511, anom=0.3916,
                 e=0.0081, omega=-2.090, inc=1.56)#, Omega=-0.436 )  #K11f
    rebound.add( primary=K11, m=9.3729e-6, a=0.463924, anom=0.6021,
                 e=0., inc=1.58)#, Omega=0.646 )  #K11g


def MigaIII():

    # Model III PRETTY GOOD!!!
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



def Yanqin():
    # Yanqin... OKAY?
    rebound.add( primary=K11, m=1.2e-5, a=0.091110, anom=0.4528,
                 e=0.0031, inc=1.54 )  #K11b       , omega=-1.861
    rebound.add( primary=K11, m=3.6e-5, a=0.106497, anom=2.9540,
                 e=0.0028, inc=1.55 )  #K11c      , omega=-1.170
    rebound.add( primary=K11, m=2.3e-5, a=0.154141, anom=-0.1530,
                 e=0.0195, inc=1.59 )  #K11d      , omega=-0.709
    rebound.add( primary=K11, m=2.8e-5, a=0.193972, anom=1.3164,
                 e=0.0160, inc=1.55 )  #K11e   , omega=-1.570
    rebound.add( primary=K11, m=1.0e-5, a=0.249498, anom=0.3916,
                 e=0.0125, inc=1.56 )  #K11f   , omega=-1.872
    rebound.add( primary=K11, m=2.8e-6, a=0.463981, anom=0.6021,
                 e=0. , inc=1.57 )  #K11g



def InitOrbitalElem(Noutputs,tmax):
    cfg.NumPlanet = rebound.N-1
    cfg.a = np.zeros((cfg.NumPlanet,Noutputs))
    cfg.e = np.zeros((cfg.NumPlanet,Noutputs))
    cfg.x = np.zeros((cfg.NumPlanet,Noutputs))
    cfg.y = np.zeros((cfg.NumPlanet,Noutputs))
    cfg.z = np.zeros((cfg.NumPlanet,Noutputs))
    cfg.times = np.linspace(0.,tmax,Noutputs)
    


def InitializeInteg(tmax,model,Noutputs):

    # Resets any values stored by rebound
    rebound.reset()

    # Sets the time step to tenth of Inner Period
    rebound.dt=0.1*np.sqrt(0.09100**3)

    # Chooses whfast as an integrator
    rebound.integrator='whfast'

    # Sets G so units are AU, years, and solar masses
    rebound.G=4.*np.pi**2


    if model=='simple':
        # Add particles with values from http://kepler.nasa.gov/Mission/discoveries/
        NASAPlanets()
    elif model=='MigaI':
        MigaI()
    elif model=='MigaII':
        MigaII()
    elif model=='MigaIII':
        MigaIII()
    elif model=='Yanqin':
        Yanqin()
    else:
        print '\nNot a valid model. Please choose from:\nsimple, MigaI, MigaII, MigaIII, or Yanqin\n'
        sys.exit()
    
    rebound.move_to_com()

    global ps; ps = rebound.particles


def OutputOrbit():

    step = 0
    for i,time in enumerate(cfg.times):
        step += 1
        if step%100 == 0:
            print time

        rebound.integrate(time)
        
        for j in range(cfg.NumPlanet):
            cfg.a[j][i] = ps[j+1].calculate_orbit().a
            cfg.e[j][i] = ps[j+1].calculate_orbit().e
            cfg.x[j][i] = ps[j+1].x
            cfg.y[j][i] = ps[j+1].y
            cfg.z[j][i] = ps[j+1].z

def Plot_a():
    plt.figure()
    plt.xlabel('Time (y)')
    plt.ylabel('Semi-Major Axis (AU)')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    for Planet in range(cfg.NumPlanet):
        plt.plot(cfg.times,cfg.a[Planet],label=Names[Planet])
    plt.legend()
    


def Plot_e():
    plt.figure()
    plt.xlabel('Time (y)')
    plt.ylabel('Eccentricity')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    for Planet in range(cfg.NumPlanet):
        plt.plot(cfg.times,cfg.e[Planet],label=Names[Planet])
    plt.legend()

def Plot_Orbit():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.xlabel('x')
    plt.ylabel('y')
    Names = ['k11b','k11c','k11d','k11e','k11f','k11g','Jup1','Jup2']
    Colours = ['r','y','g','b','c','k','m']
    for Planet in range(cfg.NumPlanet):
        ax.scatter(cfg.x[Planet],cfg.y[Planet],cfg.z[Planet],'o',s=5,
                   c=Colours[Planet],edgecolors='none', label=Names[Planet])
    plt.legend()

def PlotLatex():
    plt.rc('text', usetex=True)
    plt.rc('font', **{'family':'serif','serif':['Helvetica']})

if __name__=="__main__":
    
    model='simple'
    tmax = 100. #years
    Noutputs = 1000
    time1 = time.time()

    InitializeInteg(tmax,model,Noutputs)

    dt = rebound.dt
    InitOrbitalElem(Noutputs,tmax)

    OutputOrbit()

    print time.time() - time1


    PlotLatex()
    Plot_a()
    Plot_e()
    Plot_Orbit()

    plt.show()

    
