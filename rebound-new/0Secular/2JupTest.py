import numpy as np
import rebound
import matplotlib.pyplot as plt
import time

def integration(tmax):

    # Reset the integration
    #--------------------------------------------------------------

    rebound.reset()
    rebound.dt = 0.1*np.sqrt(0.09100**3)
    #rebound.integrator = 'ias15'
    rebound.integrator = 'whfast'
    rebound.G = 4.*np.pi**2
    
    # Add all the particles
    #---------------------------------------------------------------

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )
    rebound.add( primary=K11, m=5.70467e-6, a=0.091089, e=0.045, anom=1.831598 )  #K11b
    rebound.add( primary=K11, m=8.70713e-6, a=0.106522, e=0.026, anom=4.92541063 )  #K11c
    rebound.add( primary=K11, m=2.19179e-5, a=0.154241, e=0.004, anom=1.42486222 )  #K11d
    rebound.add( primary=K11, m=2.40197e-5, a=0.193937, e=0.012, anom=1.81425975 )  #K11e
    rebound.add( primary=K11, m=6.00492e-6, a=0.249489, e=0.013, anom=1. )  #K11f
    rebound.add( primary=K11, m=7.50615e-5, a=0.463918, e=0., anom=1. )  #K11g
    rebound.add( primary=K11, m=1e-3, a=5., e=0.05, anom=2. ) #K11Jup1
    rebound.add( primary=K11, m=1e-3, a=7.7, e=0.01, anom=3. )

    # Set up the integration 
    #---------------------------------------------------------------

    rebound.move_to_com()
    ps = rebound.particles


    # Set up the variables
    #--------------------------------------------------------------

    time = 0
    NumPlanet = rebound.N-1
    Planets = ['k11b','k11c','k11d','k11e','k11f','k11g']
    a = []
    for i in range(NumPlanet): 
        a.append(ps[i+1].calculate_orbit().a)
    aLowerLimit = np.array(a)*0.8
    aUpperLimit = np.array(a)*1.2
    a0 = np.array(a)
    StepSize = 0.01*np.sqrt(a[0]**3)


    # While still 'stable' and before end time run through the integrations
    #---------------------------------------------------------------

    while checkSMA(aUpperLimit, aLowerLimit, a)==True and time<tmax:

        time = rebound.t

        for j in range(NumPlanet):
            a[j] = ps[j+1].calculate_orbit().a
        check = checkSMA(aUpperLimit, aLowerLimit, a)

        rebound.integrate(time+StepSize)
        

    # Returns time it went unstable, which planet, and final eccentricity
    #---------------------------------------------------------------

    return time,np.array(a)/a0

# Function to check if Semi-Major axis has changed by 20%
#-------------------------------------------------------------------
def checkSMA(au,al,a):

    for i in range(len(a)):
        if (a[i]<al[i] or a[i]>au[i]) and i != len(a)-1:
            return False
        elif a[i]>al[i] and a[i]<au[i] and i == len(a)-1:
            return True

time1 = time.time()

t,per = integration(1000.)

comptime = time.time() - time1

print t,per, comptime
