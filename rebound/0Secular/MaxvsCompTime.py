import numpy as np
import rebound
import matplotlib.pyplot as plt
import numpy.random
import time
import reboundxf

def integration(tmax):  


    # Set up the integrator
    #---------------------------------------------------------------------

    rebound.reset()

    rebound.dt=0.1*np.sqrt(0.09100**3) 
    rebound.integrator='whfast'
    rebound.G=4.*np.pi**2
   
    # Add particles
    #---------------------------------------------------------------------

    K11 = rebound.Particle( m=0.950 )
    rebound.add( K11 )

    rebound.add( primary=K11, m=5.70467e-6, a=0.09100, e=0.045, anom=0 )  #K11b
    rebound.add( primary=K11, m=8.70713e-6, a=0.10700, e=0.026, anom=0 )  #K11c
    rebound.add( primary=K11, m=2.19179e-5, a=0.15500, e=0.004, anom=0 )  #K11d
    rebound.add( primary=K11, m=2.40197e-5, a=0.19500, e=0.012, anom=0 )  #K11e
    rebound.add( primary=K11, m=6.00492e-6, a=0.25000, e=0.013, anom=0 )  #K11f
    rebound.add( primary=K11, m=7.50615e-5, a=0.46600, e=0., anom=0 )  #K11g
    rebound.add( primary=K11, m=1.00000e-3, a=5.00000, e=0.050, anom=np.random.rand(1) )  #K11Jup1
    rebound.add( primary=K11, m=1.00000e-3, a=7.50000, e=0.036, anom=np.random.rand(1))  #K11Jup2
   
    rebound.move_to_com()
    ps = rebound.particles

    # Set up variables 
    #---------------------------------------------------------------------

    dt = rebound.dt 
    Noutputs = 10000
    NumPlanet = rebound.N-1
    times = np.linspace(0.,tmax,Noutputs)
    a = np.zeros((NumPlanet,Noutputs))
    x = np.zeros((NumPlanet,Noutputs))
    y = np.zeros((NumPlanet,Noutputs))
    step = 0

    # Set up Migration
    #rebound.post_timestep_modifications = reboundxf.modify_elements()
    rebound.additional_forces = reboundxf.forces()

    tauas = [0.,0.,0.,0.,0.,0.,0.,0.,-2.55e6]
    reboundxf.set_migration(tauas)
    

    # Write the semi-major axis and coordinates to a list for each planet
    #---------------------------------------------------------------------

    for i,time in list(enumerate(times)):
        step += 1
        if step%1000 == 0:
            print time

        rebound.integrate(time)


        for j in range(NumPlanet):
            a[j][i] = ps[j+1].calculate_orbit().a
            x[j][i] = ps[j+1].x
            y[j][i] = ps[j+1].y
    return a,x,y,times


# Call and time the function
#---------------------------------------------------------------------

comptime = []

tmax = np.logspace(0,6,num=11)

for times in tmax:
    print str(times)+':'
    time1 = time.time()
    a,x,y,times = integration(times)
    comptime.append(time.time()-time1)

print times
print comptime

plt.figure()
plt.set_xscale('log')
plt.set_yscale('log')
plt.plot(tmax,comptime)
plt.show()
