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
    rebound.post_timestep_modifications = reboundxf.modify_elements()
    #rebound.additional_forces = reboundxf.forces()

    tauas = [0.,0.,0.,0.,0.,0.,0.,0.,-1.92e6]
    reboundxf.set_migration(tauas)
    

    # Write the semi-major axis and coordinates to a list for each planet
    #---------------------------------------------------------------------

    for i,time in enumerate(times):
        step += 1
        if step%100 == 0:
            print time
        
        rebound.integrate(time)

        for j in range(NumPlanet):
            a[j][i] = ps[j+1].calculate_orbit().a
            x[j][i] = ps[j+1].x
            y[j][i] = ps[j+1].y
    return a,x,y,times


# Call and time the function
#---------------------------------------------------------------------

time1 = time.time()

a,x,y,times = integration(100000.)

print time.time()-time1

print np.array(a[len(a)-1])


# plot
#---------------------------------------------------------------------
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif','serif':['Helvetica']})


plt.figure()
plt.xlabel('Time (y)')
plt.ylabel('Semi-Major Axis (AU)')
plt.plot(times,a[0],label='k11b')
plt.plot(times,a[1],label='k11c')
plt.plot(times,a[2],label='k11d')
plt.plot(times,a[3],label='k11e')
plt.plot(times,a[4],label='k11f')
plt.plot(times,a[5],label='k11g')
plt.plot(times,a[6],label='k11Jup1')
plt.plot(times,a[7],label='k11Jup2')
plt.legend()

plt.figure()
plt.xlabel('x')
plt.ylabel('y')
plt.title('Orbit of the 6 Kepler-11 Planets and a Jupiter in the xy-Plane')
plt.plot(x[0],y[0],'.',ms=1, label='k11b')
plt.plot(x[1],y[1],'.',ms=1, label='k11c')
plt.plot(x[2],y[2],'.',ms=1, label='k11d')
plt.plot(x[3],y[3],'.',ms=1, label='k11e')
plt.plot(x[4],y[4],'.',ms=1, label='k11f')
plt.plot(x[5],y[5],'.',ms=1, label='k11g')
plt.plot(x[6],y[6],'.',ms=1, label='k11Jup1')
plt.plot(x[7],y[7],'.',ms=1, label='k11Jup2')
plt.legend()

plt.show()