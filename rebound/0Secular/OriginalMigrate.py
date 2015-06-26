import numpy as np
import rebound
import matplotlib.pyplot as plt
import numpy.random
import time
import reboundxf
from numpy import ndarray
from mpl_toolkits.mplot3d import Axes3D


def integration(tmax,tau):
    
    
    # Set up the integrator
    #---------------------------------------------------------------------
    rebound.reset()
    rebound.dt=0.01*np.sqrt(0.09111**3)
    rebound.integrator='ias15'
    rebound.G=4.*np.pi**2
    
    
    # Add particles
    #---------------------------------------------------------------------

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
    rebound.add( m=1e-3, a=5.0, e=0.05, inc=89*np.pi/180 )  #Jup1
    rebound.add( m=1e-3, a=7.5, e=0.03, inc=88.5*np.pi/180 )  #Jup2
    
    rebound.move_to_com()
    ps = rebound.particles
    
    
    # Set up variables
    #---------------------------------------------------------------------
    dt = rebound.dt
    Noutputs = 10000
    NumPlanet = rebound.N-1
    times = np.linspace(0.,tmax,Noutputs)
    a = np.zeros((NumPlanet,Noutputs))
    e = np.zeros((NumPlanet,Noutputs))
    x = np.zeros((NumPlanet,Noutputs))
    y = np.zeros((NumPlanet,Noutputs))
    z = np.zeros((NumPlanet,Noutputs))
    step = 0
    
    
    # Set up Migration
    rebound.post_timestep_modifications = reboundxf.modify_elements()
    tauas = [0.,0.,0.,0.,0.,0.,0.,0.,tau]
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
            e[j][i] = ps[j+1].calculate_orbit().e
            x[j][i] = ps[j+1].x
            y[j][i] = ps[j+1].y
            z[j][i] = ps[j+1].z

    return a,e,x,y,z,times
        


# Call and time the function
#---------------------------------------------------------------------
tmax = 1000000.
tau = -2.55e7

time1 = time.time()
a,e,x,y,z,times = integration(tmax,tau)
print time.time()-time1

print len(a)
print np.array(a[len(a)-1])

'''

# Write to File
#---------------------------------------------------------------------
FileName = 'data/tmax'+str(int(tmax))+'_tau'+str(int(tau))+'.txt'
Data = open(FileName,'a')
for i,time in list(enumerate(times,start=0)):
    for j in (range(len(a))):
        Data.write('\t'+' '+str(a[j][i])+' '+str(e[j][i])+' '+str(x[j][i])+' '+str(y[j][i])+'\n')
        
Data.write('\n')
Data.close()
'''

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
#plt.plot(times,a[6],label='k11Jup1')
#plt.plot(times,a[7],label='k11Jup2')
plt.legend()

plt.figure()
plt.xlabel('Time (y)')
plt.ylabel('Eccentricity')
plt.plot(times,e[0],label='k11b')
plt.plot(times,e[1],label='k11c')
plt.plot(times,e[2],label='k11d')
plt.plot(times,e[3],label='k11e')
plt.plot(times,e[4],label='k11f')
plt.plot(times,e[5],label='k11g')
#plt.plot(times,a[6],label='k11Jup1')
#plt.plot(times,a[7],label='k11Jup2')
plt.legend()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Orbit of the 6 Kepler-11 Planets and a Jupiter in the xy-Plane')
ax.scatter(x[0],y[0],z[0],'o',s=5,c='r',edgecolors='none', label='k11b')
ax.scatter(x[1],y[1],z[1],'o',s=5,c='c',edgecolors='none', label='k11c')
ax.scatter(x[2],y[2],z[2],'o',s=5,c='y',edgecolors='none', label='k11d')
ax.scatter(x[3],y[3],z[3],'o',s=5,c='g',edgecolors='none', label='k11e')
ax.scatter(x[4],y[4],z[4],'o',s=5,c='b',edgecolors='none', label='k11f')
ax.scatter(x[5],y[5],z[5],'o',s=5,c='m',edgecolors='none', label='k11g')
plt.legend()


plt.show()
