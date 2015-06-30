import numpy as np
import rebound
import matplotlib.pyplot as plt
import numpy.random
import time
import reboundxf
from scipy import optimize

# Call and time the function
#---------------------------------------------------------------------

comptime = []

tmax = np.logspace(0,4,num=20)

for times in tmax:
    rebound.reset()

    # Set up the integrator
    #---------------------------------------------------------------------

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
    print str(times)+':'
    time1 = time.time()
    rebound.integrate(times)
    comptime.append(time.time()-time1)

para, pcov = optimize.curve_fit(lambda x,a,b: a*x**b,tmax,comptime)

print para

x = np.logspace(0,4,100)
y = para[0]*x**para[1]

plt.figure()
plt.loglog(tmax,comptime,'.')
plt.loglog(x,y)
plt.show()
