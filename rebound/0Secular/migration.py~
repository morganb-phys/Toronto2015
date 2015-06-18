import numpy as np
import rebound
from numpy import ndarray 
import sys
import integrator as integ
import cfg
import reboundxf
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D


def migration(planet,tau):
    rebound.post_timestep_modifications = reboundxf.modify_elements()

    tauas = np.zeros(cfg.NumPlanet+1)
    tauas[planet] = tau
    reboundxf.set_migration(tauas)
    

if __name__=="__main__":
    
    model='simple'
    tmax = 100. #years
    tau = -2.55e-6
    Noutputs = 1000
    time1 = time.time()

    integ.InitializeInteg(tmax,model,Noutputs)

    dt = rebound.dt
    integ.InitOrbitalElem(Noutputs,tmax)

    migration(6,-2.55e2)

    integ.OutputOrbit()

    print time.time() - time1


    integ.PlotLatex()
    integ.Plot_a()
    integ.Plot_e()
    integ.Plot_Orbit()

    plt.show()

    
