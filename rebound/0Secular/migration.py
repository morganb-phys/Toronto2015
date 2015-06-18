import numpy as np
import rebound
from numpy import ndarray 
import sys
import integration.integrator as integ
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
    

if __name__=="__main__":
    
    model='simple'
    tmax = 100. #years
    tau = -2.55e3
    Noutputs = 10000
    time1 = time.time()

    integ.InitRebound()

    integ.AddPlanets(model)
    rebound.add( m=1e-3, a=5.0, e=0.05 )  #Jup1
    rebound.add( m=1e-3, a=7.5, e=0.03 )  #Jup2

    integ.ArrangeSys()

    integ.InitOrbitalElem(Noutputs,tmax)

    migration(8,tau)

    integ.OutputOrbit(Noutputs)

    print time.time() - time1

    integ.PlotLatex()
    integ.Plot_a()
    integ.Plot_e()
    integ.Plot_xy()

    plt.show()

    
