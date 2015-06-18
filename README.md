# Toronto2015

####Project
The Kepler satellite has found several compact solar systems containing earth mass planets. Yet, these systems do not appear to conatain any giant planets. This purpose of this project is to investigate whether or not this lack of giant planets is due to the bias in Kepler to find planets closer to the star or if we don't find them because their presence would destabilize the system. 

Destabilization due to an MMR has been ruled out because for a giant planet to be in a 2:1 resonance with the outermost planet would require it to have a semi-major axis of ~0.7 AU at which point it would have been detected by Kepler. Therefore, the most likely cause for destabilization would be a secular perturbation due to two giant planets further out, ~5 AU, with an eigenfrequency resonant with an eigenfrequency of one of the inner planets.

The majority of the code uses [Rebound](https://github.com/hannorein/rebound) which is an open-source multi-purpose N-bode code developed by:

* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>
* Akihiko Fujii, National Astronomical Observatory of Japan/University of Tokyo, Tokyo, <akihiko.fujii@nao.ac.jp>
* Dan Tamayo, University of Toronto, <dtamayo@cita.utoronto.ca>

#####integrator.py

Example code used for a simple integration: 
```python
  # Choose the model for your initial conditions
  model='simple'
    
  # Choose the final time in years
  tmax = 100. # Final time in years

  # Choose the Number of Outputs
  # It will save semi-major axis, eccentricity, and coordinates every (tmax/Noutputs) years
  Noutputs = 10000
    
  # Initializes the Rebound code, adds the planets and creates a pointer to the planets. 
  InitializeInteg(model)
  
  # Builds Numpy arrays to hold the values for the orbital elements, the coordinates and time 
  InitOrbitalElem(Noutputs,tmax)
  
  #Integrates from t=0 to t=tmax
  OutputOrbit(Noutputs)

  # Allows plot captions to be written in LateX
  PlotLatex()
  
  #Plots semi-major axis vs. time for all planets
  Plot_a()
  
  #Plots eccentricity vs. time for all the planets
  Plot_e()
  
  #Plots position of the planets first in 3D and then in the xy-plane
  Plot_Orbit()
  Plot_xy()

  # Show the plots
    plt.show()
```


*Simple()* 



#####Secular
