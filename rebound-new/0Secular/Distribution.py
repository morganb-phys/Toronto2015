import numpy as np
import matplotlib.pyplot as plt

def pdf(x,xmin,xmax,mean):

    pi = np.pi
    norm = 6/np.sqrt(2*pi)/(1+(xmax-mean)**2/(mean-xmin))

    if x<mean:
        pdf = norm/(mean-xmin)*np.exp(-9./2.*(x-mean)**2/(mean-xmin)**2)
    else:
        pdf = norm/(mean-xmin)*np.exp(-9./2.*(x-mean)**2/(xmax-mean)**2)

    return pdf
 
def randomvariate(pdf,n=1000,xmin=3.,xmax=12.,mean=6.):  
 
  # Calculates the minimal and maximum values of the PDF in the desired  
  # interval. The rejection method needs these values in order to work  
  # properly.  
  x=np.linspace(xmin,xmax,1000)
  y=np.zeros(len(x))
  for i,xi in list(enumerate(x)):
      y[i]=pdf(xi,xmin,xmax,mean)  
  pmin=y.min()  
  pmax=y.max() 

  plt.figure()
  plt.plot(x,y)
  plt.show()

  # Counters  
  naccept=0  
  ntrial=0  

  # Keeps generating numbers until we achieve the desired n  
  ran=np.array([]) # output list of random numbers  
 
  while naccept<n:  
      x=np.random.uniform(xmin,xmax) # x'  
      y=np.random.uniform(pmin,pmax) # y'  
      
      if y<pdf(x,xmin,xmax,mean):  
          ran = np.append(ran,x)
          naccept=naccept+1  
      ntrial=ntrial+1  
    
  ran=np.asarray(ran)  
    
  return ran,ntrial  

print randomvariate(pdf,xmin=0.995,mean=1.,xmax=1.01)
