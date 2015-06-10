import numpy as np
import rebound
from string import maketrans
import matplotlib.pyplot as plt
import scipy.stats

NumInit = 10000
NumBins = 20



# Read in the data from triangular distribution
File = open('data/Freq_triangular10000.txt')
Freq1 = np.zeros([NumInit,6])


for i,line in list(enumerate(File)):
    if i<NumInit:
        line = line.split()
        for j in range(6):
            Freq1[i,j] = float(line[j])
File.close()


# Read in the data from uniform distribution
File = open('data/Freq_Gaussian10000.txt')
Freq2 = np.zeros([NumInit,6])


for i,line in list(enumerate(File)):
    if i<NumInit:
        line = line.split()
        for j in range(6):
            Freq2[i,j] = float(line[j])

File.close()

print 'Triangular:', Freq1.max(), Freq1.min()
print 'Uniform:', Freq2.max(), Freq2.min()

plt.figure()
plt.hist(Freq1[:,0], bins = 50)
plt.hist(Freq2[:,0], bins = 50, fill=False)

plt.figure()
plt.hist(Freq1[:,1], bins = 50)
plt.hist(Freq2[:,1], bins = 50, fill=False)

plt.figure()
plt.hist(Freq1[:,2], bins = 50)
plt.hist(Freq2[:,2], bins = 50, fill=False)

plt.figure()
plt.hist(Freq1[:,3], bins = 50)
plt.hist(Freq2[:,3], bins = 50, fill=False)

plt.figure()
plt.hist(Freq1[:,4], bins = 50)
plt.hist(Freq2[:,4], bins = 50, fill=False)


plt.figure()
plt.hist(Freq1[:,5], bins = 50)
plt.hist(Freq2[:,5], bins = 50, fill=False)


plt.show()
'''
# Read in the data from whfast
wh001 = open('data/'+str(NumInit)+'angles_whfast.txt')
Twh001 = np.zeros(NumInit)
Fwh001 = np.zeros(NumInit)

for i,line in list(enumerate(wh001)):
    if i<NumInit:
        line = line.translate(maketrans("",""), '[]')
        line = line.split()
        Twh001[i] = float(line[0])
        #Fwh001[i] = [float(anom) for anom in line[i:]]



# Make a histogram of whfast data and a CDF
HistWH001, BinWH001 = np.histogram(Twh001,bins=NumBins,range=[0,Twh001.max()])
CDFWH001 = np.cumsum(HistWH001)/float(NumInit)
widwh001 = BinWH001[0]-BinWH001[1]



# Make a histogram of whfast data and a CDF
HistWH01, BinWH01 = np.histogram(Twh01,bins=NumBins,range=[0,Twh01.max()])
CDFWH01 = np.cumsum(HistWH01)/float(NumInit)
widwh01 = BinWH01[0]-BinWH01[1]




test = scipy.stats.anderson_ksamp([Twh001,Twh01])
test2 = scipy.stats.ks_2samp(Twh001, Twh01)
print test
print test2


plt.rc('text', usetex=True)
plt.rc('font', **{'family':'serif','serif':['Helvetica']})


plt.figure()
plt.xlabel('Time to Destabilize')
plt.ylabel('Frequency')
#plt.title('Kepler11 System Scaled by 0.396 and integrated with WHfast')
plt.bar(BinWH001[1:],HistWH001,width=widwh001,label='dt=0.01*T')
#plt.figure()
#plt.xlabel('Time to Destabilize')
#plt.title('Kepler11 System Scaled by 0.396 and integrated with ias15')
plt.bar(BinWH01[1:],HistWH01,width=widwh01,color='none', label='dt=0.1*T')
plt.legend()

plt.figure()
plt.xlabel('Time to Destabilize')
#plt.title('Cumulative Distribution Function - WHfast')
plt.plot(BinWH001[:NumBins],CDFWH001,label='dt=0.01*T')
plt.ylabel('Cumulative Probability')
#plt.xlabel('Time to Destabilize')
#plt.title('Cumulative Distribution Function - ias15')
plt.plot(BinWH01[:NumBins],CDFWH01,label='dt=0.1*T')
plt.legend(loc=4)

plt.show()

'''
