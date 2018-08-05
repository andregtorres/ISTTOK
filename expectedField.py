from field import *
from getMirnov import *
from scipy.constants import mu_0
#SDAS
shotnr=42952;
#Coil signals
vert, times, tbs = getSignal(ch_vert, shotnr )
#mirnov signals
times, data = getMirnovs(shotnr,mirnv,True)

radius=9.35 #cm
angle=345. - 30.*np.arange(12)
geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
geometryR=radius*np.cos(np.radians(angle))

Hr=np.zeros(len(angle))
Hz=np.zeros(len(angle))
for r,z in zip(geometryR,geometryZ):
    for Rw,Zw, sign in zip([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.]):
        coilHr, coilHz= biotsavart(r+46., z, Rw,Zw,1.0) #46.
        Hr+=sign*coilHr
        Hz+=sign*coilHz
Hr=np.asarray(Hr)
Hz=np.asarray(Hz)
Hp=Hr*np.sin(np.radians(angle))-Hz*np.cos(np.radians(angle))
Htot=np.sqrt(Hr**2+Hz**2)
Htot*340*5
plt.plot(np.arange(12)+1,Htot)

plt.plot(np.arange(12)+1,np.array(data)[:,5001])
Hz*340*5
Hp*340

Hp=Hp*50.*49e-6
Vf=5.*np.outer(np.asarray(Hp).T,vert)


for mirnov in Vf[0:3]:
    plt.plot(times,mirnov)
    print np.max(mirnov)
for mirnov in [data[0], data[11], data[5], data[1], data[4]]:
    plt.plot(times,mirnov)
    print np.max(mirnov)
