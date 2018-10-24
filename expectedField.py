from field import *
from getMirnov import *
from scipy.constants import mu_0
#SDAS
shotV=42952
shotH=42966
shotP=43066
#Coil signals
vert, times, tbs = getSignal(ch_vert, shotV )
hor, times, tbs = getSignal(ch_hor, shotH )
prim, times, tbs = getSignal(ch_prim, shotP )
#mirnov signals
times, dataV = getMirnovs(shotV,mirnv,True)
times, dataH = getMirnovs(shotH,mirnv,True)
times, dataP = getMirnovs(shotP,mirnv,True)

#computes the flux on each mirnov normalized for a 1 amp current running on coils in Rw,Zw
def getMirnovFlux(Rw_,Zw_,polarity,windings,biotSavart=True):
    #mirnov positions
    radius=9.35 #cm
    angle=345. - 30.*np.arange(12)
    geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
    geometryR=radius*np.cos(np.radians(angle))
    #loop on the mirnovs
    Hr=np.zeros(len(angle))
    Hz=np.zeros(len(angle))
    i=0
    for r,z in zip(geometryR,geometryZ):
        #loop on the PFCs
        for Rw,Zw, sign in zip(Rw_,Zw_,polarity):
            if biotSavart:
                coilHr, coilHz= biotsavart((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2,1.0) #46.
            else:
                coilHr, coilHz= Hcoil((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2) #46.
            Hr[i]+=sign*coilHr
            Hz[i]+=sign*coilHz
        i+=1
    Hr=np.asarray(Hr)
    Hz=np.asarray(Hz)
    Hp=-Hr*np.sin(np.radians(angle))+Hz*np.cos(np.radians(angle))
    return Hp*windings*50*49e-6

V=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5)*340
H=getMirnovFlux([58.,58.],[-7.,7.],[+1.,-1.],4)*260
P=getMirnovFlux([62.,62.],[-13.,13.],[+1.,1.],14)*157

#get te mirnov flat top value with heaviside shots
def flatTops (data):
    return np.asarray([np.mean(np.array(i)[4000:6000]) for i in data])

plt.figure()
plt.plot(np.arange(12)+1,-V*1e8)
plt.plot(np.arange(12)+1,flatTops(dataV)*1e8)

plt.figure()
plt.plot(np.arange(12)+1,H*1e8)
plt.plot(np.arange(12)+1,flatTops(dataH)*1e8)

plt.figure()
plt.plot(np.arange(12)+1,P*1e8)
plt.plot(np.arange(12)+1,flatTops(dataP)*1e8)
