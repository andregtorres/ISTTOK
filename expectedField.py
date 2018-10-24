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

ftV=flatTops(dataV)
squareSums=[]
for dx in np.arange(-0.5,0.5,0.1):
    err1=[]
    for dy in np.arange(-1.,1.,0.1):
        V=getMirnovFlux([58.+dx,58.+dx,35.+dx,35.+dx],[-7.+dy,7.+dy,-7.+dy,7.+dy],[-1.,-1.,1.,1.],5)*340
        err=((ftV-(-V))**2)
        err1.append(err.sum())
    squareSums.append(np.asarray(err1))

np.asarray(squareSums)
%matplotlib qt4
plt.figure()
plt.contourf(np.arange(-0.5,0.5,0.1),np.arange(-1.,1.,0.1),np.transpose(np.asarray(squareSums)))
plt.colorbar()
plt.xlabel("dR")
plt.ylabel("dZ")

V2=getMirnovFlux([58.-0.25,58.-0.25,35.-0.25,35.-0.25],[-7.,7.+1,-7.,7.+1],[-1.,-1.,1.,1.],5)*340
H2=getMirnovFlux([58-3.5,58-5],[-7-1,7.+5],[+1.,-1.],4)*260
H0=getMirnovFlux([58.,58.],[-7.,7.],[+1.,-1.],4, biotSavart=False)*260


ftH=flatTops(dataH)
xx1=np.arange(-4,4,0.5)
xx2=np.arange(-4,4,0.5)
yy1=np.arange(-4,4,0.5)
yy2=np.arange(-4,4,0.5)
ii=np.arange(-10,10,2)
error_min=[1,0,0,0,0,0]
for di in ii:
    print di
    for dx1 in xx1:
        for dx2 in xx2:
            for dy1 in yy1:
                for dy2 in yy2:
                    H=getMirnovFlux([58.+dx1,58.+dx2],[-7+dy1,7.+dy2],[+1.,-1.],4, biotSavart=False)*(260+di)
                    err=np.sqrt(((ftH-H)**2).sum())
                    if err < error_min[0]:
                        error_min=[err,dx1,dx2,dy1,dy2,di]
error_min
H2=getMirnovFlux([58+error_min[1],58+error_min[2]],[-7+error_min[3],7.+error_min[4]],[+1.,-1.],4)*(260+error_min[5])

'''
squareSums=[]
for dx1 in xx1:
    err1=[]
    H=getMirnovFlux([58.+dx1,58.+dx2],[-7-1,7.+4],[+1.,-1.],4, biotSavart=False)*260
    err=((ftH-H)**2)
    err1.append(err.sum())
    for dx2 in xx2:
    squareSums.append(np.asarray(err1))

plt.figure()
plt.contourf(xx1,xx2,np.transpose(np.log(np.asarray(squareSums))))
plt.colorbar()
plt.xlabel("dx1")
plt.ylabel("dx2")
'''


'''
plt.figure()
plt.plot(np.arange(12)+1,-V*1e8)
plt.plot(np.arange(12)+1,-V2*1e8)
plt.plot(np.arange(12)+1,flatTops(dataV)*1e8)
'''
plt.figure()
plt.plot(np.arange(12)+1,H0*1e8)
plt.plot(np.arange(12)+1,flatTops(dataH)*1e8)
plt.plot(np.arange(12)+1,H2*1e8)

'''
plt.figure()
plt.plot(np.arange(12)+1,P*1e8)
plt.plot(np.arange(12)+1,flatTops(dataP)*1e8)
'''
