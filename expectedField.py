from field import *
from getMirnov import *
from scipy.constants import mu_0
#SDAS
shotV=42952
shotH=44330 #44123 175A, 42966 260A, 44330 XA
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
def flatTops (data,from_=4000, to_=6000):
    return np.asarray([np.mean(np.array(i)[from_:to_]) for i in data])

#Horizontal coils position calculation
'''
np.asarray(squareSums)
%matplotlib qt4
plt.figure()
plt.contourf(np.arange(-0.5,0.5,0.1),np.arange(-1.,1.,0.1),np.transpose(np.asarray(squareSums)))
plt.colorbar()
plt.xlabel("dR")
plt.ylabel("dZ")
'''


ftH=flatTops(dataH,2000,4000)
xx1=np.arange(-4,0,0.2)
xx2=np.arange(-7,-2,0.2)
yy1=np.arange(-4,0,0.2)
yy2=np.arange(2,6,0.2)
ii=np.arange(0,1,1)
error_min=[1,0,0,0,0,0]
for di in ii:
    for dx1 in xx1:
        print dx1
        for dx2 in xx2:
            for dy1 in yy1:
                for dy2 in yy2:
                    H=getMirnovFlux([58.+dx1,58.+dx2],[-7+dy1,7.+dy2],[+1.,-1.],4, biotSavart=False)*(180+di)
                    err=np.sqrt(((-ftH-H)**2).sum())
                    if err < error_min[0]:
                        error_min=[err,dx1,dx2,dy1,dy2,di]
error_min
#ErrorMIn (-4,-5.2,-3.2,4.8)

H0=getMirnovFlux([58.,58.],[-7.,7.],[+1.,-1.],4, biotSavart=False)*175 #260
H2=getMirnovFlux([58+error_min[1],58+error_min[2]],[-7+error_min[3],7.+error_min[4]],[+1.,-1.],4, biotSavart=False)*(175+error_min[5])
H3=getMirnovFlux([58-2.52,58-5.3632],[-7-2.756,7.+4.1782],[+1.,-1.],4,biotSavart=False)*(175)


#Variation of I
def varyCurrent(ii,dx1=0,dx2=0,dy1=0,dy2=0):
    var_i=[]
    for di in ii:
        H=getMirnovFlux([58.+dx1,58.+dx2],[-7+dy1,7.+dy2],[+1.,-1.],4, biotSavart=False)*(175+di)
        err=np.sqrt(((ftH+H)**2).sum())
        var_i.append(err)
    return(var_i)
ftH=flatTops(dataH,2000,4000)
ii=np.arange(-50,51,1)
var_i0=varyCurrent(ii)
var_i=varyCurrent(ii,-2.52,-5.36,-2.756,4.1782)
plt.figure()
plt.plot(ii, np.array(var_i0)*1e6, label="original position")
plt.plot(ii, np.array(var_i)*1e6, label="optimized position")
plt.xlabel("dI on active coil [A]")
plt.ylabel("RMS error [uV s]")
plt.title ("Pulse #44330 - Variation of Hfield current")
plt.legend()
plt.grid()

%matplotlib qt4
#Variation of I 2D
def varyCurrent2D(ii,dx1=0,dx2=0,dy1=0,dy2=0):
    var_i2D=[]
    for di1 in ii:
        var_i=[]
        H1=getMirnovFlux([58.+dx1],[-7+dy1],[+1.],4, biotSavart=False)*(175+di1)
        for di2 in ii:
            H2=getMirnovFlux([58.+dx2],[7.+dy2],[-1.],4, biotSavart=False)*(175+di2)
            err=np.sqrt(((-ftH-(H1+H2))**2).sum())
            var_i.append(err)
        var_i2D.append(np.asarray(var_i))
    return(var_i2D)
ftH=flatTops(dataH,2000,4000)
ii=np.arange(-30,72,2)
#var_i0=varyCurrent(ii)
var_i=varyCurrent2D(ii,-2.52,-5.36,-2.756,4.1782)
plt.figure()
plt.contourf(ii,ii,np.transpose(np.asarray(var_i)*1e6))
plt.colorbar()
plt.xlabel("dI1")
plt.ylabel("dI2")


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


H0=getMirnovFlux([58.,58.],[-7.,7.],[+1.,-1.],4, biotSavart=False)*180 #260
H2=getMirnovFlux([58-2.52,58-5.3632],[-7-2.756,7.+4.1782],[+1.,-1.],4,biotSavart=False)*180
H3=getMirnovFlux([58-2.52,58-5.3632],[-7-2.756,7.+4.1782],[+1.,-1.],4,biotSavart=False)*(175+30)

H11=getMirnovFlux([58.+dx1],[-7.+dy1],[1.],4, biotSavart=False)*(175+32)
H12=getMirnovFlux([58.+dx2],[7.+dy2],[-1.],4, biotSavart=False)*(175+17)
plt.figure()
plt.plot(np.arange(12)+1,-ftH*1e6, label="Measured")
plt.plot(np.arange(12)+1,H2*1e6,label="Optimized, 180A")
plt.plot(np.arange(12)+1,(H11+H12)*1e6,label="Optimized, 208,192A")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44330 - Mirnov flux with optimized coil position")
plt.legend()




plt.figure()
plt.plot(np.arange(12)+1,-ftH*1e6, label="Measured")
#plt.plot(np.arange(12)+1,H0*1e6,label="Original, 180A")
plt.plot(np.arange(12)+1,H2*1e6,label="Optimized, 180A")
plt.plot(np.arange(12)+1,H3*1e6,label="Optimized, 205A")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44330 - Mirnov flux with optimized coil position")
plt.legend()


'''
plt.figure()
plt.plot(np.arange(12)+1,P*1e8)
plt.plot(np.arange(12)+1,flatTops(dataP)*1e8)
'''
