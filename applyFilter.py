#Andre Torres - 29/07/2018
#applies the filter to a real plasma pulse

from getMirnov import *
from filters import CSfilter

shotNr=43740

#mirnov signals
times, data = getMirnovs(shotNr,mirnv,True)

#Coil signals
vert, times, tbs = getSignal(ch_vert, shotNr )
hor, times, tbs = getSignal(ch_hor, shotNr )
prim, times, tbs = getSignal(ch_prim, shotNr )


plt.figure()
plt.plot(times, vert, label= "vertical")
plt.plot(times, hor, label= "horizontal")
plt.plot(times, prim, label= "primary")
plt.xlim(60000, 80000)
plt.legend()


#filter
coilNr=0
filteredV=[]
correctionV=[]
Vgains=np.array([1.73E-08,3.83E-08,1.62E-08,-2.44E-08,-2.88E-08,-1.67E-08,-1.86E-08,-3.78E-08,-2.28E-08,1.84E-08,2.75E-08,1.45E-08])
Vfc=np.array([6.36		,7.08		,9.23		,7.08		,6.47		,4.06		,4.35		,6.88		,5.48		,5.73		,5.63		,3.72])
for mirnov in data:
    coilNr+=1
    filteredV.append(CSfilter(vert,Vfc[coilNr-1]*1e-6, tbs))
    correctionV.append(Vgains[coilNr-1]*filteredV[coilNr-1])
plt.figure()
plt.plot(times, data[7])
plt.plot(times, data[7]-correctionV[7])
plt.xlim(60000, 80000)

#centroid
radius=9.35 #cm

data=np.asarray(data) #from list to ndarray
corrected=np.asarray(data-correctionV)
Hsum=np.sum(data, axis=0)
Hsum_corr=np.sum(corrected, axis=0)
angle=345. - 30.*np.arange(12)
geometryZ=radius*np.sin(np.radians(angle)) #positions
geometryR=radius*np.cos(np.radians(angle))
z=np.dot(geometryZ,data)/Hsum
R=np.dot(geometryR,data)/Hsum
zV=np.dot(geometryZ,corrected)/Hsum_corr
RV=np.dot(geometryR,corrected)/Hsum_corr

plt.figure()
plt.plot(times, z, label="z")
plt.plot(times, R, label="R")
plt.plot(times, zV, label="zV")
plt.plot(times, RV, label="RV")
plt.xlim(60000, 80000)
plt.ylim(-2, 2)
leg = plt.legend()
plt.show()
