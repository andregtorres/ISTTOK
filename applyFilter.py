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
correctedV=[]
Vgains=np.array([1.73E-08,3.83E-08,1.62E-08,-2.44E-08,-2.88E-08,-1.67E-08,-1.86E-08,-3.78E-08,-2.28E-08,1.84E-08,2.75E-08,1.45E-08])
Vfc=np.array([6.36		,7.08		,9.23		,7.08		,6.47		,4.06		,4.35		,6.88		,5.48		,5.73		,5.63		,3.72])
for mirnov in data:
    coilNr+=1
    filteredV.append(CSfilter(vert,Vfc[coilNr-1]*1e-6, tbs))
    correctedV.append(mirnov-Vgains[coilNr-1]*filteredV[coilNr-1])
plt.figure()
plt.plot(times, data[7])
plt.plot(times, correctedV[7])
plt.xlim(60000, 80000)

print (correctedV[0][2000]-data[0][2000])/data[0][2000]

#centroid
coilNr=0
Hsum=np.zeros(len(data[0]))
z=np.zeros(len(data[0]))
R=np.zeros(len(data[0]))
HsumV=np.zeros(len(data[0]))
zV=np.zeros(len(data[0]))
RV=np.zeros(len(data[0]))
for coil in data:
    coilNr+=1
    angle=345. - 30.*(coilNr-1.)
    print "coil: " + str(coilNr) + " angle: " + str(angle) + " pos: (" + str(9.35*(np.cos(np.radians(angle)))) + str(9.35*(np.sin(np.radians(angle))))
    Hsum += coil
    z += coil*9.35*(np.sin(np.radians(angle)))
    R += coil*9.35*(np.cos(np.radians(angle)))
    HsumV += correctedV[coilNr-1]
    zV += correctedV[coilNr-1]*9.35*(np.sin(np.radians(angle)))
    RV += correctedV[coilNr-1]*9.35*(np.cos(np.radians(angle)))
zt=z/Hsum
Rt=R/Hsum
ztV=zV/HsumV
RtV=RV/HsumV

print (ztV[2000]-zt[2000])/zt[2000]

plt.figure()
plt.plot(times, zt, label="z")
plt.plot(times, Rt, label="R")
plt.plot(times, ztV, label="zV")
plt.plot(times, RtV, label="RV")
plt.xlim(60000, 80000)
plt.ylim(-2, 2)
leg = plt.legend()
#plt.plot(times, data[2]);
plt.show()




plt.show()
