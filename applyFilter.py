#Andre Torres - 29/07/2018
#applies the filter to a real plasma pulse

from getMirnov import *
from filters import CSfilter

shotNr=43740

#mirnov signals
times, data = getMirnovs(shotNr,mirnv,True)
data=np.asarray(data) #from list to ndarray

#Coil signals
vert, times, tbs = getSignal(ch_vert, shotNr )
hor, times, tbs = getSignal(ch_hor, shotNr )
prim, times, tbs = getSignal(ch_prim, shotNr )

#filter
Vgains=np.array([1.73E-08,3.83E-08,1.62E-08,-2.44E-08,-2.88E-08,-1.67E-08,-1.86E-08,-3.78E-08,-2.28E-08,1.84E-08,2.75E-08,1.45E-08])
Vfc=np.array([6.36		,7.08		,9.23		,7.08		,6.47		,4.06		,4.35		,6.88		,5.48		,5.73		,5.63		,3.72])
filteredV=[]
for fc in Vfc:
    filteredV.append(CSfilter(vert,fc*1e-6, tbs))
filteredV=np.asarray(filteredV)
correctionV=(Vgains*filteredV.T).T

plt.figure()
plt.plot(times*1e-3, data[7]*1e6, label="Original")
plt.plot(times*1e-3, (data[7]-correctionV[7])*1e6, label="PF corrected")
plt.xlim(60, 80)
plt.xlabel("Time (ms)")
plt.ylabel("uV.s")
plt.title("MIRNOV 8 signal")
plt.legend()

plt.figure()
plt.plot(times*1e-3, vert, label= "vertical")
plt.plot(times*1e-3, filteredV[7], label= "v. filter w mirnov 8 params")
plt.plot(times*1e-3, hor, label= "horizontal")
plt.plot(times*1e-3, prim, label= "primary")
plt.xlim(60, 80)
plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("A")
plt.title("PF COILS")


#centroid
radius=9.35 #cm
corrected=np.asarray(data-correctionV)
Hsum=np.sum(data, axis=0)
Hsum_corr=np.sum(corrected, axis=0)
angle=345. - 30.*np.arange(12)
geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
geometryR=radius*np.cos(np.radians(angle))
z=np.dot(geometryZ,data)/Hsum
R=np.dot(geometryR,data)/Hsum
zV=np.dot(geometryZ,corrected)/Hsum_corr
RV=np.dot(geometryR,corrected)/Hsum_corr

plt.figure()
plt.plot(times*1e-3, z, label="z")
plt.plot(times*1e-3, R, label="R")
plt.plot(times*1e-3, zV, label="zV")
plt.plot(times*1e-3, RV, label="RV")
plt.xlim(60, 80)
plt.ylim(-2, 2)
plt.xlabel("Time (ms)")
plt.ylabel("Displacement (cm)")
leg = plt.legend()
plt.show()

#slicing
slice_start=np.where(times==60000)[0][0]
slice_end=np.where(times==80000)[0][0]
decimation=10
timesS=times[slice_start:slice_end:decimation]
RS=R[slice_start:slice_end:decimation]
zS=z[slice_start:slice_end:decimation]
RVS=RV[slice_start:slice_end:decimation]
zVS=zV[slice_start:slice_end:decimation]

#Thanks to:
#https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
from matplotlib.animation import FuncAnimation
RLimiter=8.5
fig, ax = plt.subplots(figsize=(4, 4))
limiter = plt.Circle((0, 0), RLimiter, color='b', fill=False)
ax.add_artist(limiter)
ax.set(xlim=(-9,9), ylim=(-9,9))
plt.xlabel("R (cm)")
plt.ylabel("z (cm)")
ax.set_aspect('equal')
scat1 = ax.scatter(np.array(RS[0]), zS[0])
scat2 = ax.scatter(RVS[0], zVS[0])
def animate(i):
    # Must pass scat.set_offsets an N x 2 array
    ax.set_title('t=' + str(timesS[i]*1e-3)+" ms" )
    scat1.set_offsets(np.c_[zS[i], RS[i]])
    scat2.set_offsets(np.c_[zVS[i], RVS[i]])

anim = FuncAnimation(fig, animate, interval=100, frames=len(timesS)-1, repeat=True)
anim.save('centroid.gif', writer='imagemagick')
