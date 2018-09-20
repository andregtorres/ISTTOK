#Andre Torres - 29/07/2018
#applies the filter to a real plasma pulse

from getMirnov import *
from filters import CSfilter
from scipy.signal import savgol_filter

Vgains=np.array([1.73E-08,3.83E-08,1.62E-08,-2.44E-08,-2.88E-08,-1.67E-08,-1.86E-08,-3.78E-08,-2.28E-08,1.84E-08,2.75E-08,1.45E-08])
Pgains=np.array([7.435532510975972786e-08,7.810303078411286263e-08,2.560861005186963348e-08,-3.252130268641640995e-08,-7.464898517923825386e-08,-8.264170300351341432e-08,-8.223951515040492875e-08,-6.606898430625495534e-08,-3.312167971342630853e-08,2.401873722827468487e-08,8.358529401801281549e-08,7.283223519795826231e-08])
Hgains=np.array([4.036905791462669659e-08,7.544069073198127677e-08,3.278532096899874485e-08,2.192334653621645880e-08,0.0,0.0,0.0,0.0,-2.947408155500183389e-08,-4.780829102874075150e-08,-3.287957774601663319e-08,-2.172399946407347938e-08])
Vfc=np.array([6.36		,7.08		,9.23		,7.08		,6.47		,4.06		,4.35		,6.88		,5.48		,5.73		,5.63		,3.72])
Pfc=np.array([8.313019851691642259,6.106860538591612375,4.887008770245405920,3.694572564114741997,7.229626404862063538,8.588720809064279038,6.875714847099883009,5.917367134596915434,6.191243890425031537,2.245538136116884687,5.651611285594944967,5.773617831023106994])
Hfc=np.array([27.07,17.02,9.66,27.57,0.0,0.0,0.0,0.0,26.63,13.42,18.52,29.01])


#filter currents
def getCorrections(vert_,prim_,hor_):
    #create the filters
    filteredV=[]
    filteredP=[]
    filteredH=[]
    for vfc, pfc, hfc in zip(Vfc,Pfc,Hfc):
        filteredV.append(CSfilter(vert_,vfc*1e-6, tbs))
        filteredP.append(CSfilter(prim_,pfc*1e-6, tbs))
        #account for missing fits
        if hfc != 0.0:
            filteredH.append(CSfilter(hor_,hfc*1e-6, tbs))
        else:
            filteredH.append(hor_)
    filteredV=np.asarray(filteredV)
    filteredP=np.asarray(filteredP)
    filteredH=np.asarray(filteredH)
    return filteredV, filteredP, filteredH

def getCentroid(data_,correctionV_=np.asarray([[]]),correctionH_=np.asarray([[]]),correctionP_=np.asarray([[]])):
    radius=9.35 #cm
    if correctionV_.shape != data_.shape: correctionV_=np.zeros(data_.shape)
    if correctionH_.shape != data_.shape: correctionH_=np.zeros(data_.shape)
    if correctionP_.shape != data_.shape: correctionP_=np.zeros(data_.shape)
    corrected=np.asarray(data_-correctionV_-correctionH_-correctionP_)
    Hsum_corr=np.sum(corrected, axis=0)
    angle=345. - 30.*np.arange(12)
    geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
    geometryR=radius*np.cos(np.radians(angle))
    z=np.dot(geometryZ,corrected)/Hsum_corr
    R=np.dot(geometryR,corrected)/Hsum_corr
    return R,z


shotNr=43740

#mirnov signals
times, data = getMirnovs(shotNr,mirnv_corr,False)
data=np.asarray(data) #from list to ndarray

#Coil signals
vert, times, tbs = getSignal(ch_vert, shotNr )
hor, times, tbs = getSignal(ch_hor, shotNr )
prim, times, tbs = getSignal(ch_prim, shotNr )

#refernce signals
density, times, tbs= getSignal("MARTE_NODE_IVO3.DataCollection.Channel_088",shotNr)
density_filtered=np.ndarray.flatten(savgol_filter(density,29,5))

#Filter currents
filteredV, filteredP, filteredH =getCorrections(vert,prim,hor)
#compute the corrections
correctionV=(Vgains*filteredV.T).T
correctionP=(Pgains*filteredP.T).T
correctionH=(Hgains*filteredH.T).T

plt.figure()
plt.plot(times*1e-3, data[7]*1e6, label="Original")
plt.plot(times*1e-3, (data[7]-correctionV[7])*1e6, label="V PF corrected")
plt.plot(times*1e-3, (data[7]-correctionP[7])*1e6, label="P PF corrected")
plt.plot(times*1e-3, (data[7]-correctionH[7])*1e6, label="H PF corrected")
plt.xlabel("Time (ms)")
plt.ylabel("uV.s")
plt.title("MIRNOV 8 signal")
plt.legend()

plt.figure()
plt.plot(times*1e-3, density)
plt.plot(times*1e-3, density_filtered)
plt.title("Interferometer Density")
plt.xlim(60, 80)
plt.ylim(3000, 5000)
plt.xlabel("Time (ms)")
plt.grid()

plt.figure()
plt.plot(times*1e-3, vert, label= "vertical")
plt.plot(times*1e-3, filteredV[2], label= "v. filter w mirnov 3 params")
plt.plot(times*1e-3, hor, label= "horizontal")
plt.plot(times*1e-3, filteredH[2], label= "h. filter w mirnov 3 params")
plt.plot(times*1e-3, prim, label= "primary")
plt.plot(times*1e-3, filteredP[2], label= "p. filter w mirnov 3 params")
plt.xlim(60, 80)
plt.legend()
plt.xlabel("Time (ms)")
plt.ylabel("A")
plt.title("PF COILS")

#COMPUTE CENTROID
R,z=getCentroid(data)
R_corr,z_corr=getCentroid(data,correctionV,correctionH)

plt.figure()
plt.plot(times*1e-3, z, label="z")
plt.plot(times*1e-3, R, label="R")
plt.plot(times*1e-3, z_corr, label="z corr.")
plt.plot(times*1e-3, R_corr, label="R corr.")
plt.xlim(60, 80)
plt.ylim(-2, 2)
plt.xlabel("Time (ms)")
plt.ylabel("Displacement (cm)")
plt.title("Current centroid position")
plt.grid()
leg = plt.legend()

RLimiter=8.5
d=np.sqrt(R**2+z**2)
plasmaRadius=RLimiter-d
d_corr=np.sqrt(R_corr**2+z_corr**2)
plasmaRadius_corr=RLimiter-d_corr
chord_corr=2.*np.sqrt(plasmaRadius_corr**2-R_corr**2)
chord=2.*np.sqrt(plasmaRadius**2-R**2)

plt.figure()
plt.plot(times*1e-3, chord, label="Original")
plt.plot(times*1e-3, chord_corr, label="Corrected")
plt.xlim(60, 80)
plt.ylim(0, 17)
plt.xlabel("Time (ms)")
plt.ylabel("Length (cm)")
plt.title("Chord Length at R=0")
plt.grid()
leg = plt.legend()
plt.show()

#slicing
slice_start=np.where(times==60000)[0][0]
slice_end=np.where(times==80000)[0][0]
decimation=10
timesS=times[slice_start:slice_end:decimation]
RS=R[slice_start:slice_end:decimation]
zS=z[slice_start:slice_end:decimation]
R_corrS=R_corr[slice_start:slice_end:decimation]
z_corrS=z_corr[slice_start:slice_end:decimation]
plasmaRadius_corrS=plasmaRadius_corr[slice_start:slice_end:decimation]
plasmaRadiusS=plasmaRadius[slice_start:slice_end:decimation]

#Thanks to:
#https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/
from matplotlib.animation import FuncAnimation
fig, ax = plt.subplots(figsize=(4, 4))
limiter = plt.Circle((0, 0), RLimiter, color='k', fill=False)
plasma_corr = plt.Circle((R_corrS[0], z_corrS[0]), plasmaRadius_corrS[0], color='r', fill=False)
ax.add_artist(limiter)
ax.add_artist(plasma_corr)
ax.set(xlim=(-9,9), ylim=(-9,9))
plt.xlabel("R (cm)")
plt.ylabel("z (cm)")
ax.set_aspect('equal')
scat1 = ax.scatter(np.array(RS[0]), zS[0])
scat2 = ax.scatter(R_corrS[0], z_corrS[0])
def animate(i):
    # Must pass scat.set_offsets an N x 2 array
    global plasma_corr
    ax.set_title('t=' + str(timesS[i]*1e-3)+" ms" )
    scat1.set_offsets(np.c_[RS[i], zS[i]])
    scat2.set_offsets(np.c_[R_corrS[i], z_corrS[i]])
    plasma_corr.remove()
    plasma_corr = plt.Circle((R_corrS[i], z_corrS[i]), plasmaRadius_corrS[i], color='r', fill=False)
    ax.add_artist(plasma_corr)
anim = FuncAnimation(fig, animate, interval=100, frames=len(timesS)-1, repeat=True ) #fargs=[plasma_corr]
anim.save('centroid_noprim.gif', writer='imagemagick')
