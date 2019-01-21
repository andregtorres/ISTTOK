#Code for acessment of the external coils at equatorial plane
#Andre Torres
#17-01-19
from getMirnov import *
from scipy.optimize import curve_fit
from coilDefinitions import PF0,PF1,PF2, tripleCoil
from field import flatTops

%matplotlib qt4

#SDAS INFO
shotN=45076 #45206
#45076 Horizontal -+200A, 400ms
#45078 Vertical -+400A, 400ms
#45083 Primario +140A, 400ms
#45085 Primario -140A, 400ms
# 44835 plasma


ch_rad_u = 'MARTE_NODE_IVO3.DataCollection.Channel_141'
ch_vertical= 'MARTE_NODE_IVO3.DataCollection.Channel_142'
ch_rad_b = 'MARTE_NODE_IVO3.DataCollection.Channel_143'

#reference signals
primary, times_p,tbs=getSignal( 'MARTE_NODE_IVO3.DataCollection.Channel_093', shotN)
PF_vert, times_v,tbs=getSignal( ch_vert, shotN)
PF_hor, times_h,tbs=getSignal( ch_hor, shotN)

#triple sadle
#ADC-Vs factor
vertScale =  1.7102e-4 / 2.0e6 # LSB to Volt * Sampling Period
rad_u, times,tbs=getSignal(ch_rad_u, shotN, vertScale/10.)
rad_b, times,tbs=getSignal(ch_rad_b, shotN, vertScale/10.)
vert, times,tbs=getSignal(ch_vertical, shotN, vertScale/10.)

#correct drift
#def correctDrift(V):
#    slope=(np.mean(V[1000:1050])-np.mean(V[0:50]))/1000.
#    drift=np.mean(V[0:50])+np.arange(len(V))*slope
#    return(V-drift)
def correctDrift(V):
    drift=np.linspace(0,V[-1], num=len(V))
    return(V-drift)

#CORRECT RAD_U
rad_u=-1*correctDrift(rad_u) #INVERT SIGNAL
rad_b=correctDrift(rad_b)
vert=correctDrift(vert)
'''
#save files for offline
np.save("dataFiles/FluxCoils/"+str(shotN)+"/times", times)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/primary", primary)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/PF_vert", PF_vert)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/PF_hor", PF_hor)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/density", density)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/mirnovs0", mirnovs[0])
np.save("dataFiles/FluxCoils/"+str(shotN)+"/rad_u", rad_u)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/rad_b", rad_b)
np.save("dataFiles/FluxCoils/"+str(shotN)+"/vert",  vert)
'''
'''
#load files
import numpy as np
import matplotlib.pyplot as plt
times   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/times.npy")
primary =np.load("dataFiles/FluxCoils/"+str(shotN)+"/primary.npy")
PF_vert =np.load("dataFiles/FluxCoils/"+str(shotN)+"/PF_vert.npy")
PF_hor =np.load("dataFiles/FluxCoils/"+str(shotN)+"/PF_hor.npy")
density =np.load("dataFiles/FluxCoils/"+str(shotN)+"/density.npy")
mirnovs0=np.load("dataFiles/FluxCoils/"+str(shotN)+"/mirnovs0.npy")
rad_u   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/rad_u.npy")
rad_b   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/rad_b.npy")
vert    =np.load("dataFiles/FluxCoils/"+str(shotN)+"/vert.npy")
'''
#Plot 3 signals
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.tight_layout()
plt.legend()

#COMPUTE PF
tc=tripleCoil([primary,PF_vert,PF_hor])

plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,tc.v.PF0*1e6, label="PF0")
plt.plot(times*1e-3,tc.v.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,tc.v.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,tc.v.PF2*1e6, label="PF2 no gain")
plt.legend()

#horizontal coil top
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,tc.ht.PF0*1e6, label="PF0")
plt.plot(times*1e-3,tc.ht.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,tc.ht.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,tc.ht.PF2*1e6, label="PF2 no gain")
plt.legend()

#horizontal coil bottom
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,tc.hb.PF0*1e6, label="PF0")
plt.plot(times*1e-3,tc.hb.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,tc.hb.PF2*1e6, label="PF2")
plt.legend()
#in one plot
a=0.7
fig, axs = plt.subplots(3, 1, sharex=True)
for ax in axs:
    ax.grid()
axs[0].set_title("Vertical")
axs[1].set_title("Upper radial")
axs[2].set_title("Lower radial")
axs[0].set_ylabel("Flux [uV.s]")
axs[1].set_ylabel("Flux [uV.s]")
axs[2].set_ylabel("Flux [uV.s]")
axs[2].set_xlabel("Time [ms]")
axs[0].plot(times*1e-3,vert*1e6, label="Signal")
axs[0].plot(times*1e-3,tc.v.PF0*1e6, label="PF0")
axs[0].plot(times*1e-3,tc.v.PF1*1e6, label="PF1")
axs[1].plot(times*1e-3,rad_u*1e6, label="Upper radial")
axs[1].plot(times*1e-3,tc.ht.PF0*1e6, label="PF0")
axs[1].plot(times*1e-3,tc.ht.PF1*1e6, label="PF1")
axs[2].plot(times*1e-3,rad_b*1e6, label="Upper radial")
axs[2].plot(times*1e-3,tc.hb.PF0*1e6, label="PF0 ")
axs[2].plot(times*1e-3,tc.hb.PF1*1e6, label="PF1")
axs[0].legend()



plt.figure()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,(rad_u-rad_b)*1e6, label="U-L")
plt.plot(times*1e-3,(tc.ht.PF0-tc.hb.PF0)*1e6, label="U-L, PF0")
plt.plot(times*1e-3,(tc.ht.PF1-tc.hb.PF1)*1e6, label="U-L, PF1")
plt.legend()


from FluxCoils_vacuum_v import getFit
slice_start=np.where(times==200000)[0][0]
slice_end=np.where(times==500000)[0][0]
times_s=times[slice_start:slice_end]

popt=[]
perr=[]
Rsq=[]
for signal in [vert,rad_u,rad_b]:
    popt_, perr_, Rsq_ = getFit(signal,times,True)
    popt.append(popt_)
    perr.append(perr_)
    Rsq.append(Rsq_)
popt
Rsq

#hacking
flatTopsSignal=flatTops([vert,rad_u,rad_b])
flatTopsPF0=flatTops([tc.v.PF0,tc.ht.PF0,tc.hb.PF0])
flatTopsPF1=flatTops([tc.v.PF1,tc.ht.PF1,tc.hb.PF1])
ratiosPF0=flatTopsSignal/flatTopsPF0
ratiosPF1=flatTopsSignal/flatTopsPF1
ratiosPF0
ratiosPF1
