#Code for acessment of the external coils at equatorial plane
#Andre Torres
#10-01-19
from getMirnov import *
from scipy.optimize import curve_fit
from coilDefinitions import PF0,PF1,PF2, v,ht,hb
from field import flatTops
#%matplotlib qt4

#SDAS INFO
shotN=45078 #44409
#45076 Horizontal -+200A, 400ms
#45078 Vertical -+400A, 400ms
#45083 Primario +140A, 400ms
#45085 Primario -140A, 400ms
# 44835 plasma
'''
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
def correctDrift(V):
    slope=(np.mean(V[1000:1050])-np.mean(V[0:50]))/1000.
    drift=np.mean(V[0:50])+np.arange(len(V))*slope
    return(V-drift)

#CORRECT RAD_U
rad_u=-1*correctDrift(rad_u) #INVERT SIGNAL
rad_b=correctDrift(rad_b)
vert=correctDrift(vert)
'''

'''
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

%matplotlib tk

#Plot 3 signals
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
#plt.xlim([0,1000])
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.tight_layout()
plt.legend()

#COMPUTE PF
v.setSignals(primary,PF_vert,PF_hor)
ht.setSignals(primary,PF_vert,PF_hor)
hb.setSignals(primary,PF_vert,PF_hor)

plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,v.PF0*1e6, label="PF0")
plt.plot(times*1e-3,v.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,v.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,v.PF2*1e6, label="PF2 no gain")
plt.legend()

#horizontal coil top
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,ht.PF0*1e6, label="PF0")
plt.plot(times*1e-3,ht.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,ht.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,ht.PF2*1e6, label="PF2 no gain")
plt.legend()

#horizontal coil bottom
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,hb.PF1*1e6, label="PF0")
plt.plot(times*1e-3,hb.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,ht.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,ht.PF2*1e6, label="PF2 no gain")
plt.legend()

#FIND DECAY TIME
#FIT
def exponential(x, a, b, c):
    return c -a*np.exp(-(x-2010) / b)

def getFit(V,times,full=True):
    #slicing
    if full:
        slice_start=np.where(times==200000)[0][0]
    else:
        slice_start=np.where(times==225000)[0][0]
    slice_end=np.where(times==500000)[0][0]
    times_s=times[slice_start:slice_end]
    V_s=V[slice_start:slice_end]
    guess=[-1.8e-5, 3e4, -3e-5]
    popt, pcov = curve_fit(exponential, times_s, V_s, p0=guess,  maxfev=50000)
    #Calculate R squared
    residuals = V_s - exponential(times_s, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((V_s-np.mean(V_s))**2)
    Rsq = 1.0 - ss_res/ss_tot
    perr=np.sqrt(np.absolute(np.diag(pcov)))
    return popt, perr, Rsq

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
flatTopsPF0=flatTops([v.PF0,ht.PF0,hb.PF0])
flatTopsPF1=flatTops([v.PF1,ht.PF1,hb.PF1])
ratiosPF0=flatTopsSignal/flatTopsPF0
ratiosPF1=flatTopsSignal/flatTopsPF1
ratiosPF0
ratiosPF1


from filters import CSfilter
def cheatFilter(currents, PFset=0):
    #COMPUTE PF
    v.setSignals(*currents)
    ht.setSignals(*currents)
    hb.setSignals(*currents)

    ratiosPF0=np.asarray([0.22671773, 0.12373415, 0.17530242])
    ratiosPF1=np.asarray([0.60805321, 0.9006564 , 1.27601997])
    popt=[[-7.61731222e-03,  3.73277163e+04, -3.64464177e-05], [ 1.05060429e-04,  3.97652858e+04, -9.12632314e-07], [-4.73349226e-05,  6.12193708e+04, -1.29608285e-06]]
    tbs=100.
    if PFset == 0:
        amplified=(np.asarray([v.PF0,hb.PF0,hb.PF0]).T*ratiosPF0).T
    if PFset == 1:
        amplified=(np.asarray([v.PF1,hb.PF1,hb.PF1]).T*ratiosPF1).T

    filtered=[[],[],[]]
    for i in range(3):
        fc1=1./popt[i][1]/2./np.pi
        filtered[i]=CSfilter(amplified[i],fc1, tbs)
    return filtered
filtered0=cheatFilter([primary,PF_vert,PF_hor],PFset=0)
#%matplotlib qt4
plt.figure()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,rad_u*1e6, label="Radial Up")
plt.plot(times*1e-3,rad_b*1e6, label="Radial Down")
plt.plot(times*1e-3,filtered0[0]*1e6)
plt.plot(times*1e-3,filtered0[1]*1e6)
plt.plot(times*1e-3,filtered0[2]*1e6)
plt.legend()
