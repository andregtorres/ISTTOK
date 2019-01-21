#Code for acessment of the external coils at equatorial plane
#Andre Torres
#16-01-19
from getMirnov import *
from scipy.optimize import curve_fit
from coilDefinitions import PF0,PF1,PF2, v,ht,hb
from field import flatTops
%matplotlib qt4

#SDAS INFO
shotN=45210
#45076 Horizontal -+200A, 400ms
#45078 Vertical -+400A, 400ms
#45083 Primario +140A, 400ms
#45085 Primario -140A, 400ms
#44833 plasma 44835 no plasma
#45210 primary slope

ch_rad_u = 'MARTE_NODE_IVO3.DataCollection.Channel_141'
ch_vertical= 'MARTE_NODE_IVO3.DataCollection.Channel_142'
ch_rad_b = 'MARTE_NODE_IVO3.DataCollection.Channel_143'

#reference signals
primary, times_p,tbs=getSignal( 'MARTE_NODE_IVO3.DataCollection.Channel_093', shotN)
PF_vert, times_v,tbs=getSignal( ch_vert, shotN)
PF_hor, times_h,tbs=getSignal( ch_hor, shotN)

#triple sadle
#ADC-Vs factor
vertScale =  1.7102e-4 / 2.0e6  # LSB to Volt * Sampling Period
rad_u, times,tbs=getSignal(ch_rad_u, shotN, vertScale/10.)
rad_b, times,tbs=getSignal(ch_rad_b, shotN, vertScale/10.)
vert, times,tbs=getSignal(ch_vertical, shotN, vertScale/10.)

#correct drift
def correctDrift(V):
    slope=(np.mean(V[1000:1050])-np.mean(V[0:50]))/1000.
    drift=np.mean(V[0:50])+np.arange(len(V))*slope
    return(V-drift)

#def correctDrift(V):
#    drift=np.linspace(0,V[-1], num=len(V))
#    return(V-drift)
#

#CORRECT RAD_U
rad_u=-1*correctDrift(rad_u) #INVERT SIGNAL
rad_b=correctDrift(rad_b)
vert=correctDrift(vert)

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
plt.xlim([0,1000])
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

#vertical coil
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.xlim([0,1000])
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

plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.xlim([0,1000])
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,(vert-v.PF0)*1e6, label="Vertical - PF0")
plt.plot(times*1e-3,(vert-v.PF1)*1e6, label="Vertical - PF1")
plt.legend()

#np.max(vert)/np.max(Hz0)
#1.59
#1.70

#horizontal coil top
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.xlim([0,1000])
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
plt.xlim([0,1000])
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,hb.PF0*1e6, label="PF0")
plt.plot(times*1e-3,hb.PF1*1e6, label="PF1")
#plt.plot(times*1e-3,hb.PF2g*1e6, label="PF2")
#plt.plot(times*1e-3,hb.PF2*1e6, label="PF2 no gain")
plt.legend()


plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.xlim([0,1000])
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Uper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,(rad_u-rad_b)*1e6, label="U-L")
plt.legend()

%matplotlib tk
#FIND DECAY TIME

#FIT
def exponential(x, a, b, c):
    return c -a*np.exp(-(x-200000) / b)

def getFit(V,times,full=True):
    #slicing
    if full:
        slice_start=np.where(times==200000)[0][0]
    else:
        slice_start=np.where(times==225000)[0][0]
    slice_end=np.where(times==500000)[0][0]
    times_s=times[slice_start:slice_end]
    V_s=V[slice_start:slice_end]
    guess=[1.8e-4, 5.1e3, -1.7e-5]
    popt, pcov = curve_fit(exponential, times_s, V_s, p0=guess,  maxfev=50000)
    #Calculate R squared
    residuals = V_s - exponential(times_s, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((V_s-np.mean(V_s))**2)
    Rsq = 1.0 - ss_res/ss_tot
    perr=np.sqrt(np.absolute(np.diag(pcov)))
    return popt, perr, Rsq, times_s

popt=[]
perr=[]
Rsq=[]
for signal in [vert,rad_u,rad_b]:
    popt_, perr_, Rsq_, times_s= getFit(signal,times,True)
    popt.append(popt_)
    perr.append(perr_)
    Rsq.append(Rsq_)
popt
perr
Rsq


plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3, vert*1e6, label="vertical")
#plt.plot(times*1e-3, exponential(times, *guess))
plt.plot(times_s*1e-3, exponential(times_s, *popt[0])*1e6, label=r'Fit, $\tau$=24.36 $\pm$ 0.14 ms')
#plt.plot(times*1e-3,v.PF0*1e6, label="PF0")
plt.legend()

#hacking
flatTopsSignal=flatTops([vert,rad_u,rad_b])
flatTopsPF0=flatTops([v.PF0,ht.PF0,hb.PF0])
flatTopsPF1=flatTops([v.PF1,ht.PF1,hb.PF1])
ratiosPF0=flatTopsSignal/flatTopsPF0
ratiosPF1=flatTopsSignal/flatTopsPF1
ratiosPF0
ratiosPF1

plt.figure()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,v.PF0*ratiosPF0[0]*1e6, label="PF0")
plt.plot(times*1e-3,v.PF1*ratiosPF1[0]*1e6, label="PF1")
plt.legend()


plt.figure()
plt.tight_layout()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,ht.PF0*ratiosPF0[1]*1e6, label="PF0")
plt.plot(times*1e-3,ht.PF1*ratiosPF1[1]*1e6, label="PF1")
plt.legend()


from filters import CSfilter
import scipy.signal as signal

def bandFilter(x,cutoffs, tbs_):
    fs=1./tbs_
    nyq = 0.5 * fs
    filterCutOff=np.asarray(cutoffs)*1./nyq
    b, a = signal.butter(1, filterCutOff, 'bandstop')
    return signal.lfilter(b, a, x)

amplified=(np.asarray([v.PF0,hb.PF0,hb.PF0]).T*ratios).T
fc1=1./2.436e+04/2./np.pi
fc2=1./3.212e+04/2./np.pi
filtered=[[],[],[]]
for i in range(3):
    fc1=1./popt[i][1]/2./np.pi
    filtered[i]=CSfilter(amplified[i],fc1, tbs)

#PF0_filter1=CSfilter(v.PF0*ratios[0],fc1, tbs)
#PF0_filter2=CSfilter(v.PF0*ratios[0],fc2, tbs)
#PF0_filter1=bandFilter(v.PF0*ratios[0],[fc1,fc1*500], tbs)
#PF0_filter2=bandFilter(v.PF0*ratios[0],[fc2,fc2*10], tbs)
popt



def cheatFilter(currents, PFset=0):
    #COMPUTE PF
    v.setSignals(*currents)
    ht.setSignals(*currents)
    hb.setSignals(*currents)

    ratiosPF0=np.asarray([1.12773403, 0.39033158, 1.40338548])
    ratiosPF1=np.asarray([1.46461457, 0.38986613, 1.401712  ])
    popt=[[7.89588047e-05, 2.43622814e+04, 6.55629815e-05],[-2.06196145e-06,  4.38844884e+04, -1.60822894e-06],[-8.12553610e-06,  1.27860286e+04, -6.02436304e-06]]
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

plt.figure()
for a in filtered0:
    plt.plot(times*1e-3,a*1e6)

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
