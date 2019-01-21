#removes the PFcoils contribuitions from the external coil setSignals
# Andre Torres
# 19.01.19
from coilDefinitions import PF0, PF1, PF2, tripleCoil
from getMirnov import *

#SDAS INFO
shotN=44835 #44833
# 44835 no plasma; 44833 plasma
def getSDAS(shotN):
    ch_rad_u = 'MARTE_NODE_IVO3.DataCollection.Channel_141'
    ch_vertical= 'MARTE_NODE_IVO3.DataCollection.Channel_142'
    ch_rad_b = 'MARTE_NODE_IVO3.DataCollection.Channel_143'

    #reference signals
    primary, times_p,tbs=getSignal( ch_prim, shotN)
    PF_vert, times_v,tbs=getSignal( ch_vert, shotN)
    PF_hor, times_h,tbs=getSignal( ch_hor, shotN)

    #triple coil
    #ADC-Vs factor
    vertScale =  1.7102e-4 / 2.0e6 /10. # LSB to Volt * Sampling Period
    rad_u, times,tbs=getSignal(ch_rad_u, shotN, vertScale)
    rad_b, times,tbs=getSignal(ch_rad_b, shotN, vertScale)
    vert, times,tbs=getSignal(ch_vertical, shotN, vertScale)

    def correctDrift(V):
        if shotN==44833:
            drift=np.linspace(np.mean(V[0:3]),np.mean(V[10900:10903]), num=len(V))
        else:
            slope=(np.mean(V[-4:-1]))/(len(V))
            #slope=0.
            drift=np.arange(len(V))*slope
        return(V-drift)

    #CORRECT RAD_U
    rad_u=-1*correctDrift(rad_u) #INVERT SIGNAL
    rad_b=correctDrift(rad_b)
    vert=correctDrift(vert)


    return times, primary, PF_vert, PF_hor, rad_u, rad_b, vert

def saveSignals(shotN, times, primary, PF_vert, PF_hor, rad_u, rad_b, vert):
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/times", times)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/primary", primary)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/PF_vert", PF_vert)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/PF_hor", PF_hor)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/rad_u", rad_u)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/rad_b", rad_b)
    np.save("dataFiles/FluxCoils/"+str(shotN)+"/vert",  vert)

def loadSignals(shotN):
    times   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/times.npy")
    primary =np.load("dataFiles/FluxCoils/"+str(shotN)+"/primary.npy")
    PF_vert =np.load("dataFiles/FluxCoils/"+str(shotN)+"/PF_vert.npy")
    PF_hor =np.load("dataFiles/FluxCoils/"+str(shotN)+"/PF_hor.npy")
    rad_u   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/rad_u.npy")
    rad_b   =np.load("dataFiles/FluxCoils/"+str(shotN)+"/rad_b.npy")
    vert    =np.load("dataFiles/FluxCoils/"+str(shotN)+"/vert.npy")
    return times, primary, PF_vert, PF_hor, rad_u, rad_b, vert

#
times, primary, PF_vert, PF_hor, rad_u, rad_b, vert = getSDAS(shotN)
#saveSignals(shotN, times, primary, PF_vert, PF_hor, rad_u, rad_b, vert)
#times, primary, PF_vert, PF_hor, rad_u, rad_b, vert=loadSignals(shotN)

#%matplotlib qt4
#Plot 3 signals
plt.figure(figsize=(8, 6), dpi=100)
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.grid()
#plt.xlim([0,1000])
#plt.ylim([-100,110])
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.tight_layout()
plt.legend()

plt.plot(primary)
plt.plot(PF_vert)
plt.plot(PF_hor)
#COMPUTE PF
tc=tripleCoil([primary,PF_vert,PF_hor],True, True)
tc0=tripleCoil([primary,PF_vert,PF_hor],False, False)
tc1=tripleCoil([primary,PF_vert,PF_hor],True, False)

#With gains
plt.figure(figsize=(8, 6), dpi=100)
plt.grid()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.plot(times*1e-3,tc1.v.PF1*1e6, label="PF1 with gain")
plt.tight_layout()
plt.legend()

plt.figure()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,tc1.ht.PF0*1e6, label="PF0")
plt.plot(times*1e-3,tc1.ht.PF1*1e6, label="PF1")
plt.tight_layout()
plt.legend()

plt.figure()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,tc1.hb.PF0*1e6, label="PF0")
plt.plot(times*1e-3,tc1.hb.PF1*1e6, label="PF1")
plt.tight_layout()
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
axs[0].plot(times*1e-3,tc1.v.PF1*1e6, alpha=a, label="PF1 with gain")
axs[1].plot(times*1e-3,rad_u*1e6, label="Upper radial")
axs[1].plot(times*1e-3,tc1.ht.PF1*1e6, alpha=a,label="PF1 with gain")
axs[2].plot(times*1e-3,rad_b*1e6, label="Upper radial")
axs[2].plot(times*1e-3,tc1.hb.PF1*1e6,alpha=a, label="PF1 with gain")
axs[0].legend()


#With NO gains
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
axs[0].plot(times*1e-3,tc0.v.PF0*1e6, alpha=a, label="PF0")
axs[0].plot(times*1e-3,tc0.v.PF1*1e6, alpha=a, label="PF1")
axs[1].plot(times*1e-3,rad_u*1e6, label="Upper radial")
axs[1].plot(times*1e-3,tc0.ht.PF0*1e6, alpha=a,label="PF0")
axs[1].plot(times*1e-3,tc0.ht.PF1*1e6, alpha=a,label="PF1")
axs[2].plot(times*1e-3,rad_b*1e6, label="Upper radial")
axs[2].plot(times*1e-3,tc0.hb.PF0*1e6,alpha=a, label="PF0")
axs[2].plot(times*1e-3,tc0.hb.PF1*1e6,alpha=a, label="PF1")
axs[0].legend()


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
axs[0].plot(times*1e-3,tc1.v.PF1*1e6, alpha=a, label="PF1 with gain")
axs[0].plot(times*1e-3,tc0.v.PF1*1e6, alpha=a, label="PF1")
axs[1].plot(times*1e-3,rad_u*1e6, label="Upper radial")
axs[1].plot(times*1e-3,tc1.ht.PF1*1e6, alpha=a,label="PF0")
axs[1].plot(times*1e-3,tc0.ht.PF1*1e6, alpha=a,label="PF1")
axs[2].plot(times*1e-3,rad_b*1e6, label="Upper radial")
axs[2].plot(times*1e-3,tc1.hb.PF1*1e6,alpha=a, label="PF0")
axs[2].plot(times*1e-3,tc0.hb.PF1*1e6,alpha=a, label="PF1")
axs[0].legend()


#Diferences
from scipy.signal import savgol_filter
a=0.3
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

axs[0].plot(times*1e-3,(vert-tc0.v.PF1)*1e6, color="b", alpha=a)
axs[0].plot(times*1e-3,savgol_filter((vert-tc0.v.PF1)*1e6,101,2), color="b", label="Signal - PF1")
axs[0].plot(times*1e-3,(vert-tc1.v.PF1)*1e6, color="orange", alpha=a)
axs[0].plot(times*1e-3,savgol_filter((vert-tc1.v.PF1)*1e6,101,2), color="orange", label="Signal - PF1 w/ gain")

axs[1].plot(times*1e-3,(rad_u-tc0.ht.PF1)*1e6,color="b", alpha=a)
axs[1].plot(times*1e-3,savgol_filter((rad_u-tc0.ht.PF1)*1e6,101,2), color="b",label="Upper radial - PF1")
axs[1].plot(times*1e-3,(rad_u-tc1.ht.PF1)*1e6,color="orange", alpha=a)
axs[1].plot(times*1e-3,savgol_filter((rad_u-tc1.ht.PF1)*1e6,101,2), color="orange",label="Upper radial - PF1 w/ gain")

axs[2].plot(times*1e-3,(rad_b-tc0.hb.PF1)*1e6, color="b", alpha=a)
axs[2].plot(times*1e-3,savgol_filter((rad_b-tc0.hb.PF1)*1e6,101,2), color="b",label="Lower radial - PF1")
axs[2].plot(times*1e-3,(rad_b-tc1.hb.PF1)*1e6,color="orange", alpha=a)
axs[2].plot(times*1e-3,savgol_filter((rad_b-tc1.hb.PF1)*1e6,101,2), color="orange",label="Lower radial - PF1 w/ gain")

axs[0].legend()


#Plasma SHOT
shotN=44833
times, primary, PF_vert, PF_hor, rad_u, rad_b, vert = getSDAS(shotN)

#Plot 3 signals
plt.figure()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.tight_layout()
plt.legend()


#COMPUTE PF
tc0=tripleCoil([primary,PF_vert,PF_hor],False, False)
tc1=tripleCoil([primary,PF_vert,PF_hor],True, False)

#ip, times_ip,tbs=getSignal( "POST.PROCESSED.DENSITY", shotN)
ip, times_ip,tbs=getSignal( "MARTE_NODE_IVO3.DataCollection.Channel_088", shotN)


#Diferences
a=0.3
fig, axs = plt.subplots(4, 1, sharex=True)
for ax in axs:
    ax.grid()
axs[0].set_title("Vertical")
axs[1].set_title("Upper radial")
axs[2].set_title("Lower radial")
axs[0].set_ylabel("Flux [uV.s]")
axs[1].set_ylabel("Flux [uV.s]")
axs[2].set_ylabel("Flux [uV.s]")
axs[3].set_ylabel("Iplasma [A]")

axs[3].set_xlabel("Time [ms]")

axs[0].plot(times*1e-3,(vert-tc0.v.PF1)*1e6, color="b", alpha=a)
axs[0].plot(times*1e-3,savgol_filter((vert-tc0.v.PF1)*1e6,51,6), color="b", label="Signal - PF1")
axs[0].plot(times*1e-3,(vert-tc1.v.PF1)*1e6, color="orange", alpha=a)
axs[0].plot(times*1e-3,savgol_filter((vert-tc1.v.PF1)*1e6,51,6), color="orange", label="Signal - PF1 w/ gain")
axs[0].plot(times*1e-3,(vert)*1e6, color="g", alpha=a)
axs[0].plot(times*1e-3,savgol_filter((vert)*1e6,51,6), color="g", label="Signal")

axs[1].plot(times*1e-3,(rad_u-tc0.ht.PF1)*1e6,color="b", alpha=a)
axs[1].plot(times*1e-3,savgol_filter((rad_u-tc0.ht.PF1)*1e6,51,6), color="b",label="Upper radial - PF1")
axs[1].plot(times*1e-3,(rad_u-tc1.ht.PF1)*1e6,color="orange", alpha=a)
axs[1].plot(times*1e-3,savgol_filter((rad_u-tc1.ht.PF1)*1e6,51,6), color="orange",label="Upper radial - PF1 w/ gain")
axs[1].plot(times*1e-3,(rad_u)*1e6,color="g", alpha=a)
axs[1].plot(times*1e-3,savgol_filter((rad_u)*1e6,51,6), color="g",label="Upper radial - PF1 w/ gain")

axs[2].plot(times*1e-3,(rad_b-tc0.hb.PF1)*1e6, color="b", alpha=a)
axs[2].plot(times*1e-3,savgol_filter((rad_b-tc0.hb.PF1)*1e6,51,6), color="b",label="Lower radial - PF1")
axs[2].plot(times*1e-3,(rad_b-tc1.hb.PF1)*1e6,color="orange", alpha=a)
axs[2].plot(times*1e-3,savgol_filter((rad_b-tc1.hb.PF1)*1e6,51,6), color="orange",label="Lower radial - PF1 w/ gain")
axs[2].plot(times*1e-3,(rad_b)*1e6,color="g", alpha=a)
axs[2].plot(times*1e-3,savgol_filter((rad_b)*1e6,51,6), color="g",label="Lower radial - PF1 w/ gain")

axs[3].plot(times_ip*1e-3,(ip),color="k")

axs[0].legend()
axs[0].set_xlim([0,1000])
axs[2].set_ylim([-50,50])

#RADIAL COMPARISSOM
plt.figure()
plt.xlim([0,1000])
#plt.ylim([-60,60])
plt.grid()
plt.title("Pulse #"+str(shotN)+ " Radial flux difference")
plt.ylabel("Flux [uV.s]")
plt.xlabel("Time [ms]")
#plt.plot(times*1e-3,(rad_u-rad_b)*1e6, alpha=0.5, label="U - L")
plt.plot(times*1e-3,(tc0.ht.PF1-tc0.hb.PF1)*1e6, label="PF1")
plt.plot(times*1e-3,(tc1.ht.PF1-tc1.hb.PF1)*1e6, label="PF1 with gain")
plt.tight_layout()
plt.legend()

#matplotlib qt4

plt.figure()
plt.plot(primary)

ip, times_ip,tbs=getSignal( "MARTE_NODE_IVO3.DataCollection.Channel_088", shotN)
plt.figure()
plt.xlim([0,1000])
plt.title("Pulse #"+str(shotN)+ " Plasma Current")
plt.ylabel("Plasma current [A]")
plt.xlabel("Time [ms]")
plt.plot(times_ip*1e-3,ip, label="radial upper")
plt.tight_layout()
