#Code for acessment of the external coils at equatorial plane
#Andre Torres
#21-12-18
from getMirnov import *
%matplotlib qt4

#SDAS INFO
shotN=44835 #44409
ch_rad_u = 'MARTE_NODE_IVO3.DataCollection.Channel_141'
ch_vertical= 'MARTE_NODE_IVO3.DataCollection.Channel_142'
ch_rad_b = 'MARTE_NODE_IVO3.DataCollection.Channel_143'
#client.searchParametersByName("plasma")

#reference signals
primary, times_p,tbs=getSignal( 'MARTE_NODE_IVO3.DataCollection.Channel_093', shotN)
PF_vert, times_v,tbs=getSignal( ch_vert, shotN)
density, times_rho,tbs=getSignal( "POST.PROCESSED.DENSITY", shotN)
times, mirnovs = getMirnovs(shotN,mirnv_corr,False)
mirnovs0=mirnovs[0]

#triple sadle
#ADC-Vs factor
vertScale =  1.7102e-4 / 2.0e6 # LSB to Volt * Sampling Period
rad_u, times,tbs=getSignal(ch_rad_u, shotN, vertScale)
rad_b, times,tbs=getSignal(ch_rad_b, shotN, vertScale)
vert, times,tbs=getSignal(ch_vertical, shotN, vertScale)

#save files for offline
np.save("dataFiles/FluxCoils/times", times)
np.save("dataFiles/FluxCoils/primary", primary)
np.save("dataFiles/FluxCoils/PF_vert", PF_vert)
np.save("dataFiles/FluxCoils/density", density)
np.save("dataFiles/FluxCoils/mirnovs0", mirnovs[0])
np.save("dataFiles/FluxCoils/rad_u", rad_u)
np.save("dataFiles/FluxCoils/rad_b", rad_b)
np.save("dataFiles/FluxCoils/vert",  vert)

#load files
times=np.load("dataFiles/FluxCoils/times.npy")
primary=np.load("dataFiles/FluxCoils/primary.npy")
PF_vert=np.load("dataFiles/FluxCoils/PF_vert.npy")
density=np.load("dataFiles/FluxCoils/density.npy")
mirnovs0=np.load("dataFiles/FluxCoils/mirnovs0.npy")
rad_u=np.load("dataFiles/FluxCoils/rad_u.npy")
rad_b=np.load("dataFiles/FluxCoils/rad_b.npy")
vert=np.load("dataFiles/FluxCoils/vert.npy")

#Plot 3 signals
plt.figure()
plt.title("Pulse #"+str(shotN))
plt.ylabel("Flux [V.s]")
plt.xlabel("Time [ms]")
plt.plot(times*1e-3,rad_u*1e6, label="Upper radial")
plt.plot(times*1e-3,rad_b*1e6, label="Lower radial")
plt.plot(times*1e-3,vert*1e6, label="Vertical")
plt.tight_layout()
plt.legend()



plt.plot(times,-mirnovs[0]*max(vert)/max(mirnovs[0])*1e6)

plt.figure()
plt.plot(times,mirnovs[0])


plt.figure()
plt.plot(times_p, primary)
plt.plot(times_v, PF_vert)
