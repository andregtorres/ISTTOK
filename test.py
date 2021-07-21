from getMirnov import *
import numpy as np
#SDAS
shotV=42952
shotH=42966
shotP=43066
#Coil signals
vert, times, tbs = getSignal(ch_vert, shotV )
plt.plot(times,vert)
#mirnov signals
times, mirnvV = getMirnovs(shotV,mirnv,True)
plt.plot(times,mirnvV[0])
