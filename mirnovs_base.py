#Andre Torres 25-01-2019
#Base for mirnov analysis
from coilDefinitions import PF0, PF1, PF2, diagCoil
from getMirnov import *
from scipy.constants import mu_0
import keyboard
#SDAS
shotV=44278 #44278   42952
#Coil signals
vert, times, tbs = getSignal(ch_vert, shotV )
#mirnov signals
times, dataV = getMirnovs(shotV,mirnv,True)

iV=np.mean(vert[2100:5900])
print (iV)
ftV=flatTops(dataV,5000,6000)
V0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
