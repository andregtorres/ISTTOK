#Andre Torres 10-11-2018
#run in terminal as root
#computes the 'real' position of the horizontal field coils
from __future__ import print_function
from field import *
from getMirnov import *
from scipy.constants import mu_0
import keyboard

#SDAS
shotH=44330
#Coil signals
hor, times, tbs = getSignal(ch_hor, shotH)
#mirnov signals
times, dataH = getMirnovs(shotH,mirnv,True)
iH=np.mean(hor[1100:3900])
ftH=-flatTops(dataH,2900,3900) #MINUS SIGN!
initial=[[58.,58.],[-7.,7.],[+1.,-1.]]
H0=getMirnovFlux([58.,58.],[-7.,7.],[+1.,-1.],4,biotSavart=False)*iH

print("START")
plt.ion()
x=np.arange(12)+1
fig=plt.figure()
plt.tight_layout()
ax=fig.add_subplot(211)
plt.ylabel("Mirnov Flux [uV s]")
line0,=ax.plot(x,ftH*1e6, label="Measurement")
line1,=ax.plot(x,H0*1e6, label="Original positions")
line2,=ax.plot(x,H0*1e6, label="Optimized positions")
plt.legend()

ax2=fig.add_subplot(212, sharex=ax)
plt.xlabel("Mirnov probe")
plt.ylabel("Flux difference [uV s]")
line20,=ax2.plot(x,(ftH-ftH)*1e6)
line21,=ax2.plot(x,(ftH-H0)*1e6)
line22,=ax2.plot(x,(ftH-H0)*1e6)
d=[[0,0],[0,0],[0,0],[1.,1.]]

coil=1
while( not keyboard.is_pressed("q")):
    if keyboard.is_pressed("1"):
        coil=1
    if keyboard.is_pressed("2"):
        coil=2
    if keyboard.is_pressed("3"):
        coil=3
    if keyboard.is_pressed("4"):
        coil=4
    if keyboard.is_pressed("5"):
        coil=5

    if keyboard.is_pressed("a"):
        d[coil-1][0]+=0.1
    if keyboard.is_pressed("z"):
        d[coil-1][0]-=0.1
    if keyboard.is_pressed("s"):
        d[coil-1][0]+=0.5
    if keyboard.is_pressed("x"):
        d[coil-1][0]-=0.5
    if keyboard.is_pressed("d"):
        d[coil-1][1]+=0.1
    if keyboard.is_pressed("c"):
        d[coil-1][1]-=0.1
    if keyboard.is_pressed("f"):
        d[coil-1][1]+=0.5
    if keyboard.is_pressed("v"):
        d[coil-1][1]-=0.5
    if keyboard.is_pressed(" "):
        print (d)
    H=getMirnovFluxCorrected([initial[0][0]+d[0][0],initial[0][1]+d[1][0]],[initial[1][0]+d[0][1],initial[1][1]+d[1][1]],[+1.,-1.],(4.+d[2][1]),d[3],biotSavart=False)*(iH+d[2][0])

    line22.set_ydata((ftH-H)*1e6)
    line2.set_ydata(H*1e6)
    fig.canvas.draw()
    fig.canvas.flush_events()
print (d)
d=[[-2.0, -3.0], [-4.3, 5.7], [47.5, 0]]
d2=[[-0.6, -5.1], [-1.7, 9.0], [167.5, 0]]
