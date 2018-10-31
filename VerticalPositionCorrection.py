#Andre Torres 31-10-2018
#run in terminal as root
#computes the 'real' position of the vertical field coils
from field import *
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
ftV=flatTops(dataV,5000,6000)
V0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
#guesses
Hguess=[[55.5,-9.75],[52.6,11.2]]
V1=getMirnovFlux([Hguess[0][0],Hguess[1][0],35.,35.],[Hguess[0][1],Hguess[1][1],-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
V2=getMirnovFlux([Hguess[0][0],Hguess[1][0],35.-4,35.-3],[Hguess[0][1],Hguess[1][1],-7-2.,7.-2],[-1.,-1.,1.,1.],5,biotSavart=False)*iV

print("START")
plt.ion()
fig=plt.figure()
ax=fig.add_subplot(111)
line0,=ax.plot(ftV)
line1,=ax.plot(V0)
line2,=ax.plot(V0)
d=[[0,0],[0,0],[0,0],[0,0]]
while( not keyboard.is_pressed("q")):
    if keyboard.is_pressed("1"):
        coil=1
    if keyboard.is_pressed("2"):
        coil=2
    if keyboard.is_pressed("3"):
        coil=3
    if keyboard.is_pressed("4"):
        coil=4

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
    #V2=getMirnovFlux([Hguess[0][0]+d[0][0],Hguess[1][0]+d[1][0],35.+d[2][0],35.+d[3][0]],[Hguess[0][1]+d[0][1],Hguess[1][1]+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
    V2=getMirnovFlux([58.+d[0][0],58.+d[1][0],35.+d[2][0],35.+d[3][0]],[-7.+d[0][1],7.+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
    line2.set_ydata(V2)
    fig.canvas.draw()
    fig.canvas.flush_events()
print d
#[-0.1, -3.9], [2.7, 5.9], [3.5, -9.1], [2.4, 7.1] EM relacao a optimizacao horizontal
#[-2.9, -6.2], [-2.6, 9.7], [3.5, -8.2], [2.2, 5.7]

#plots
d=[[-2.9, -6.2], [-2.6, 9.7], [3.5, -8.2], [2.2, 5.7]]

V2=getMirnovFlux([58.+d[0][0],58.+d[1][0],35.+d[2][0],35.+d[3][0]],[-7.+d[0][1],7.+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
plt.figure()
plt.plot(np.arange(12)+1,ftV*1e6, label="Measured flux")
plt.plot(np.arange(12)+1,V0*1e6, label="Original PF positions")
plt.plot(np.arange(12)+1,V2*1e6, label="Optimized PF positions")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44278 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()
plt.savefig("plots/VerticalOptim_44278.png")

shotV=42952 #shot antigo
vert, times, tbs = getSignal(ch_vert, shotV )
times, dataV = getMirnovs(shotV,mirnv,True)
iV=np.mean(vert[1100:5900])
ftV=flatTops(dataV,4000,6000)
V0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
V2=getMirnovFlux([58.+d[0][0],58.+d[1][0],35.+d[2][0],35.+d[3][0]],[-7.+d[0][1],7.+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV

plt.figure()
plt.plot(np.arange(12)+1,-ftV*1e6, label="Measured flux")
plt.plot(np.arange(12)+1,V0*1e6, label="Original PF positions")
plt.plot(np.arange(12)+1,V2*1e6, label="Optimized PF positions")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #42952 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()
plt.savefig("plots/VerticalOptim_42952.png")

shotV=44473
vert, times, tbs = getSignal(ch_vert, shotV )
times, dataV = getMirnovs(shotV,mirnv,True)
plt.plot(vert)
iV=np.mean(vert[1100:5900])
ftV=flatTops(dataV,4000,6000)
V0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
V2=getMirnovFlux([58.+d[0][0],58.+d[1][0],35.+d[2][0],35.+d[3][0]],[-7.+d[0][1],7.+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV

plt.figure()
plt.plot(np.arange(12)+1,ftV*1e6, label="Measured flux")
plt.plot(np.arange(12)+1,V0*1e6, label="Original PF positions")
plt.plot(np.arange(12)+1,V2*1e6, label="Optimized PF positions")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44473 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()
plt.savefig("plots/VerticalOptim_44473.png")
