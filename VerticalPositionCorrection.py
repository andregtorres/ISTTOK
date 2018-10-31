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
line1,=ax.plot(V1)
line2,=ax.plot(V2)
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
    V2=getMirnovFlux([Hguess[0][0]+d[0][0],Hguess[1][0]+d[1][0],35.+d[2][0],35.+d[3][0]],[Hguess[0][1]+d[0][1],Hguess[1][1]+d[1][1],-7+d[2][1],7.+d[3][1]],[-1.,-1.,1.,1.],5,biotSavart=False)*iV
    line2.set_ydata(V2)
    fig.canvas.draw()
    fig.canvas.flush_events()
print d
#[-0.1, -3.9], [2.7, 5.9], [3.5, -9.1], [2.4, 7.1]
