#Andre Torres 5-11-2018
#run in terminal as root
#computes the 'real' position of the primary field coils
from field import *
from getMirnov import *
from scipy.constants import mu_0
import keyboard
#SDAS
shotP=44501
#Coil signals
prim, times, tbs = getSignal(ch_prim, shotP)
#mirnov signals
times, dataP = getMirnovs(shotP,mirnv,True)

iV=np.mean(prim[2100:5900])
ftP=-flatTops(dataP,5000,6000) #MINUS SIGN!
P0=getMirnovFlux([62.,62.],[-13.,13.],[+1.,1.],14,biotSavart=False)*iV

print("START")
plt.ion()
x=np.arange(12)+1
fig=plt.figure()
ax=fig.add_subplot(211)
plt.tight_layout()
plt.ylabel("Mirnov Flux [uV s]")
line0,=ax.plot(x,ftP*1e6, label="Measurement")
line1,=ax.plot(x,P0*1e6, label="Original positions")
line2,=ax.plot(x,P0*1e6, label="Optimized positions")
plt.legend()

ax2=fig.add_subplot(212, sharex=ax)
plt.xlabel("Mirnov probe")
plt.ylabel("Flux difference [uV s]")
line20,=ax2.plot(x,(ftP-ftP)*1e6)
line21,=ax2.plot(x,(ftP-P0)*1e6)
line22,=ax2.plot(x,(ftP-P0)*1e6)

d=[[0,-1],[0,+1],[0,0],[1.,1.]]
initial=[[62.,62.],[-13.,13.]]
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
        if coil in [1,2]:
            d[0][1]+=0.1
            d[1][1]+=0.1
        else:
            d[coil-1][1]+=0.1
    if keyboard.is_pressed("c"):
        if coil in [1,2]:
            d[0][1]-=0.1
            d[1][1]-=0.1
        else:
            d[coil-1][1]-=0.1
    if keyboard.is_pressed("f"):
        if coil in [1,2]:
            d[0][1]+=0.5
            d[1][1]+=0.5
        else:
            d[coil-1][1]+=0.5
    if keyboard.is_pressed("v"):
        if coil in [1,2]:
            d[0][1]-=0.5
            d[1][1]-=0.5
        else:
            d[coil-1][1]-=0.5
    if keyboard.is_pressed(" "):
        print d
    P=getMirnovFluxCorrected([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,d[3],biotSavart=False)*(iV+d[2][0])
    line22.set_ydata((ftP-P)*1e6)
    line2.set_ydata(P*1e6)
    fig.canvas.draw()
    fig.canvas.flush_events()
print d
#[-0.1, -3.9], [2.7, 5.9], [3.5, -9.1], [2.4, 7.1] EM relacao a optimizacao horizontal
#[-2.9, -6.2], [-2.6, 9.7], [3.5, -8.2], [2.2, 5.7]
'''
#plots

d=[[-1.3, 1.7], [-2.6, 0.7], [-38., 0]] #6/11
Pold=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=False)*(iV+d[2][0])

d=[[-1.4, -1.1], [-1.4, 0.9], [-24.7, 0]]
Pnew=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=False)*(iV+d[2][0])
Pnew2=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=True)*(iV+d[2][0])

d=[[-0.5, -1.0], [-0.5, 1.0], [-15, 0]]
Pprof=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=False)*(iV+d[2][0])

d=[[-0, -0], [-0, 0], [0, 0]]
P0=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=False)*(iV+d[2][0])
P0BS=getMirnovFlux([62.+d[0][0],62.+d[1][0]],[-13.+d[0][1],13.+d[1][1]],[+1.,1.],14,biotSavart=True)*(iV+d[2][0])

plt.figure()
plt.plot(np.arange(12)+1,Pold*1e6 -ftP*1e6, label="Old optimization -Measured flux")
plt.plot(np.arange(12)+1,Pnew*1e6 -ftP*1e6, label="New optimization -Measured flux")
plt.plot(np.arange(12)+1,Pprof*1e6 -ftP*1e6, label="Professors optimization -Measured flux")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44501 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()
plt.savefig("plots/PrimOptim_44501_comparison_rel2.png")

plt.figure()
plt.plot(np.arange(12)+1,ftP*1e6, "-", label="Measured flux")
plt.plot(np.arange(12)+1,Pnew*1e6,"-", label="New optimization")
plt.plot(np.arange(12)+1,Pnew2*1e6,"-", label="New optimization BS")
plt.plot(np.arange(12)+1,P0BS*1e6,"-", label="NO optimization BS")
plt.plot(np.arange(12)+1,P0*1e6,"-", label="NO optimization")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44501 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()



plt.figure()
plt.plot(np.arange(11)+1,np.diff(ftP)*1e6, "x", label="Measured flux")
plt.plot(np.arange(11)+1,np.diff(Pold)*1e6,"x", label="Old optimization")
plt.plot(np.arange(11)+1,np.diff(Pnew)*1e6,"x", label="New optimization")
plt.plot(np.arange(11)+1,np.diff(Pprof)*1e6, "x",label="Professors optimization")
plt.xlabel("Mirnov probe")
plt.ylabel("Mirnov Flux [uV s]")
plt.title ("Pulse #44501 - Mirnov flux with optimized coil position")
plt.legend()
plt.tight_layout()
plt.savefig("plots/PrimOptim_44501_comparison_diff2.png")


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
'''
