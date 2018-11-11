#Andre Torres 11-11-2018
#Plots the comparissoms f the optimized positions
from field import *
from getMirnov import *

#Optimizations
originalV=[[58.,58.,35.,35.],[-7.,7.,-7.,7.],[1.,5],[1.,1.]]
optimV1=[[55.1,55.4,38.5,37.2],[-13.2,16.7,-15.2,12.7],[1.,5],[1.,1.]]
optimV2=[[55.4752,52.6368,39.05,38.027],[-11.25,12.6782,-12.0378,10.8862],[1.3852,3],[1.,1.]]
optimV3=[[55.9,57.,37.8,36.8],[-14.,17.1,-16.3,13.8],[1.,5],[1.6,0.8]]

originalH=[[58.,58.],[-7.,7.],[1.,4],[1.,1.]]
optimH1=[[54.4,54.5],[-10.8,13.4],[1.,4],[1.,1.]]
#optimH2=[[56.0,55.0],[-2.7,12.7],[1.2669,4],[1.,1.]]
optimH2=[[55.9,53.8],[-10.0,14.2],[1.324,4],[1.,1.]]
optimH3=[[52.3,55.0],[-10.9,13.5],[0.910,4],[1.6,0.8]]

originalP=[[62.,62.],[-13.,13.],[1.,14],[1.,1.]]
optimP1=[[61.5,61.5],[-14.4,14.5],[1.,14],[1.,1.]]
optimP2=[[60.6,60.6],[-14.1,13.9],[0.842,14],[1.,1.]]
optimP3=[[61.8,62.2],[-13.9,15.2],[0.978,14],[1.6,0.8]]


#SDAS
shotsV=[44278, 44473, 44475]
iV=[393.96,339.61,-334.21]
ftV=[]
#mirnov signals
for shot in shotsV:
    times, dataV = getMirnovs(shot,mirnv,True)
    ftV.append(flatTops(dataV,5000,6000)*1e6)

shotsH=[44330, 44480, 44481]
iH=[178.19,178.63,-173.95]
ftH=[]
#mirnov signals
for shot in shotsH:
    times, dataH = getMirnovs(shot,mirnv,True)
    ftH.append(-flatTops(dataH,3100,3900)*1e6)

shotsP=[44501, 44499, 44503]
iP=[156.8,156.5,-159.6]
ftP=[]
#mirnov signals
for shot in shotsP:
    times, dataP = getMirnovs(shot,mirnv,True)
    ftP.append(-flatTops(dataP,5100,5900)*1e6)



x=np.arange(12)+1

#Vertical
for shot, i, ft in zip(shotsV,iV,ftV):
    V0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.],5,biotSavart=False)*i*1e6
    V1=getMirnovFluxCorrected(optimV1[0],optimV1[1],[-1.,-1.,1.,1.],optimV1[2][1],optimV1[3],biotSavart=False)*i*optimV1[2][0]*1e6
    V2=getMirnovFluxCorrected(optimV2[0],optimV2[1],[-1.,-1.,1.,1.],optimV2[2][1],optimV2[3],biotSavart=False)*i*optimV2[2][0]*1e6
    V3=getMirnovFluxCorrected(optimV3[0],optimV3[1],[-1.,-1.,1.,1.],optimV3[2][1],optimV3[3],biotSavart=False)*i*optimV3[2][0]*1e6
    plt.figure()
    plt.plot(x,ft,"o", label="Measured flux")
    plt.plot(x,V0, label="Original")
    plt.plot(x,V1, label="Optimization 1")
    plt.plot(x,V2, label="Optimization 2")
    plt.plot(x,V3, label="Optimization 3")
    plt.xlabel("Mirnov probe")
    plt.ylabel("Mirnov Flux [uV s]")
    plt.title ("Pulse "+str(shot)+" - Mirnov flux with optimized Vertical field coil positions")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/OptimV_"+str(shot)+".png",dpi=200)

#Horizontal
for shot, i, ft in zip(shotsH,iH,ftH):
    H0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[1.,-1.],14,biotSavart=False)*i*1e6
    H1=getMirnovFluxCorrected(optimH1[0],optimH1[1],[1.,-1.],optimH1[2][1],optimH1[3],biotSavart=False)*i*optimH1[2][0]*1e6
    H2=getMirnovFluxCorrected(optimH2[0],optimH2[1],[1.,-1.],optimH2[2][1],optimH2[3],biotSavart=False)*i*optimH2[2][0]*1e6
    H3=getMirnovFluxCorrected(optimH3[0],optimH3[1],[1.,-1.],optimH3[2][1],optimH3[3],biotSavart=False)*i*optimH3[2][0]*1e6
    plt.figure()
    plt.plot(x,ft,"o", label="Measured flux")
    plt.plot(x,H0, label="Original")
    plt.plot(x,H1, label="Optimization 1")
    plt.plot(x,H2, label="Optimization 2")
    plt.plot(x,H3, label="Optimization 3")
    plt.xlabel("Mirnov probe")
    plt.ylabel("Mirnov Flux [uV s]")
    plt.title ("Pulse "+str(shot)+" - Mirnov flux with optimized Horizontal field coil positions")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/OptimH_"+str(shot)+".png",dpi=200)

#Primary
for shot, i, ft in zip(shotsP,iP,ftP):
    P0=getMirnovFlux([58.,58.,35.,35.],[-7.,7.,-7.,7.],[1.,1.],14,biotSavart=False)*i*1e6
    P1=getMirnovFluxCorrected(optimP1[0],optimP1[1],[1.,1.],optimP1[2][1],optimP1[3],biotSavart=False)*i*optimP1[2][0]*1e6
    P2=getMirnovFluxCorrected(optimP2[0],optimP2[1],[1.,1.],optimP2[2][1],optimP2[3],biotSavart=False)*i*optimP2[2][0]*1e6
    P3=getMirnovFluxCorrected(optimP3[0],optimP3[1],[1.,1.],optimP3[2][1],optimP3[3],biotSavart=False)*i*optimP3[2][0]*1e6
    plt.figure()
    plt.plot(x,ft,"o", label="Measured flux")
    plt.plot(x,P0, label="Original")
    plt.plot(x,P1, label="Optimization 1")
    plt.plot(x,P2, label="Optimization 2")
    plt.plot(x,P3, label="Optimization 3")
    plt.xlabel("Mirnov probe")
    plt.ylabel("Mirnov Flux [uV s]")
    plt.title ("Pulse "+str(shot)+" - Mirnov flux with optimized Primary coils positions")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/OptimP_"+str(shot)+".png", dpi=200)
