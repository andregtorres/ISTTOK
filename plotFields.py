# PLOTS THE FIELD GENERATED BY THE PF COILS
# Andre Torres
#18.01.19

from field import getPFFlux2
from coilDefinitions import PF0, PF1, PF2, v, ht, hb
import matplotlib.pyplot as plt
import numpy as np

#%matplotlib qt4

def plotPF(PF, title ="", probes=True):
    nr=50
    nz=50
    if probes:
        rlim=[0.5,0.75]
        zlim=[-0.1,0.1]
    else:
        rlim=[0.25,0.65]
        zlim=[-0.2,0.2]
    R=np.linspace(rlim[0],rlim[1], nr)
    Z=np.linspace(zlim[0],zlim[1], nz)
    Br=np.zeros(shape=(nz,nr))
    Bz=np.zeros(shape=(nz,nr))

    for i in range(len(R)):
        for j in range(len(Z)):
            Br[j][i], Bz[j][i] = getPFFlux2(R[i],Z[j],PF, biotSavart=False)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    plt.xlabel("R [m]")
    plt.ylabel("Z [m]")
    plt.xlim(rlim)
    plt.ylim(zlim)
    plt.title(title)
    limiter = plt.Circle((0.46, 0), 0.0935, color='k', linestyle="-", fill=False, linewidth=3)
    ax.add_artist(limiter)
    if probes:
        vcoil=plt.plot([v.r-v.w/2., v.r+v.w/2.], [v.z,v.z], "r", linewidth=2)
        htcoil=plt.plot([ht.r, ht.r], [ht.z-ht.w/2.,ht.z+ht.w/2.], "g",linewidth=2)
        hbcoil=plt.plot([hb.r, hb.r], [hb.z-hb.w/2.,hb.z+hb.w/2.], "b",linewidth=2)
    ctr=plt.contourf(R,Z,np.sqrt(Bz**2+Br**2)*1e3, cmap="GnBu")#, levels=np.linspace(-1e-2,0,200)
    plt.streamplot(R,Z,Br,Bz, color="k")
    cbar=plt.colorbar(ctr)
    cbar.set_label("[mT/A]")
    ax.set_aspect('equal')

plotPF(PF0[0],"Primary, PF0")
plotPF(PF0[1],"Vertical, PF0")
plotPF(PF0[2],"Horizontal, PF0")

plotPF(PF0[0],"Primary", False)
plotPF(PF0[1],"Vertical PF", False)
plotPF(PF0[2],"Horizontal PF", False)

if __name__ == '__main__':
    plt.show()
