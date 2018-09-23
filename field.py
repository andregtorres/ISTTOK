# Andre' Torres - 2/08/18
#computes the H poloidal field and flux on a point R,Z caused by a coil in Rw,Zw
#adapted from http://baco.ipfn.ist.utl.pt/magnetic/hRZ_espira.m

import numpy as np
from scipy.special import ellipk, ellipe
from scipy.constants import mu_0

def Hcoil (R, Z, Rw, Zw):
    dZ=Z-Zw
    ss=(Rw+R)**2+dZ**2
    ss_sqr=np.sqrt(ss)
    num_r=Rw**2 + R**2 + dZ**2
    num_z=2*Rw**2 - num_r
    den= ss - 4*Rw * R #dR**2+dz**2
    m = 4 * Rw * R / ss
    K=ellipk(m)
    E=ellipe(m)
    flux=ss_sqr*((1.-m/2.)*K-E)
    hR  = dZ  / R * (num_r / den * E - K) / ss_sqr  / 2. / np.pi
    hZ = (K + num_z / den * E ) / ss_sqr  / 2. / np.pi
    return hR*mu_0,hZ*mu_0


def biotsavart(r,z,Rw,zw,I,N=100):
    point=np.array((r,0,z))
    B=np.zeros(3)
    for phi in np.linspace(0,2*np.pi,N):
        w=np.array([Rw*np.cos(phi),Rw*np.sin(phi),zw])
        w1=np.array([Rw*np.cos(phi+(2*np.pi)/N),Rw*np.sin(phi+(2*np.pi)/N),zw])
        dist=point-w
        d=np.sqrt(np.sum([i**2 for i in dist]))
        dl=np.array([-np.sin(phi), np.cos(phi), 0])*2.*np.pi*Rw/N
        dl2=w1-w
        cross=np.cross(dl,dist)
        B=B+(cross/d**3)
    B=B*mu_0/4./np.pi*I
    return B[0], B[2]

def biotSavartCyl(r,z,Rw,Zw,I):
    d=np.sqrt((r-Rw)**2+ (z-Zw)**2)
    B=mu_0/4./np.pi*I/(d**2)
    alpha=np.arctan((z-Zw)/(r-Rw))
    Bz=B*np.sin(alpha)
    Br=-B*np.cos(alpha)
    return Br, Bz
'''
R=np.linspace(30,60, 20)
Z=np.linspace(-10,10, 20)
Br=np.zeros(shape=(20,20))
Bz=np.zeros(shape=(20,20))

for Rw,Zw, sign in zip([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.]):
    for i in range(20):
        print i
        for j in range(20):
            br,bz=biotsavart(R[i]*1e-2, Z[j]*1e-2, Rw*1e-2,Zw*1e-2,340*sign) #46.
            Br[i][j] = Br[i][j] + br
            Bz[i][j] = Bz[i][j] + bz


Br2=np.zeros(shape=(20,20))
Bz2=np.zeros(shape=(20,20))
for Rw,Zw, sign in zip([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.]):
    for i in range(20):
        print i
        for j in range(20):
            br2,bz2 =Hcoil(R[i]*1e-2, Z[j]*1e-2, Rw*1e-2,Zw*1e-2) #46.
            Br2[i][j] = Br2[i][j] + br2*340*sign
            Bz2[i][j] = Bz2[i][j] + bz2*340*sign

Br3=np.zeros(shape=(20,20))
Bz3=np.zeros(shape=(20,20))
for Rw,Zw, sign in zip([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.]):
    for i in range(20):
        print i
        for j in range(20):
            br3,bz3=biotSavartCyl(R[i]*1e-2, Z[j]*1e-2, Rw*1e-2,Zw*1e-2,340.*sign) #46.
            Br3[i][j] = Br3[i][j] + br3
            Bz3[i][j] = Bz3[i][j] + bz3
print Bz
print Bz2*mu_0
import matplotlib.pyplot as plt
print Bz2*mu_0
plt.contourf(R,Z,Bz2*mu_0, levels=np.linspace(-1e-2,0,200))
plt.plot(Bz[10]*1000)
plt.plot(Bz2[10]*mu_0*1000)
plt.contourf(R,Z,np.sqrt(Bz**2+Br**2))#, levels=np.linspace(-1e-2,0,200)
plt.quiver(Br2,Bz2)

plt.plot(Z,np.sqrt(Br**2+Bz**2)*1000)
plt.plot(Z,Br*1000)
plt.plot(Z,Bz*1000)
'''
