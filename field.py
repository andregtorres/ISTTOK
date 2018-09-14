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
    return hR,hZ,flux

def biotsavart(r,z,Rw,zw,I,N=400):
    point=np.array((r,0,z))
    B=np.zeros(3)
    for phi in np.linspace(0,2*np.pi,N):
        w=np.array([Rw*np.cos(phi),Rw*np.sin(phi),zw])
        dist=point-w
        d=np.sqrt(np.sum([i**2 for i in dist]))
        dl=np.array([-np.sin(phi), np.cos(phi), 0])*2.*np.pi*Rw/N
        cross=np.cross(dl,dist)
        B=B+(cross/d**3)
    B=B*mu_0/4./np.pi*I
    return np.sqrt(B[0]**2+B[1]**2), B[2]

def biotSavartCyl(r,z,Rw,Zw,I):
    d=np.sqrt((r-Rw)**2+ (z-Zw)**2)
    B=mu_0/4./np.pi*I/(d**2)
    alpha=np.arctan2((z-Zw),(r-Rw))
    Br=B*np.sin(alpha)
    Bz=B*np.cos(alpha)
    return Br, Bz

R=np.linspace(20,70)
Z=np.linspace(-10,10)
Br=np.zeros(50)
Bz=np.zeros(50)
for Rw,Zw, sign in zip([58.,58.,35.,35.],[-7.,7.,-7.,7.],[-1.,-1.,1.,1.]):
    for i in range(50):
        br,bz=biotSavartCyl(50*1e-2, Z[i]*1e-2, Rw*1e-2,Zw*1e-2,340*sign) #46.
        Br[i] = Br[i] + br
        Bz[i] = Bz[i] + bz

import matplotlib.pyplot as plt
plt.plot(Z,np.sqrt(Br**2+Bz**2)*1000)
plt.plot(Z,Br*1000)
plt.plot(Z,Bz*1000)
