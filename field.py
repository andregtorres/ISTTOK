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

def getMirnovFlux(Rw_,Zw_,polarity,windings,biotSavart=True):
    #mirnov positions
    radius=9.35 #cm
    angle=345. - 30.*np.arange(12)
    geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
    geometryR=radius*np.cos(np.radians(angle))
    #loop on the mirnovs
    Hr=np.zeros(len(angle))
    Hz=np.zeros(len(angle))
    i=0
    for r,z in zip(geometryR,geometryZ):
        #loop on the PFCs
        for Rw,Zw, sign in zip(Rw_,Zw_,polarity):
            if biotSavart:
                coilHr, coilHz= biotsavart((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2,1.0) #46.
            else:
                coilHr, coilHz= Hcoil((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2) #46.
            Hr[i]+=sign*coilHr
            Hz[i]+=sign*coilHz
        i+=1
    Hr=np.asarray(Hr)
    Hz=np.asarray(Hz)
    Hp=-Hr*np.sin(np.radians(angle))+Hz*np.cos(np.radians(angle))
    return Hp*windings*50*49e-6

#magnetic flux per unit area a set of PF coil at Rw_, Zw_ (array) with windings and polarity creates on the positions r,z
# coordinates in [m]
def getPFFlux(r,z,Rw_,Zw_,polarity,windings):
    #loop on the PFCs
    Hr=0
    Hz=0
    for Rw,Zw, sign in zip(Rw_,Zw_,polarity):
        coilHr, coilHz= Hcoil(r, z, Rw,Zw)
        Hr+=sign*coilHr
        Hz+=sign*coilHz
    return Hr*windings, Hz*windings
#assuming input as array
def getPFFlux2(r,z,PF, biotSavart=False):
    Rw_=PF[0]
    Zw_=PF[1]
    polarity=PF[2]
    windings=PF[3]
    #loop on the PFCs
    Hr=0
    Hz=0
    for Rw,Zw, sign in zip(Rw_,Zw_,polarity):
        if not biotSavart:
            coilHr, coilHz= Hcoil(r, z, Rw,Zw)
        else:
            coilHr, coilHz=biotsavart(r, z, Rw,Zw,1.0)
        Hr+=sign*coilHr
        Hz+=sign*coilHz
    return Hr*windings, Hz*windings

def getMirnovFluxCorrected(Rw_,Zw_,polarity,windings,correction,biotSavart=True):
    #mirnov positions
    radius=9.35 #cm
    angle=345. - 30.*np.arange(12)
    geometryZ=radius*np.sin(np.radians(angle)) #positions of the mirnovs
    geometryR=radius*np.cos(np.radians(angle))
    #loop on the mirnovs
    Hr=np.zeros(len(angle))
    Hz=np.zeros(len(angle))
    i=0
    for r,z in zip(geometryR,geometryZ):
        #loop on the PFCs
        for Rw,Zw, sign in zip(Rw_,Zw_,polarity):
            if biotSavart:
                coilHr, coilHz= biotsavart((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2,1.0) #46.
            else:
                coilHr, coilHz= Hcoil((r+46.)*1e-2, z*1e-2, Rw*1e-2,Zw*1e-2) #46.
            Hr[i]+=sign*coilHr*correction[0]
            Hz[i]+=sign*coilHz*correction[1]
        i+=1
    Hr=np.asarray(Hr)
    Hz=np.asarray(Hz)
    Hp=-Hr*np.sin(np.radians(angle))+Hz*np.cos(np.radians(angle))
    return Hp*windings*50*49e-6

#get the mirnov flat top value with heaviside shots
def flatTops (data,from_=4000, to_=6000):
    return np.asarray([np.mean(np.array(i)[from_:to_]) for i in data])
