#Definitions of coils
#Andre Torres 16.01.19

import numpy as np
from field import getPFFlux2
from filters import CSfilter

#define the active coils
H0=[[58.*1e-2,58.*1e-2],[-7*1e-2,7*1e-2],[-1.,1.],4]
V0=[[0.58,0.58,0.35,0.35],[-0.07,0.07,-0.07,0.07],[-1.,-1.,1.,1.],5]
P0=[[0.62,0.62],[-0.13,0.13],[-1.,-1.],14]

H1=[[.544,.545],[-0.108,0.134],[-1.,1.],4]
V1=[[0.551,0.554,0.385,0.372],[-0.132,0.167,-0.152,0.127],[-1.,-1.,1.,1.],5]
P1=[[0.615,0.615],[-0.144,0.145],[-1.,-1.],14]

H2=[[.559,.588],[-0.100,0.146],[-1.,1.],4,1.324]
V2=[[0.5547,0.5264,0.3905,0.3803],[-0.1125,0.1268,-0.1204,0.1089],[-1.,-1.,1.,1.],5,0.83]
P2=[[0.606,0.606],[-0.141,0.139],[-1.,-1.],14,0.842]

PF0=[P0,V0,H0]
PF1=[P1,V1,H1]
PF2=[P2,V2,H2]

#define the FluxCoils
class diagCoil:
    #angle of the normal defined in deg with the vertical
    def __init__(self, r=0,z=0,l=0,w=0, turns=1, angle=0):
        self.r=r
        self.z=z
        self.l=l
        self.w=w
        self.area=l*w
        self.turns=turns
        self.angle=angle
        self.PFcontribution=[[],[],[]]

    def getSignal(self,PF=PF0,Iprim=0,Ivert=0,Ihor=0):
        HrP, HzP= getPFFlux2(self.r,self.z,PF[0])
        HrV, HzV= getPFFlux2(self.r,self.z,PF[1])
        HrH, HzH= getPFFlux2(self.r,self.z,PF[2])
        Hz=(HzV*Ivert+ HzP*Iprim + HzH*Ihor)
        Hr=(HrV*Ivert+ HrP*Iprim + HrH*Ihor)
        signal=Hz*np.cos(np.radians(self.angle))-Hr*np.sin(np.radians(self.angle))
        return signal*self.area*self.turns


    #computes and sets the signal for the PF optimizations
    def setSignals(self,Iprim=0,Ivert=0,Ihor=0):
        self.PF0=self.getSignal(PF0,Iprim,Ivert,Ihor)
        self.PF1=self.getSignal(PF1,Iprim,Ivert,Ihor)
        self.PF2=self.getSignal(PF2,Iprim,Ivert,Ihor)
        self.PF2g=self.getSignal(PF2,Iprim,Ivert,Ihor)*P2[4]

    def setSignals2(self, gains, Iprim,Ivert,Ihor):
        self.PF0=self.getSignal(PF0,Iprim*gains[0][0],Ivert*gains[0][1],Ihor*gains[0][2])
        self.PF1=self.getSignal(PF1,Iprim*gains[1][0],Ivert*gains[1][1],Ihor*gains[1][2])

class tripleCoil:
    def __init__(self, currents, amplify=False, applyFilter=False):
        self.v=diagCoil(r=0.706, z=0.0,l=0.18483,w=0.028, turns=10,angle=0)
        self.ht=diagCoil(r=0.7145, z=0.014,l=0.15588,w=0.014, turns=10,angle=-90)
        self.hb=diagCoil(r=0.7145, z=-0.014,l=0.15588,w=0.014, turns=10,angle=-90)

        currents_m=[currents,currents,currents] #probes
        if applyFilter:
            tbs=100.
            popt_P=[[7.89588047e-05, 2.43622814e+04, 6.55629815e-05],[-2.06196145e-06,  4.38844884e+04, -1.60822894e-06],[-8.12553610e-06,  1.27860286e+04, -6.02436304e-06]]
            for i in range(3): #over the probes
                fc1=1./popt_P[i][1]/2./np.pi
                currents_m[i][0]=np.asarray(CSfilter(currents[0],fc1, tbs))
            popt_V=[[-7.61731222e-03,  3.73277163e+04, -3.64464177e-05], [ 1.05060429e-04,  3.97652858e+04, -9.12632314e-07], [-4.73349226e-05,  6.12193708e+04, -1.29608285e-06]]
            for i in range(3): #over the probes
                fc1=1./popt_V[i][1]/2./np.pi
                currents_m[i][1]=np.asarray(CSfilter(currents[1],fc1, tbs))
            popt_H=[[-2.53755150e-02,  2.28065608e+04, -7.28439496e-06],[7.30941300e+05, 8.05266976e+03, 2.64429245e-05],[1.66489804e+06, 7.74906454e+03, 2.30214087e-05]]
            for i in range(3): #over the probes
                fc1=1./popt_H[i][1]/2./np.pi
                currents_m[i][2]=np.asarray(CSfilter(currents[2],fc1, tbs))

        if not amplify:
            self.v.setSignals(*currents_m[0])
            self.ht.setSignals(*currents_m[1])
            self.hb.setSignals(*currents_m[2])
        else:
            ratios_P_PF0=np.asarray([ 1.06403604, -0.53502599,  1.90483049])
            ratios_P_PF1=np.asarray([ 1.37882102, -0.57252196,  1.90218157])
            ratios_V_PF0=np.asarray([ 0.27341249, -0.19930837,  0.28213556])
            ratios_V_PF1=np.asarray([1.14112469, 0.38179557, 1.0416212 ])
            ratios_H_PF0=np.asarray([-796.92962351,    1.75420108,    1.52989164])
            ratios_H_PF1=np.asarray([2.06270312, 2.03482828, 1.82948384])

            gains_v=[[ratios_P_PF0[0],ratios_V_PF0[0],ratios_H_PF0[0]],[ratios_P_PF1[0],ratios_V_PF1[0],ratios_H_PF1[0]]]
            gains_ht=[[ratios_P_PF0[1],ratios_V_PF0[1],ratios_H_PF0[1]],[ratios_P_PF1[1],ratios_V_PF1[1],ratios_H_PF1[1]]]
            gains_hb=[[ratios_P_PF0[2],ratios_V_PF0[2],ratios_H_PF0[2]],[ratios_P_PF1[2],ratios_V_PF1[2],ratios_H_PF1[2]]]

            self.v.setSignals2(gains_v,*currents_m[0])
            self.ht.setSignals2(gains_ht,*currents_m[1])
            self.hb.setSignals2(gains_hb,*currents_m[2])

v=diagCoil(r=0.706, z=0.0,l=0.18483,w=0.028, turns=10,angle=0)
ht=diagCoil(r=0.7145, z=0.014,l=0.15588,w=0.014, turns=10,angle=-90)
hb=diagCoil(r=0.7145, z=-0.014,l=0.15588,w=0.014, turns=10,angle=-90)
