#Andre Torres - 23/07/2018
#computes horizontal and vertical displacement from the mirnovs

import matplotlib.pyplot as plt
from getSdasSignal import *

mirnv=['MARTE_NODE_IVO3.DataCollection.Channel_129',
'MARTE_NODE_IVO3.DataCollection.Channel_130',
'MARTE_NODE_IVO3.DataCollection.Channel_131',
'MARTE_NODE_IVO3.DataCollection.Channel_132',
'MARTE_NODE_IVO3.DataCollection.Channel_133',
'MARTE_NODE_IVO3.DataCollection.Channel_134',
'MARTE_NODE_IVO3.DataCollection.Channel_135',
'MARTE_NODE_IVO3.DataCollection.Channel_136',
'MARTE_NODE_IVO3.DataCollection.Channel_137',
'MARTE_NODE_IVO3.DataCollection.Channel_138',
'MARTE_NODE_IVO3.DataCollection.Channel_139',
'MARTE_NODE_IVO3.DataCollection.Channel_140']

mirnv_raw=['MARTE_NODE_IVO3.DataCollection.Channel_145',
'MARTE_NODE_IVO3.DataCollection.Channel_146',
'MARTE_NODE_IVO3.DataCollection.Channel_147',
'MARTE_NODE_IVO3.DataCollection.Channel_148',
'MARTE_NODE_IVO3.DataCollection.Channel_149',
'MARTE_NODE_IVO3.DataCollection.Channel_150',
'MARTE_NODE_IVO3.DataCollection.Channel_151',
'MARTE_NODE_IVO3.DataCollection.Channel_152',
'MARTE_NODE_IVO3.DataCollection.Channel_153',
'MARTE_NODE_IVO3.DataCollection.Channel_154',
'MARTE_NODE_IVO3.DataCollection.Channel_155',
'MARTE_NODE_IVO3.DataCollection.Channel_156']

mirnv_corr=['MARTE_NODE_IVO3.DataCollection.Channel_166',
'MARTE_NODE_IVO3.DataCollection.Channel_167',
'MARTE_NODE_IVO3.DataCollection.Channel_168',
'MARTE_NODE_IVO3.DataCollection.Channel_169',
'MARTE_NODE_IVO3.DataCollection.Channel_170',
'MARTE_NODE_IVO3.DataCollection.Channel_171',
'MARTE_NODE_IVO3.DataCollection.Channel_172',
'MARTE_NODE_IVO3.DataCollection.Channel_173',
'MARTE_NODE_IVO3.DataCollection.Channel_174',
'MARTE_NODE_IVO3.DataCollection.Channel_175',
'MARTE_NODE_IVO3.DataCollection.Channel_176',
'MARTE_NODE_IVO3.DataCollection.Channel_177']

ch_prim='MARTE_NODE_IVO3.DataCollection.Channel_093';
ch_hor='MARTE_NODE_IVO3.DataCollection.Channel_091';
ch_vert='MARTE_NODE_IVO3.DataCollection.Channel_092';

ch_Ip_rog='MARTE_NODE_IVO3.DataCollection.Channel_088';
ch_chopper='MARTE_NODE_IVO3.DataCollection.Channel_141';

#SAVES MIRNOV DATA IN A LIST OF NP.ARRAYS
def getMirnovs(shot_, node=mirnv_corr, correct=True):
    coilNr=0
    data=[]
    for coil in node:
        coilNr+=1
        coilData, times, tbs = getSignal(coil, shot_)
        if correct:
            f=100 #correction length for slope calculation
            slope=0.0
            if node==mirnv:
                slope=np.linspace(np.mean(coilData[0:f]), np.mean(coilData[-f-1:-1]), num=len(coilData))
            if coilNr in [1,2,4,11]:
                coilData=(coilData-slope) #positive polarity
            else:
                coilData=-(coilData-slope) #negative polarity
        if shot_ > 44000:     # Check this exact number
            data.append(coilData * 0.85e-10 /10.0 ) # Return values in V.s units
        else:
            data.append(coilData * 0.85e-10 ) # Return values in V.s units
    return times, data

#SAVES MIRNOV DATA IN A LIST OF NP.ARRAYS
def getMirnovs2(shot_, node=mirnv_corr, ADC=False, polarity=True, slope=False):
    coilNr=0
    data=[]
    pol=1.
    adc_factor=1.
    if ADC:
        adc_factor=0.85e-10
    for coil in node:
        coilNr+=1
        #polarity correction
        if coilNr  not in [1,2,4,11] and polarity:
            pol=-1.
        else:
            pol=1.
        coilData, times, tbs = getSignal(coil, shot_)
        #correct slope
        sl=np.zeros(len(coilData))
        if slope and node==mirnv:
            f=100 #correction length for slope calculation
            sl=np.linspace(np.mean(coilData[0:f]), np.mean(coilData[-f-1:-1]), num=len(coilData))
        coilData=pol*(coilData-sl)*adc_factor #positive polarity
        if shot_ > 44000:     # Check this exact number
            data.append(coilData * 0.85e-10 /11.0 ) # Return values in V.s units
        else:
            data.append(coilData * 0.85e-10 ) # Return values in V.s units
    return times, data

#PLOTS ALL DATA FROM MIRNOVS
def plotAll(times_, data_, show=True, title=''):
    plt.figure()
    coilNr=0
    plt.suptitle(title)
    ax=[]
    for coil in data_:
        coilNr+=1
        ax.append( plt.subplot(3, 4, coilNr))
        ax[-1].set_title("MIRNOV #"+str(coilNr))
        ax[-1].ticklabel_format(style='sci',axis='y', scilimits=(0,0))
        plt.plot(times_*1e-3, coil)

    ax[0].get_shared_x_axes().join(ax[4], ax[8])
    if show:
        plt.show()

#PLOTS ONE MIRNOV
def plotMirnov(times_, data_, show=True, title=''):
    plt.figure()
    plt.title(title)
    #ax[-1].ticklabel_format(style='sci',axis='y', scilimits=(0,0))
    plt.plot(times_*1e-3, data_)
    if show:
        plt.show()

if __name__ == "__main__":
    #horizontal coils
    plotAll(*getMirnovs(42966,mirnv,True), show=True, title="Horizontal Field Coils")

'''
times_old, data_old = getMirnovs(42966,mirnv,True)
times_new, data_new = getMirnovs(44330,mirnv,True)
%matplotlib qt4
n=1
sty="-"
for coil_old, coil_new in zip(data_old, data_new):
    if n>10: sty="+"
    plt.figure(1)
    plt.plot(times_old,coil_old, sty, label="c"+str(n))
    plt.figure(2)
    plt.plot(times_new,coil_new, sty, label="c"+str(n))
    n+=1
plt.figure(1)
plt.legend()
plt.figure(2)
plt.legend()
'''
