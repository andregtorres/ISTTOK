#Andre Torres - 27/07/2018
#computes horizontal and vertical displacement from the mirnovs

from getMirnov import *
import scipy.signal as signal
from scipy.optimize import *
from shieldingfilter import CSfilter

shotnr=42952;

#Coil signals
prim, times, tbs = getSignal(ch_prim, shotnr )
hor, times, tbs = getSignal(ch_hor, shotnr )
vert, times, tbs = getSignal(ch_vert, shotnr )

#slicing
slice_start=np.where(times==100000)[0][0]
slice_end=np.where(times==200000)[0][0]


#mirnov signals
times, data = getMirnovs(shotnr,mirnv_corr,True)

def exponential(x, a, b, c):
    return a*(1.-np.exp(-x / b))+c

guess = [1.88230074e-04,  2.68313919e+04, -1.82439561e-04]
guessM = [[ 0.00018823010565143905 , 26831.390780271355 , -0.0001824395933670592 ],
[ 0.00045526020482325745 , 27988.340464523386 , -0.0004411626641315281 ],
[ 0.0008593054691817632 , 20059.108499765076 , -0.0008539245081783592 ],
[ -0.0015214804939659302 , 19456.892995528884 , 0.0015136440080124175 ],
[ -0.00011821245784153389 , 40438.64770346272 , 0.00010671136091748272 ],
[ -1.5818271879615867e-05 , 108672.80602646185 , 7.637094041139406e-06 ],
[ -1.7592712656177852e-05 , 118204.53771271136 , 7.976603686725276e-06 ],
[ -0.0007603947988702048 , 24709.405624924613 , 0.0007473401741001225 ],
[ -3.9961069925505015e-05 , 55387.59724011559 , 3.110028335283982e-05 ],
[ 0.00022674870872464836 , 28515.41411518887 , -0.00022041532927244783 ],
[ 0.00026508140311219 , 31355.542905561575 , -0.00025524675518042225 ],
[ 7.352997274425692e-05 , 35438.751251357215 , -6.89680383018832e-05 ]]

plt.figure()
coilNr=0
for coil in data:
    coilNr+=1
    #FIT
    popt, pcov = curve_fit(exponential, times[slice_start:slice_end], coil[slice_start:slice_end], p0=guessM[coilNr-1],  maxfev=50000)
    #Calculate R squared
    residuals = coil[slice_start:slice_end] - exponential(times[slice_start:slice_end], *popt)
    #Sum of the residuals squared
    ss_res = np.sum(residuals**2)
    #Total sum of squares
    ss_tot = np.sum((coil[slice_start:slice_end]-np.mean(coil[slice_start:slice_end]))**2)
    #R-Squared
    Rsq = 1.0 - ss_res/ss_tot
    ax = plt.subplot(3, 4, coilNr)
    plt.plot(times, coil)
    plt.plot(times[slice_start:slice_end], exponential(times[slice_start:slice_end], *popt), label="fit" )
    #guess=popt  #guess based on last fit
    print "MIRNOV "+str(coilNr)+ " tau: "+str(popt[1]*1e-3)+" ms"+" fc="+str(1./popt[1]/2./np.pi*1e6)+" Hz R2="+str(Rsq)
    filtered=CSfilter(vert,1./popt[1]/2./np.pi, tbs)
    if popt[0] > 0:
        scale=max(exponential(times[slice_start:slice_end], *popt))/max(filtered)
    else:
        scale=min(exponential(times[slice_start:slice_end], *popt))/max(filtered)
    plt.plot(times, filtered*scale)

plt.show()
