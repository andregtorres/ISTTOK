#Andre Torres - 27/07/2018
#computes horizontal and vertical displacement from the mirnovs

from getMirnov import *
import scipy.signal as signal
from scipy.optimize import *
from filters import CSfilter

#43066 e 43063
shotnr=43066;

#Coil signals
prim, times, tbs = getSignal(ch_prim, shotnr )
prim2, times2, tbs = getSignal(ch_prim, 43063 )
hor, times, tbs = getSignal(ch_hor, shotnr )
vert, times, tbs = getSignal(ch_vert, shotnr )

#slicing
slice_start=np.where(times==100000)[0][0]
slice_end=np.where(times==400000)[0][0]


#mirnov signals
times, data = getMirnovs(shotnr,mirnv,True)

plt.figure()
plt.plot(times, prim)
plt.plot(times2, prim2)
plt.show

def exponential(x, a, b, c):
    return a*(1.-np.exp(-x / b))+c

guess = [5.94e-6,  1.0313919e+05, 0]
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
params=np.zeros([12,7])
ax=[]
for coil in data:
    coilNr+=1
    ax.append( plt.subplot(3, 4, coilNr))
    ax[-1].set_title("MIRNOV #"+str(coilNr))
    ax[-1].ticklabel_format(style='sci',axis='y', scilimits=(0,0))
    plt.plot(times*1e-3, coil)

    #FIT
    popt, pcov = curve_fit(exponential, times[slice_start:slice_end], coil[slice_start:slice_end], p0=guessM[coilNr-1],  maxfev=50000) #guessM[coilNr-1]
    #Calculate R squared
    residuals = coil[slice_start:slice_end] - exponential(times[slice_start:slice_end], *popt)
    #Sum of the residuals squared
    ss_res = np.sum(residuals**2)
    #Total sum of squares
    ss_tot = np.sum((coil[slice_start:slice_end]-np.mean(coil[slice_start:slice_end]))**2)
    #R-Squared
    Rsq = 1.0 - ss_res/ss_tot

    plt.plot(times[slice_start:slice_end]*1e-3, exponential(times[slice_start:slice_end], *popt), label="fit" )
    #guess=popt  #guess based on last fit
    print "MIRNOV "+str(coilNr)+ " tau: "+str(popt[1]*1e-3)+" ms"+" fc="+str(1./popt[1]/2./np.pi*1e6)+" Hz R2="+str(Rsq)
    filtered=CSfilter(prim,1./popt[1]/2./np.pi, tbs)
    scale=(popt[0]+popt[2])/max(filtered)
    params[coilNr-1]=np.array([popt[0], popt[1], popt[2], 1./popt[1], 1./popt[1]/2./np.pi*1e6, Rsq, scale])
    plt.plot(times*1e-3, filtered*scale)

print "AVG:"
print "p0, p1, p2, 1/tau (us^-1), fc (Hz), R2, scale (Vs/A)"
print np.average(params, axis=0)
np.savetxt('primFits.txt', params)
plt.show()
