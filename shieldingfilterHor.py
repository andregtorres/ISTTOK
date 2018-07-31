#Andre Torres - 27/07/2018
#computes horizontal and horical displacement from the mirnovs

from getMirnov import *
import scipy.signal as signal
from scipy.optimize import curve_fit
from filters import CSfilter

shotVer=42952
shotHor=42966
#Coil signals
prim, times, tbs = getSignal(ch_prim, shotHor )
vert, times, tbs = getSignal(ch_vert, shotHor )
hor, times, tbs = getSignal(ch_hor, shotHor )

#slicing
slice_start=np.where(times==100000)[0][0]
slice_end=np.where(times==200000)[0][0]

#mirnov signals
times, data = getMirnovs(shotHor,mirnv,True)

def exponential(x, a, b, c):
    return a*(1.-np.exp(-x / b))+c

guess = [1.88230074e-04,  2.68313919e+04, -1.82439561e-04]
guessM = [[ 80.7919370901674 , 6590.632343915145 , -80.79191856826026 ],
[ 2.0077302113187816 , 8899.817967723295 , -2.007706843416215 ],
[ 4.208050151502539 , 7322.987193204733 , -4.208045228622198 ],
[ 57.501582301707586 , 6516.179752199343 , -57.50157218465615 ],
[ 57.501582301707586 , 6516.179752199343 , -57.50157218465615 ],
[ 57.501582301707586 , 6516.179752199343 , -57.50157218465615 ],
[ 57.501582301707586 , 6516.179752199343 , -57.50157218465615 ],
[ 57.501582301707586 , 6516.179752199343 , -57.50157218465615 ],
[ -131.40718157126113 , 6189.509534048303 , 131.4071711802073 ],
[ -0.0013073426039006074 , 27831.2345744908 , 0.0012732434791092925 ],
[ -12.559165108129614 , 7164.670256164515 , 12.559155376979342 ],
[ -156.7688623998181 , 6143.6995925094525 , 156.76885071336318 ]]

plt.figure()
coilNr=0
params=np.zeros([12,7])
ax=[]
for coil in data:
    coilNr+=1
    ax.append(plt.subplot(3, 4, coilNr))
    ax[-1].set_title("MIRNOV #"+str(coilNr))
    ax[-1].ticklabel_format(style='sci',axis='y', scilimits=(0,0))
    plt.plot(times*1e-3, coil)
    #FIT
    if coilNr not in [5,6,7,8]:
        popt, pcov = curve_fit(exponential, times[slice_start:slice_end], coil[slice_start:slice_end], p0=guessM[coilNr-1],  maxfev=50000)
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
        #if coilNr == 4: guess[0]=-guess[0]
    #print "[",popt[0],",",popt[1],",",popt[2],"],"
    print "MIRNOV "+str(coilNr)+ " tau: "+str(popt[1]*1e-3)+" ms"+" fc="+str(1./popt[1]/2./np.pi*1e6)+" Hz R2="+str(Rsq)
    filtered=CSfilter(hor,1./popt[1]/2./np.pi, tbs)
    scale=(popt[0]+popt[2])/max(filtered)
    params[coilNr-1]=np.array([popt[0], popt[1], popt[2], 1./popt[1], 1./popt[1]/2./np.pi*1e6, Rsq, scale])
    if coilNr not in [5,6,7,8]: plt.plot(times*1e-3, filtered*scale)

print "AVG:"
print "p0, p1, p2, 1/tau (us^-1), fc (Hz), R2, scale (Vs/A)"
print np.average(params, axis=0)
np.savetxt('horFits.txt', params)
plt.show()
