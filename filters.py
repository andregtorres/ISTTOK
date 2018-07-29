#Andre Torres - 23/07/2018
#computes horizontal and vertical displacement from the mirnovs

from getMirnov import *
import scipy.signal as signal
from scipy.optimize import curve_fit

def CSfilter(x,cutoff, tbs_):
    fs=1./tbs_
    nyq = 0.5 * fs
    filterCutOff=cutoff/nyq
    b, a = signal.butter(1, filterCutOff, 'low')
    return signal.lfilter(b, a, x)

if __name__ == "__main__":
    shotnr=42952;

    prim, times, tbs = getSignal(ch_prim, shotnr )
    hor, times, tbs = getSignal(ch_hor, shotnr )
    vert, times, tbs = getSignal(ch_vert, shotnr )

    #slicing
    slice_start=np.where(times==100000)[0][0]
    slice_end=np.where(times==200000)[0][0]

    times, data=getMirnovs(shotnr,mirnv_corr,True)

    def exponential2(x, a, b, c):
        return a*np.exp(-b*x)+c

    guess = [-7.35300726e-05,2.82177103e-05,4.56193437e-06]

    plt.figure()
    plt.plot(times, vert)
    plt.figure()
    coilNr=0
    for coil in data:
        coilNr+=1
        popt, pcov = curve_fit(exponential2, times[slice_start:slice_end], coil[slice_start:slice_end], p0=guess,  maxfev=50000)
        ax = plt.subplot(3, 4, coilNr)
        plt.plot(times, coil)
        plt.plot(times[slice_start:slice_end], exponential2(times[slice_start:slice_end], *popt), label="fit" )
        guess=popt
        print "MIRNOV "+str(coilNr)+ " lambda: "+str(popt[1]*1e3)+" ms^-1"
        print popt
        filtered=CSfilter(vert,popt[1]/2./np.pi, tbs)
        if popt[0] < 0:
            scale=max(exponential2(times[slice_start:slice_end], *popt))/max(filtered)
        else:
            scale=min(exponential2(times[slice_start:slice_end], *popt))/max(filtered)
            plt.plot(times, filtered*scale)
    plt.show()
