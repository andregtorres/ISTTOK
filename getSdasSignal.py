from sdas.core.SDAStime import Date, Time, TimeStamp
import numpy as np
from StartSdas import StartSdas

client = StartSdas()

#signal + time
def getSignal(channelID, shotNr, gain=1.):
    signalStructArray=client.getData(channelID,'0x0000', shotNr)
    signalStruct=signalStructArray[0]
    signal=signalStruct.getData()*gain
    tstart = signalStruct.getTStart()
    tend = signalStruct.getTEnd()
    #Calculate the time between samples
    tbs = (tend.getTimeInMicros() - tstart.getTimeInMicros())/(len(signal)*1.0)
    #Get the events  associated with this data
    events = signalStruct.get('events')
    tevent = TimeStamp(tstamp=events[0].get('tstamp'))
    #The delay of the start time relative to the event time
    delay = tstart.getTimeInMicros() - tevent.getTimeInMicros()
    #Finally create the time array
    times = np.linspace(delay,delay+tbs*(len(signal)-1),len(signal))
    return (signal, times, tbs)
