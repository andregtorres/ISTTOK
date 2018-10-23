#Code for acessment of a new external coil at equatorial plane
#Andre Torres
#23-10-18
from getMirnov import *
%matplotlib qt4
#44411, 44412, 44413
shotN=44412
plotAll(*getMirnovs(44412,mirnv_corr,True), show=False)
primary, times,tbs=getSignal( 'MARTE_NODE_IVO3.DataCollection.Channel_093', shotN)
vert, times,tbs=getSignal( ch_vert, shotN)
ip, times_ip,tbs=getSignal( 'MARTE_NODE_IVO3.DataCollection.Channel_088', shotN)
plt.plot(times_ip, ip)

times, mirnovs = getMirnovs(shotN,mirnv_corr,False)
times, mirnovs_raw = getMirnovs(shotN,mirnv_raw,False)


plt.figure()
plt.plot(times, mirnovs_raw[6])

client.searchParametersByName("plasma")

#thanks @ https://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin
def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

fig, ax1=plt.subplots()
ax1.axhline(xmin=-10, xmax=10000, color="lightgrey")
ax1.plot(times, mirnovs[6], color="r", label="Vprobe")
ax1.set_ylabel("Integrated coil singal (V.s)", color="r")
ax2=ax1.twinx()
ax2.plot(times, primary, color="k", label="Iprim")
ax2.plot(times_ip, ip, color="b", label="Ip")
ax2.plot(times, vert, color="g", label="Ivert")
ax2.set_ylabel("(kA)", color="k")
align_yaxis(ax1,0,ax2,0)
plt.tight_layout()
plt.legend()
