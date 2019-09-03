import numpy as np
import linecache
import matplotlib
matplotlib.use('Agg')
import pywt
import astropy.io.fits as pyfits
import os
import datetime
import ephem
import time
import sys
from decimal import Decimal
from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

secperday = 3600 * 24
path=sys.argv[1]
#day=sys.argv[3]
yemo=sys.argv[2]
os.environ['path']=str(path)
os.environ['yemo']=str(yemo)
#os.environ['day']=str(day)
#os.system('mkdir -p /data1/home/yuanmao/wideband_re/$yemo/$day')
#os.system('mkdir -p /data1/home/yuanmao/wideband_re/fits_file/$yemo/$day')
#os.system('mkdir -p /data1/home/yuanmao/wideband_re/png/$yemo')
os.system('mkdir -p /home/yuanmao/M31_rfi/$yemo')
channel_bad=[]
day_ra_dec=[]
day_mjd=[]

##os.system('find /data2/wideband/psr_2017/$path/$day/ -name "*drifting*_0-1GHz_*.fits" -exec basename {} \; >/data1/home/yuanmao/wideband_re/shell/bandpass_${yemo}${day}.list')
os.system('find $path -name "*M01*.fits" -exec basename {} \; >/home/yuanmao/M31_rfi/fitslist/${yemo}.list')
os.system('sort -d /home/yuanmao/M31_rfi/fitslist/${yemo}.list -o /home/yuanmao/M31_rfi/fitslist/${yemo}.list')
#print ('Begin to get bandpass as a 32768*2800 data set')

def smooth(sig,threshold = 3, level=8, wavelet='db8'):
    sigma = sig.std()
    dwtmatr = pywt.wavedec(data=sig, wavelet=wavelet, level=level)
    denoised = dwtmatr[:]
    denoised[1:] = [pywt.threshold(i, value=threshold*sigma, mode='soft') for i in dwtmatr[1:]]
    smoothed_sig = pywt.waverec(denoised, wavelet, mode='sp1')[:sig.size]
    noises = sig - smoothed_sig
    return smoothed_sig, noises
count=0
for name in open('/home/yuanmao/M31_rfi/fitslist/%s.list'%(yemo), 'r'):
    print (name)
    filename = '%s'%(path) + name
    filename=filename.replace("\n","")
    try:
        hdulist = pyfits.open(filename)
        hdu0 = hdulist[0]
        hdu1 = hdulist[1]
        data1 = hdu1.data['data']
        tsamp = hdu1.header['TBIN']
        a,b,c,d,e = data1.shape
        if c > 1:
                 data = data1[:,:,1,:,:].squeeze().reshape((-1,d))
        else:

                 data = data1.squeeze().reshape((-1,d))
        l, m = data.shape
        data = data.reshape(l/8, 8, d).sum(axis=1)
        #data = data[:,400:3600]
        data1 = np.sum(data,axis=0)
        if count==0:
            subintoffset = hdu1.header['NSUBOFFS']   
            samppersubint = int(hdu1.header['NSBLK'])
            tstart = "%.18f" % (Decimal(hdu0.header['STT_IMJD']) + Decimal(hdu0.header['STT_SMJD'] + tsamp * samppersubint * subintoffset )/secperday )
            jd=float(tstart)+2400000.5
            date=ephem.julian_date('1899/12/31 12:00:00')
            djd=jd-date
            str1=ephem.Date(djd).datetime()
            #str2=str1.strftime('%Y-%m-%d %H:%M:%S')
            str2=str1.strftime('%H:%M:%S')
            t_total = tsamp*l
    except:
            print ('Error')
    else:
        basename = name[:name.find(".fits")]
        sig,nos = smooth(data1,level=5)###level can change
        idxarr = np.arange(data1.size)
        y2 = abs((sig-data1)/sig)### shreshold for RFI
        nan=np.isnan(y2)
        inf=np.isinf(y2)
        y2[nan]=0
        y2[y2==1]=0
        idxgood = idxarr[y2<0.05]### good channel
        idxbad = idxarr[y2>=0.05]
        channel_list=np.arange(4096)
        channel_list[idxgood]=0
        channel_list[idxbad]=1
        channel_bad.append(channel_list)
        badchannel = idxarr[idxbad]
        np.savetxt('/home/yuanmao/M31_rfi/%s/%s.txt'%(yemo,basename),badchannel)
        count +=1

count =count
print (count)
gs = gridspec.GridSpec(2, 2,width_ratios=[14,3],height_ratios=[3,10])
gs.update(hspace=0.07, wspace=0.05, top=0.90, bottom=0.10)
fig=plt.figure()
ax1 = fig.add_subplot(gs[0])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])
ax1.set_xticks([])
ax4.set_yticks([])
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height
data=np.array(channel_bad)
print (data.shape)
nos_f=data.sum(axis=1)
nos_f=nos_f/4096.0
nos_t_=data.sum(axis=0)/float(count)
nos_t=nos_t_.T
x_f=np.arange(len(nos_f))
x_t=np.arange(len(nos_t))
norm1 = matplotlib.colors.Normalize(vmin=0, vmax=1)

    #imshow
ax3.set_yticks([1,819.2,1638.4,2457.6,3276,4095])
ax3.set_yticklabels([1000,1100, 1200,1300,1400,1500])
new_test_xt=pd.date_range(start=str2, periods=count,freq='%s s'%(int(t_total)))
ax3.set_xticks(np.linspace(0,len(new_test_xt),5))
t_arr0 = new_test_xt.to_pydatetime()
t_arr1 = np.vectorize(lambda s: s.strftime('%H:%M:%S'))(t_arr0)
t_arr =  list(pd.Series(t_arr1))
x_ticks = np.arange(0,count-1,(count-1)//4)
x_tarr = [t_arr[i] for i in x_ticks]
ax3.set_xticklabels(x_tarr)
ax3.set_ylabel('Frequency(MHz)')
ax3.imshow(data.T,aspect='auto',origin='lower',norm=norm1)
    #top
ax1.plot(x_f,nos_f,linewidth =0.3)
ax1.text(0.5* (left + right),1.1* (top+bottom),'The proportion of RFI',
horizontalalignment='center',
verticalalignment='center', rotation='0',
        transform=ax1.transAxes)
ax1.set_xlim(0,len(nos_f))
    #right
ax4.plot(nos_t,x_t,linewidth =0.3)
ax4.text(1.1 * (left + right),0.5 * (top+bottom),'Apperence rate of RFI',
            horizontalalignment='center',
            verticalalignment='center',
            rotation='270',
            transform=ax4.transAxes)
#ax4.set_xticks([0,0.5,1])
#ax4.set_xticklabels([0,0.5,1])
ax4.set_ylim(0,len(nos_t))
np.savetxt('/home/yuanmao/M31_rfi/%s/a_badchannenl_rate.txt'%(yemo),nos_t)
plt.savefig('/home/yuanmao/M31_rfi/%s/a_rfi_stat.png'%(yemo),dpi=800)

