#!/usr/bin/env python3
#fft.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def fft_data(filename,n):
    #input parameters
    D = 0.0508
    Dj = D
    Uj = 762.389735
    # read data
    indata = np.loadtxt(filename+'/'+filename+'_00.dat',usecols=(0,1))
    length=indata.shape[0]
    time = indata[:,0]
    dt = time[1] - time[0]

    # calculate data length
    Fs=1./dt
    df = Fs/length
    St = df * np.arange(length//2+1)*Dj/Uj #half spectrum St array
    TotalPSD = np.zeros(length//2+1) #half spectrum PSD
    for i in range(n):
        pressure=np.loadtxt(filename+'/'+filename+'_{:02d}.dat'.format(i),usecols=(0,1))[:,1]
        pressure=pressure-np.mean(pressure)
        pressure = (fftpack.fft(pressure*np.hanning(length)) *
                        np.sqrt(dt/np.linalg.norm(np.hanning(length), 2)**2))[0:length//2+1]
        PSD = (abs(pressure) ** 2)  # Calculate PSD/Hz
        PSD[1:length-(length//2+1)+1]=2* PSD[1:length-(length//2+1)+1]                                                 
        TotalPSD += PSD
    TotalPSD = TotalPSD/n
    SPL = 10*np.log10(TotalPSD/(Dj/Uj)/4e-10)#-10*np.log10(4)
    return St,SPL

St,SPL = fft_data('smc_35_az_110',5) 
indata2 = np.loadtxt('D210_1557.dat',usecols=(0,3),skiprows=8193*(23-19),max_rows=8193)
St_exp = indata2[:,0]
SPL_exp = indata2[:,1]
fig = plt.figure()
plt.semilogx(St_exp[17:7374],SPL_exp[17:7374],'k',label='Experimental')
plt.semilogx(St[1:-1],SPL[1:-1],'r',label='FWH Prediction')
plt.xlim(left=0.014,right=3)
plt.ylim(bottom=80,top=155)
plt.legend( loc='best', numpoints = 1 )
plt.xlabel('St')
plt.ylabel('PSD, dB/St')
plt.grid()
plt.title('$\\theta = {}^\circ$'.format(110))
plt.savefig('{0:d}_deg'.format(110)+'.png')
plt.close()
# plt.show()
