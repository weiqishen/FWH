#!/usr/bin/env python3
#fft.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def fft_data(filename):
    #input parameters
    D = 0.0508
    #M_j = 1.46969575
    #gamma = 1.4
    #M_d = 1
    Dj = D
    # Dj = D*((1 + (gamma-1)*M_j**2/2)/(1 + (gamma-1)*M_d**2/2))**((gamma+1)/(4*gamma-4))*(M_d/M_j)**0.5
    Uj = 762.389735
    overlap = 0 #number of overlaping
    #length=8000 #number of points per segment
    # read data
    indata = np.loadtxt(filename,usecols=(0,1))
    length=indata.shape[0]
    time = indata[:,0]
    pressure = indata[:,1]
    dt = time[1] - time[0]
    N_tot = len(time)
    # calculate data length
    Fs=1./dt
    df = Fs/length
    St = df * np.arange(length//2+1)*Dj/Uj #half spectrum St array
    TotalPSD = np.zeros(length//2+1) #half spectrum PSD
    nseg = (N_tot-overlap)//(length-overlap)
    print('Number of segments: {0:d}'.format(nseg))
    for i in range(nseg):
        seg_pressure = pressure[i*(length-overlap):i*(length-overlap) + length].copy()
        seg_pressure = seg_pressure-np.mean(seg_pressure)
        seg_pressure = fftpack.fft(seg_pressure*np.hanning(length))[0:length//2+1] #half spectrum
        #seg_pressure = fftpack.fft(seg_pressure)[0:length//2+1]
        PSD = (dt/np.linalg.norm(np.hanning(length), 2)**2)*abs(seg_pressure) ** 2  # Calculate PSD/Hz
       #PSD = dt*8/(3*length)*abs(seg_pressure) ** 2  # Calculate PSD/Hz
        PSD[1:length-(length//2+1)+1]=2* PSD[1:length-(length//2+1)+1]                                                 
        TotalPSD += PSD
    TotalPSD = TotalPSD/nseg
    SPL = 10*np.log10(TotalPSD/(Dj/Uj)/4e-10)#-10*np.log10(4)
    return St,SPL

for i in range(0,24):
    St,SPL = fft_data('smc_35_24/smc_35_24_{0:02d}.dat'.format(i)) 
    indata2 = np.loadtxt('D210_1557.dat',usecols=(0,3),skiprows=8193*(23-i),max_rows=8193)
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
    plt.title('$\\theta = {}^\circ$'.format(15+5*i))
    plt.savefig('avg_35d/{}_deg'.format(15+5*i)+'.png')
    plt.close()
# plt.show()
