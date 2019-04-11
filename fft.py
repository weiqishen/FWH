#!/usr/bin/env python3
#fft.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def fft_data(filename,filename2,id):
    #input parameters
    D = 0.0508
    #M_j = 1.46969575
    #gamma = 1.4
    #M_d = 1
    Dj = D
    # Dj = D*((1 + (gamma-1)*M_j**2/2)/(1 + (gamma-1)*M_d**2/2))**((gamma+1)/(4*gamma-4))*(M_d/M_j)**0.5
    Uj = 762.389735
    overlap = 0 #number of overlaping
    length=8000 #number of points per segment
    # read data
    indata = np.loadtxt(filename,usecols=(0,1),skiprows=0)
    time = indata[:,0]
    pressure = indata[:,1]
    dt = time[1] - time[0]
    N_tot = len(time)
    # calculate data length
    Fs=1./dt
    df = Fs/length
    St = df * np.arange(length//2+1)*Dj/Uj #half spectrum St array
    TotalPSD = np.zeros(length//2+1) #half spectrum PSD
    nseg = int((N_tot-overlap)/(length-overlap))
    print('Number of segments: {0:d}'.format(nseg))
    for i in range(nseg):
        seg_pressure = pressure[i*(length-overlap):i*(length-overlap) + length]
        seg_pressure=seg_pressure-np.mean(seg_pressure)
        seg_pressure = (fftpack.fft(seg_pressure*np.hanning(length))*np.sqrt(8./3*dt/length))[0:length//2+1] 
        PSD = (abs(seg_pressure) ** 2)  # Calculate PSD/Hz
        PSD[1:length-(length//2+1)+1]=2* PSD[1:length-(length//2+1)+1]                                                 
        TotalPSD += PSD
    TotalPSD = TotalPSD/nseg
    SPL = 10*np.log10(TotalPSD/(Dj/Uj)/4e-10)#-10*np.log10(4)
    
    indata2 = np.loadtxt(filename2,usecols=(0,3),skiprows=8193*(23-id),max_rows=8193)
    time_a = indata2[:,0]
    pressure_a = indata2[:,1]
    
    return time_a,pressure_a,St,SPL,time,pressure

for i in range(0,24):
    time_a,pressure_a,St,SPL,time,pressure = fft_data('smc_24/smc_24_{0:02d}.dat'.format(i),'D210_1557.dat',i)

    fig = plt.figure()
    plt.semilogx(time_a,pressure_a,'k',label='Experimental')
    plt.semilogx(St[1:-1],SPL[1:-1],'r',label='FWH Prediction')
    plt.xlim(left=0.01,right=3)
    plt.ylim(bottom=80,top=150)

    plt.legend( loc='best', numpoints = 1 )
    plt.xlabel('Strouhal Number')
    plt.ylabel('PSD / St')
    plt.grid()
    plt.title('$\\theta = {}^\circ$'.format(15+5*i))
    plt.savefig('non_avg_result/{}_deg'.format(15+5*i)+'.png')
    plt.close()
# plt.show()