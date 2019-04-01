#!/usr/bin/env python3
#fft.py

import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

def fft_data(filename,filename2):
    indata = np.loadtxt(filename,usecols=(0,1))
    time = indata[:,0]
    pressure = indata[:,1]
    
    D = 0.0508
    M_j = 1.46969575
    gamma = 1.4
    M_d = 1
    Dj = D
    # Dj = D*((1 + (gamma-1)*M_j**2/2)/(1 + (gamma-1)*M_d**2/2))**((gamma+1)/(4*gamma-4))*(M_d/M_j)**0.5
    Uj = 762.389735
    nseg = 1
    
    dt = time[1] - time[0]
    
    N = int(len(pressure) / nseg)
    if (nseg == 1):
        overlap = 0
    else:
        overlap = int(N/4)
        
    print (N,overlap)
    length=N+overlap
    Fs=1./dt
    df = Fs/length
    f = df * np.arange(length//2+1) #frequency array
    St = f*Dj/Uj
    TotalPSD = np.zeros(length//2+1)
    
    for i in range(nseg):
        seg_pressure = pressure[i*(N-overlap):i*(N-overlap) + length]
        print(i*(N-overlap), i*(N-overlap) + length)
        seg_pressure = (fftpack.fft(
            seg_pressure*np.hanning(length))*np.sqrt(8./3*dt/length))[0:length//2+1] 
        PSD = (abs(seg_pressure) ** 2)  # Calculate PSD/Hz
        PSD[1:length+1-length//2]=2* PSD[1:length+1-length//2]                                                 
        TotalPSD += PSD
    TotalPSD = TotalPSD/nseg
    SPL = 10*np.log10(TotalPSD/(Dj/Uj)/4e-10)#-10*np.log10(4)
    
    indata2 = np.loadtxt(filename2,usecols=(0,1))
    time_a = indata2[:,0]
    pressure_a = indata2[:,1]
    
    return time_a,pressure_a,St,SPL,time,pressure

time_a,pressure_a,St,SPL,time,pressure = fft_data('100D/microphonePressure1.dat','exp/micPressure_experimental_70.dat')

fig = plt.figure()
plt.semilogx(time_a,pressure_a,'k',label='Experimental')
plt.semilogx(St[1:-1],SPL[1:-1],'r',label='FWH Prediction')
plt.xlim(right=3)
plt.legend( loc='best', numpoints = 1 )
plt.xlabel('Strouhal Number')
plt.ylabel('PSD / St')
plt.title('$\\theta = 70^\circ$')
plt.savefig('70degree'+'.png')
# plt.show()

#fig = plt.figure()
#plt.plot(time,pressure,linewidth=0.75,label='from farassat\'s 2b formulation')
#plt.plot(time2,pressure2,linewidth=0.75,label='from farassat\'s 2b formulation')
#plt.plot(time3,pressure3,linewidth=0.75,label='from farassat\'s 2b formulation')
#plt.legend( loc='best', numpoints = 1 )
#plt.show()
# fig.savefig("data/monopoleSphere.pdf")
