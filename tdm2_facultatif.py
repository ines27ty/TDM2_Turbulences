import numpy as np
import matplotlib.pyplot as plt
from random import *
from math import * 


# Calcul de la transformée de Fourier
def fourier(u_fluctuations, temps):
    Ntot = len(u_fluctuations)  
    dt = temps[1]-temps[0]
    freq = np.fft.fftfreq(Ntot, dt)
    freq = np.fft.fftshift(freq)
    u_fluctuations_fft = np.fft.fft(u_fluctuations)
    u_fluctuations_fft = np.fft.fftshift(u_fluctuations_fft)
    return freq, u_fluctuations_fft


def fft(u_fluctuations,temps) : 
    Ntot = len(u_fluctuations) 
    N = Ntot//10
    T = temps[-1]-temps[0]
    u_fft = np.zeros(N,dtype=complex)
    freq = np.fft.fftfreq(N, T/Ntot)
    freq = np.fft.fftshift(freq)
    for i in range(N) : 
        for j in range(N) : 
            u_fft[j] += u_fluctuations[i]*np.exp(-2j*np.pi*freq[i]*temps[j])
    u_fft = np.fft.fftshift(u_fft)
    return u_fft,freq

#Densité spectrale de puissance
def densite_spec(u_fluctuations,time):
    Ntot = len(u_fluctuations)  
    dt = time[1]-time[0]
    freq = np.fft.fftfreq(Ntot, dt)
    freq = np.fft.fftshift(freq)
    u_fluctuations_fft = np.fft.fft(u_fluctuations)
    u_fluctuations_fft = np.fft.fftshift(u_fluctuations_fft)
    DSP = (abs(u_fluctuations_fft))**2
    return freq, DSP

def densite_spectrale(u_fluctuations,temps,freq) :
    Ntot = len(u_fluctuations) 
    N = Ntot//10
    dt = temps[1]-temps[0]
    T = temps[-1]-temps[0]
    sum = np.zeros(N,dtype=complex)
    for i in range(N) :
        sum[i] += u_fluctuations[i]*np.exp(-2j*np.pi*freq[i]*temps[i])
    return sum*dt**2/T

#Bruit blanc
bruit_blanc = np.random.normal(0, 1, 10000)
time = np.linspace(0, 100, 10000)
bruit_blanc_fft = fourier(bruit_blanc, time)[1]
freq_bruit_blanc = fourier(bruit_blanc, time)[0]

bruit_blanc_fft2= fft(bruit_blanc,time)[0]
freq_bruit_blanc2 = fft(bruit_blanc,time)[1]

densite_spec_bruit_blanc = densite_spec(bruit_blanc,time)[1]
freq_densite_spec_bruit_blanc = densite_spec(bruit_blanc,time)[0]

densite_spectrale_bruit_blanc = densite_spectrale(bruit_blanc,time,freq_densite_spec_bruit_blanc)

#FFT du bruit blanc
plt.figure(0)
plt.plot(freq_bruit_blanc, abs(bruit_blanc_fft))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT du bruit blanc"])

plt.figure(1)
plt.plot(freq_bruit_blanc2, abs(bruit_blanc_fft2))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT du bruit blanc"])

#Densité spectrale de puissance du bruit blanc
plt.figure(2)
plt.plot(freq_bruit_blanc, abs(densite_spec_bruit_blanc))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP du bruit blanc"])

#Densité spectrale de puissance du bruit blanc
plt.figure(3)
plt.plot(freq_bruit_blanc2 , abs(densite_spectrale_bruit_blanc))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP du bruit blanc"])
plt.show()


