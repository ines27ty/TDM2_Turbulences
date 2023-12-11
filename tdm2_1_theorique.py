import numpy as np
import matplotlib.pyplot as plt
from random import *
from math import * 

#Calcul du coeff d'auto-correlation
def autocorr(u_fluctuations,time):
    Ntot = len(u_fluctuations)  
    var = np.std(u_fluctuations)**2
    N = Ntot//10
    C = np.zeros(N)
    for j in range(N):
        for i in range(Ntot-j):
            C[j] += u_fluctuations[i] * u_fluctuations[i+j]
        C[j] /= (Ntot-j)
    C[j] /= var
    tau=time[0 : N]
    return tau, C


#Calcul de l'échelle intégrale
def echelle_int(tau, C):
    EI = 0
    EI = np.sum(C) * (tau[1] - tau[0])
    return EI

##Application aux données
#Bruit blanc
bruit_blanc = np.random.normal(0, 1, 10000)
time = np.linspace(0, 100, 10000)

plt.figure(0)
plt.plot(time, bruit_blanc)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')

tau_bruit_blanc, autocorr_bruitblanc = autocorr(bruit_blanc,time)
echelle_int_bruitblanc = echelle_int(tau_bruit_blanc,autocorr_bruitblanc)
print("Echelle intégrale du bruit blanc =", echelle_int_bruitblanc)

plt.figure(1)
plt.plot(tau_bruit_blanc, autocorr_bruitblanc)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')


#Sinusoide
f0=0.2
sinusoide = np.sin(2*np.pi*time*f0)

plt.figure(2)
plt.plot(time, sinusoide)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')

tau_sinus, autocorr_sinus = autocorr(sinusoide,time)
echelle_int_sinus = echelle_int(tau_sinus,autocorr_sinus)
print("Echelle intégrale de la sinusoide =", echelle_int_sinus)

plt.figure(3)
plt.plot(tau_sinus[:-1], autocorr_sinus[:-1])
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')

#Processus de Langevin
def langevin():    
    dt = 0.01               #the time step of the signal
    Nit = 10000             #the number of time step of the signal 
    sigma = 2               #the variance of the langevin signal
    T = 1.                  #tha characteristic time of the signal
    mu = 3.                 #the mean value of the signal
    sigma2 = sigma*sigma  
    k= np.sqrt(2.*dt/T) 

    #the intial condition
    dw = np.random.randn()
    x0 = sigma * np.random.randn()+mu
    t0 = 0.
    
    Y = [t0,x0,dw]
    
    i=0
    x=x0
    t=t0
    while i<Nit:
        i += 1
        dw = k * np.random.randn()
        dx = sigma * dw
        dx += (mu-x)*dt/T
        
        x += dx
        t += dt
        Y = np.vstack( (Y, [t,x,dw]) )  
         
    np.savetxt("langevin.txt",Y) 
    return Y

Y = langevin()
temps_langevin = Y[:,0]
signal_langevin = Y[:,1]
bruit_langevin = Y[:,2]

plt.figure(4) 
plt.plot(temps_langevin, signal_langevin)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')


tau_langevin, autocorr_langevin = autocorr(signal_langevin,temps_langevin)
echelle_int_langevin = echelle_int(tau_langevin, autocorr_langevin)
print("Echelle intégrale de langevin =", echelle_int_langevin)

plt.figure(5)
plt.plot(tau_langevin, autocorr_langevin)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.show()
