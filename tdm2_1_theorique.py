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
    C /= var
    tau=time[0 : N]
    return tau, C


#Calcul de l'échelle intégrale
def echelle_int(tau, C):
    N=len(tau)
    EI = 0
    for i in range(N-1):
        EI += C[i]*(tau[i+1]-tau[i])
    return EI

##Application aux données
#Bruit blanc
bruit_blanc = np.random.normal(0, 1, 10000)
time = np.linspace(0, 100, 10000)

#Comparaison avec un bruit blanc d'intensité double
bruit_blanc2 = np.random.normal(0, 2, 10000)


tau_bruit_blanc, autocorr_bruitblanc = autocorr(bruit_blanc,time)
echelle_int_bruitblanc = echelle_int(tau_bruit_blanc,autocorr_bruitblanc)
print("Echelle intégrale du bruit blanc =", echelle_int_bruitblanc)

tau_bruit_blanc2, autocorr_bruitblanc2 = autocorr(bruit_blanc2,time)
echelle_int_bruitblanc2 = echelle_int(tau_bruit_blanc2,autocorr_bruitblanc2)
print("Echelle intégrale du bruit blanc 2 =", echelle_int_bruitblanc2)

bruit_blanc05 = np.random.normal(0, 0.5, 10000)

tau_bruit_blanc05, autocorr_bruitblanc05 = autocorr(bruit_blanc05,time)
echelle_int_bruitblanc05 = echelle_int(tau_bruit_blanc05,autocorr_bruitblanc05)
print("Echelle intégrale du bruit blanc 0.5 =", echelle_int_bruitblanc05)


plt.figure(0)
plt.plot(time, bruit_blanc, zorder=2)
plt.plot(time, bruit_blanc2,  color='orange', zorder=1)
plt.plot(time, bruit_blanc05, color = 'green', zorder=3)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['Bruit blanc', 'Bruit blanc intensité double', 'Bruit blanc intensité moitié'])


plt.figure(1)
plt.plot(tau_bruit_blanc, autocorr_bruitblanc)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')


plt.figure(2)
plt.plot(tau_bruit_blanc2, autocorr_bruitblanc2, color='orange', zorder=3)
plt.plot(tau_bruit_blanc, autocorr_bruitblanc05, color = 'green',zorder=2)
plt.plot(tau_bruit_blanc, autocorr_bruitblanc, zorder=1)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['Bruit blanc intensité double', 'Bruit blanc intensité moitié','Bruit blanc'])


#Sinusoide
f0=0.1
sinusoide = np.sin(2*np.pi*time*f0)

plt.figure(3)
plt.plot(time, sinusoide)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')

tau_sinus, autocorr_sinus = autocorr(sinusoide,time)
echelle_int_sinus = echelle_int(tau_sinus,autocorr_sinus)
print("Echelle intégrale de la sinusoide =", echelle_int_sinus)
autocorr_sinus_th = np.cos(2*np.pi*tau_sinus*f0)

plt.figure(4)
plt.plot(tau_sinus, autocorr_sinus)
plt.plot(tau_sinus, autocorr_sinus_th,'--')
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['Sinusoide', 'Theorie'])

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

plt.figure(5) 
plt.plot(temps_langevin, signal_langevin)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')

tau_langevin, autocorr_langevin = autocorr(signal_langevin,temps_langevin)
echelle_int_langevin = echelle_int(tau_langevin, autocorr_langevin)
print("Echelle intégrale de langevin =", echelle_int_langevin)

#Comparaison avec un processus de Langevin d'intensité double
temps_langevin2 = Y[:,0]
signal_langevin2 = Y[:,1]*2  
tau_langevin2, autocorr_langevin2 = autocorr(signal_langevin2,temps_langevin2)
echelle_int_langevin2 = echelle_int(tau_langevin2, autocorr_langevin2)
print("Echelle intégrale de langevin 2 =", echelle_int_langevin2)


#Comparaison avec un processus de Langevin d'intensité moitié
temps_langevin05 = Y[:,0]
signal_langevin05 = Y[:,1]*0.5  
tau_langevin05, autocorr_langevin05 = autocorr(signal_langevin05,temps_langevin05)
echelle_int_langevin05 = echelle_int(tau_langevin05, autocorr_langevin05)
print("Echelle intégrale de langevin 0.5 =", echelle_int_langevin05)


plt.figure(6)
plt.plot(temps_langevin, signal_langevin,linewidth=2)
plt.plot(temps_langevin, signal_langevin2)
plt.plot(temps_langevin, signal_langevin05)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['Langevin', 'Langevin intensité double', 'Langevin intensité moitié'])

plt.figure(7) 
plt.plot(tau_langevin, autocorr_langevin,linewidth=2)
plt.plot(tau_langevin2, autocorr_langevin2,linestyle='--')
plt.plot(tau_langevin05, autocorr_langevin05,linestyle='-.')
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['Langevin', 'Langevin intensité double','langenvin intensité moitié'])

#Comparaison avec un processus de Langevin d'une autre moyenne
temps_langevin3 = Y[:,0]
signal_langevin3 = Y[:,1]+3
tau_langevin3, autocorr_langevin3 = autocorr(signal_langevin3,temps_langevin3)
echelle_int_langevin3 = echelle_int(tau_langevin3, autocorr_langevin3)
print("Echelle intégrale de langevin +3=", echelle_int_langevin3)

temps_langevinm3 = Y[:,0]
signal_langevinm3 = Y[:,1]-3
tau_langevinm3, autocorr_langevinm3 = autocorr(signal_langevinm3,temps_langevinm3)
echelle_int_langevinm3 = echelle_int(tau_langevinm3, autocorr_langevinm3)
print("Echelle intégrale de langevin -3=", echelle_int_langevinm3)

plt.figure(8)
plt.plot(temps_langevin, signal_langevin,linewidth=2)
plt.plot(temps_langevin3, signal_langevin3)
plt.plot(temps_langevinm3, signal_langevinm3)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['Langevin', 'Langevin de moyenne +3', 'Langevin de moyenne -3'])

plt.figure(9) 
plt.plot(tau_langevin, autocorr_langevin,linewidth=2)
plt.plot(tau_langevin3, autocorr_langevin3,linestyle='--')
plt.plot(tau_langevinm3, autocorr_langevinm3,linestyle='--')
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['Langevin', 'Langevin de moyenne +3','Langevin de moyenne -3'])

## Influence du temps caractéristique dans le processus de Langevin


#Processus de Langevin en fonction du temps caractéristique
def langevin_temps(T):    
    dt = 0.01               #the time step of the signal
    Nit = 10000             #the number of time step of the signal 
    sigma = 2               #the variance of the langevin signal
    #T = 1.                  #tha characteristic time of the signal
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

Y_t1 = langevin_temps(1)
temps_langevin_t1 = Y_t1[:,0]
signal_langevin_t1 = Y_t1[:,1]

Y_t10 = langevin_temps(10)
temps_langevin_t10 = Y_t10[:,0]
signal_langevin_t10 = Y_t10[:,1]

Y_t100 = langevin_temps(100)
temps_langevin_t100 = Y_t100[:,0]
signal_langevin_t100 = Y_t100[:,1]

Y_t1000 = langevin_temps(1000)
temps_langevin_t1000 = Y_t1000[:,0]
signal_langevin_t1000 = Y_t1000[:,1]


Y_t10000 = langevin_temps(10000)
temps_langevin_t10000 = Y_t1000[:,0]
signal_langevin_t10000 = Y_t1000[:,1]

tau_langevin_t1, autocorr_langevin_t1 = autocorr(signal_langevin_t1,temps_langevin_t1)
echelle_int_langevin_t1 = echelle_int(tau_langevin_t1, autocorr_langevin_t1)
print("Echelle intégrale de langevin T=1=", echelle_int_langevin_t1)

tau_langevin_t10, autocorr_langevin_t10 = autocorr(signal_langevin_t10,temps_langevin_t10)
echelle_int_langevin_t10 = echelle_int(tau_langevin_t10, autocorr_langevin_t10)
print("Echelle intégrale de langevin T=10=", echelle_int_langevin_t10)

tau_langevin_t100, autocorr_langevin_t100 = autocorr(signal_langevin_t100,temps_langevin_t100)
echelle_int_langevin_t100 = echelle_int(tau_langevin_t100, autocorr_langevin_t100)
print("Echelle intégrale de langevin T=100=", echelle_int_langevin_t100)

tau_langevin_t1000, autocorr_langevin_t1000 = autocorr(signal_langevin_t1000,temps_langevin_t1000)
echelle_int_langevin_t1000 = echelle_int(tau_langevin_t1000, autocorr_langevin_t1000)
print("Echelle intégrale de langevin T=1000=", echelle_int_langevin_t1000)

tau_langevin_t10000, autocorr_langevin_t10000 = autocorr(signal_langevin_t10000,temps_langevin_t10000)
echelle_int_langevin_t10000 = echelle_int(tau_langevin_t10000, autocorr_langevin_t10000)
print("Echelle intégrale de langevin T=10000=", echelle_int_langevin_t10000)

plt.figure(10)
plt.plot(temps_langevin_t1, signal_langevin_t1)
plt.plot(temps_langevin_t10, signal_langevin_t10)
plt.plot(temps_langevin_t100, signal_langevin_t100)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['T=1','T=10','T=100'])

plt.figure(11)
plt.plot(temps_langevin_t1, signal_langevin_t1)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['T=1'])

plt.figure(12)
plt.plot(temps_langevin_t10, signal_langevin_t10)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['T=10'])

plt.figure(13)
plt.plot(temps_langevin_t100, signal_langevin_t100)
plt.xlabel('Temps (s)')
plt.ylabel('Amplitude')
plt.legend(['T=100'])

plt.figure(14) 
plt.plot(tau_langevin_t1, autocorr_langevin_t1)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['T=1'])

plt.figure(15) 
plt.plot(tau_langevin_t10, autocorr_langevin_t10)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['T=10'])

plt.figure(16) 
plt.plot(tau_langevin_t100, autocorr_langevin_t100)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['T=100'])


plt.figure(17) 
plt.plot(tau_langevin_t1, autocorr_langevin_t1)
plt.plot(tau_langevin_t10, autocorr_langevin_t10)
plt.plot(tau_langevin_t100, autocorr_langevin_t100)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.legend(['T=1','T=10','T=100'])


#Tracer l'échelle intégrale en fonction du temps caractéristique T


X=np.array([1,10,100,1000])
Y=np.array([echelle_int_langevin_t1,echelle_int_langevin_t10,echelle_int_langevin_t100,echelle_int_langevin_t1000])
plt.figure(18)
plt.plot(X,Y,'o')
plt.xlabel('Temps caractéristique')
plt.ylabel('Echelle intégrale')


plt.show()








