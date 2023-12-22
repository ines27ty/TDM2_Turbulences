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
    return freq, DSP/len(u_fluctuations)

def densite_spectrale(u_fluctuations,temps,freq) :
    Ntot = len(u_fluctuations) 
    N = Ntot//10
    dt = temps[1]-temps[0]
    T = temps[-1]-temps[0]
    sum = np.zeros(N,dtype=complex)
    for i in range(N) :
        sum[i] += u_fluctuations[i]*np.exp(-2j*np.pi*freq[i]*temps[i])
    return abs(sum)**2*dt**2/T

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

#Sinusoide
f0=0.1
sinusoide = np.sin(2*np.pi*time*f0)

sinusoide_fft = fourier(sinusoide, time)[1]
freq_sinusoide = fourier(sinusoide, time)[0]

sinusoide_fft2= fft(sinusoide,time)[0]
freq_sinusoide2 = fft(sinusoide,time)[1]

densite_spec_sinusoide = densite_spec(sinusoide,time)[1]
freq_densite_spec_sinusoide = densite_spec(sinusoide,time)[0]

densite_spectrale_sinusoide = densite_spectrale(sinusoide,time,freq_densite_spec_sinusoide)

#FFT de la sinusoide
plt.figure(4)
plt.plot(freq_sinusoide, abs(sinusoide_fft))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT de la sinusoide"])

plt.figure(5)
plt.plot(freq_sinusoide2, abs(sinusoide_fft2))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT de la sinusoide"])

#Densité spectrale de puissance de la sinusoide
plt.figure(6)
plt.plot(freq_sinusoide, abs(densite_spec_sinusoide))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP de la sinusoide"])

#Densité spectrale de puissance de la sinusoide
plt.figure(7)
plt.plot(freq_sinusoide2 , abs(densite_spectrale_sinusoide))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP de la sinusoide"])


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

signal_langevin_fft = fourier(signal_langevin, temps_langevin)[1]
freq_signal_langevin = fourier(signal_langevin, temps_langevin)[0]

signal_langevin_fft2= fft(signal_langevin,temps_langevin)[0]
freq_signal_langevin2 = fft(signal_langevin,temps_langevin)[1]

densite_spec_signal_langevin = densite_spec(signal_langevin,temps_langevin)[1]
freq_densite_spec_signal_langevin = densite_spec(signal_langevin,temps_langevin)[0]

densite_spectrale_signal_langevin = densite_spectrale(signal_langevin,temps_langevin,freq_densite_spec_signal_langevin)

#FFT du processus de Langevin
plt.figure(8)
plt.plot(freq_signal_langevin, abs(signal_langevin_fft))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT du processus de Langevin"])

plt.figure(9)
plt.plot(freq_signal_langevin2, abs(signal_langevin_fft2))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT du processus de Langevin"])

#Densité spectrale de puissance du processus de Langevin
plt.figure(10)
plt.plot(freq_signal_langevin, abs(densite_spec_signal_langevin))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP du processus de Langevin"])

#Densité spectrale de puissance du processus de Langevin
plt.figure(11)
plt.plot(freq_signal_langevin2 , abs(densite_spectrale_signal_langevin))
plt.xlabel("Frequence (Hz)")
plt.ylabel("DSP")
plt.legend(["DSP du processus de Langevin"])

#Verification de la densité spectrale de puissance
def verif(densite,temps) :
    dt = temps[1]-temps[0]
    return 2*np.sum(densite[1:])*dt


print("Verification DSP bruit blanc = ", verif(densite_spec_bruit_blanc,time))
print("Variance bruit blanc = ", np.var(bruit_blanc))

print("Verification DSP sinusoide = ", verif(densite_spec_sinusoide,time))
print("Variance sinusoide = ", np.var(sinusoide))

print("Verification DSP processus de Langevin = ", verif(densite_spec_signal_langevin,temps_langevin))
print("Variance processus de Langevin = ", np.var(signal_langevin))


#Idem avec la vitesse

# Ouvrir le fichier
with open("signal_canal3280_48.txt", 'r') as fichier:
    # Lire les lignes du fichier
    lines = fichier.readlines()

# Initialiser des listes vides pour stocker les données
iteration = []
temps = []
u = []
v = []
w = []

# Parcourir chaque ligne du fichier
for line in lines:
    line = line.strip()         # Supprimer les espaces en début et fin de ligne
    line = line.replace('\t', ',')      # Remplacer les tabulations par des virgules
    values = line.split(',')     # Diviser la ligne en une liste de valeurs
    
    # Convertir chaque valeur en float et ajouter à la liste correspondante
    iteration.append(float(values[0]))
    temps.append(float(values[1]))
    u.append(float(values[2]))
    v.append(float(values[3]))
    w.append(float(values[4]))

#Definition des fluctuations du signal de vitesse
# Moyenne et écart-type des trois composantes du signal
u_moy =np.mean(u)
v_moy = np.mean(v)
w_moy = np.mean(w)

u_ecarttype= np.std(u)
v_ecarttype= np.std(v)
w_ecarttype= np.std(w)

print("u_moy=", u_moy)
print("v_moy=", v_moy)
print("w_moy=", w_moy)

print("u_ecartype=", u_ecarttype)
print("v_ecartype=", v_ecarttype)
print("w_ecartype=", w_ecarttype)

u_fluctuations = [u[i] - u_moy for i in range(len(u))]
v_fluctuations = [v[i] - v_moy for i in range(len(v))]
w_fluctuations = [w[i] - w_moy for i in range(len(w))]

#FFT de la vitesse
u_fft = fourier(u_fluctuations, temps)[1]
freq_u = fourier(u_fluctuations, temps)[0]

u_fft2= fft(u_fluctuations,temps)[0]
freq_u2 = fft(u_fluctuations,temps)[1]

densite_spec_u = densite_spec(u_fluctuations,temps)[1]
freq_densite_spec_u = densite_spec(u_fluctuations,temps)[0]

v_fft = fourier(v_fluctuations, temps)[1]
freq_v = fourier(v_fluctuations, temps)[0]

v_fft2= fft(v_fluctuations,temps)[0]
freq_v2 = fft(v_fluctuations,temps)[1]

densite_spec_v = densite_spec(v_fluctuations,temps)[1]
freq_densite_spec_v = densite_spec(v_fluctuations,temps)[0]

w_fft = fourier(w_fluctuations, temps)[1]
freq_w = fourier(w_fluctuations, temps)[0]

w_fft2= fft(w_fluctuations,temps)[0]
freq_w2 = fft(w_fluctuations,temps)[1]

densite_spec_w = densite_spec(w_fluctuations,temps)[1]
freq_densite_spec_w = densite_spec(w_fluctuations,temps)[0]

densite_spectrale_u = densite_spectrale(u_fluctuations,temps,freq_densite_spec_u)
densite_spectrale_v = densite_spectrale(v_fluctuations,temps,freq_densite_spec_v)
densite_spectrale_w = densite_spectrale(w_fluctuations,temps,freq_densite_spec_w)

#FFT de la vitesse
plt.figure(12)
plt.plot(freq_u, abs(u_fft))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT de la vitesse"])

plt.figure(13)
plt.plot(freq_u2, abs(u_fft2))
plt.xlabel("Frequence (Hz)")
plt.ylabel("Amplitude")
plt.legend(["FFT de la vitesse"])

print("Verification DSP vitesse u = ", verif(densite_spec_u,temps))
print("Verification DSP vitesse u = ", verif(densite_spectrale_u,temps))
print("Vérification variance vitesse u = ", np.var(u_fluctuations))

print("Verification DSP vitesse v = ", verif(densite_spec_v,temps))
print("Verification DSP vitesse v = ", verif(densite_spectrale_v,temps))
print("Variance vitesse v = ", np.var(v_fluctuations))

print("Verification DSP vitesse w = ", verif(densite_spec_w,temps))
print("Verification DSP vitesse w = ", verif(densite_spectrale_w,temps))
print("Variance vitesse w = ", np.var(w_fluctuations))

#plt.show()
