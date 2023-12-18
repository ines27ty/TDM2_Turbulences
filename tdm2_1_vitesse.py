import numpy as np
import matplotlib.pyplot as plt
from random import *
from math import * 

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

#Calcul de l'échelle de temps intégrale
def echelle_int(tau, C):
    N=len(tau)
    EI = 0
    for i in range(N-1):
        EI += C[i]*(tau[i+1]-tau[i])
    return EI


plt.figure(0)
plt.plot(temps, u)
plt.xlabel("Temps(s)")
plt.ylabel("Vitesse U(cm/s)")
plt.legend(["Vitesse U"])

## Auto-correlation et échelle de temps intégrale du signal de vitesse

tau_vitesse_u, autocorr_vitesse_u = autocorr([u[i] - u_moy for i in range(len(u))],temps)
echelle_int_vitesse_u = echelle_int(tau_vitesse_u, autocorr_vitesse_u)
print("Echelle intégrale de la vitesse u =", echelle_int_vitesse_u)

#tau_vitesse_v, autocorr_vitesse_v = autocorr([v[i] - v_moy for i in range(len(v))],temps)
#echelle_int_vitesse_v = echelle_int(tau_vitesse_v, autocorr_vitesse_v)
#print("Echelle intégrale de la vitesse v =", echelle_int_vitesse_v)

#tau_vitesse_w, autocorr_vitesse_w = autocorr([w[i] - w_moy for i in range(len(w))],temps)
#echelle_int_vitesse_w = echelle_int(tau_vitesse_w, autocorr_vitesse_w)
#print("Echelle intégrale de la vitesse w =", echelle_int_vitesse_w)

plt.figure(1)
plt.plot(tau_vitesse_u, autocorr_vitesse_u)
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')

print("Echelle de longueur u_moy =",echelle_int_vitesse_u*u_moy)

#Calcul de l'échelle intégrale de longueur
def echelle_int_long(u_fluctuations,time):
    Ntot = len(u_fluctuations)  
    N = Ntot//10
    C = np.zeros(N)
    for j in range(N):
        for i in range(Ntot-j):
            C[j] += u_fluctuations[i] * u_fluctuations[i+j]
    tau=time[0 : N]
    N=len(tau)
    ELI = 0
    for i in range(N-1):
        ELI += C[i]*(tau[i+1]-tau[i])
    return ELI


print("Echelle de longueur intégrale intégrale =",echelle_int_long(u_fluctuations,temps))
print("Echelle de longueur intégrale =",echelle_int_vitesse_u*np.linalg.norm([u_fluctuations,v_fluctuations,w_fluctuations]))
print("Echelle de longueur intégrale vitesse  =", echelle_int_vitesse_u * np.linalg.norm(u_fluctuations))

#Longueur d'enregistrement
Ntot = len(u_fluctuations)
print("Longueur d'enregistrement =", Ntot)

plt.show()