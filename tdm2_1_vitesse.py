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

#Calcul de l'échelle intégrale
def echelle_int(tau, C):
    EI = 0
    EI = np.sum(C) * (tau[1] - tau[0])
    return EI

## Auto-correlation et échelle intégrale du signal de vitesse

plt.figure(0)
plt.plot(temps, u)
plt.xlabel("Temps(s)")
plt.ylabel("Vitesse U(cm/s)")

tau_vitesse_u, autocorr_vitesse_u = autocorr(u_fluctuations,temps)
echelle_int_vitesse_u = echelle_int(tau_vitesse_u, autocorr_vitesse_u)
print("Echelle intégrale de la vitesse u =", echelle_int_vitesse_u)

plt.figure(1)
plt.plot(tau_vitesse_u[:-1], autocorr_vitesse_u[:-1])
plt.xlabel('Temps (s)')
plt.ylabel('Auto-correlation')
plt.show()
