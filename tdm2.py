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

# Tracer le graphique
plt.figure(0)
plt.plot(temps, u)
plt.plot(temps, [u_moy for i in range(len(temps))], label='moyenne', color='red')
plt.plot(temps, [u_moy+u_ecarttype for i in range(len(temps))], label='ecart-type',linestyle='--', color='orange')
plt.plot(temps, [u_moy-u_ecarttype for i in range(len(temps))],linestyle='--', color='orange')
plt.xlabel("Temps(s)")
plt.ylabel("Vitesse U(cm/s)")
plt.legend()

plt.figure(1)
plt.plot(temps, v)
plt.plot(temps, [v_moy for i in range(len(temps))], label='moyenne', color='red')
plt.plot(temps, [v_moy+u_ecarttype for i in range(len(temps))], label='ecart-type',linestyle='--', color='orange')
plt.plot(temps, [v_moy-u_ecarttype for i in range(len(temps))],linestyle='--', color='orange')
plt.xlabel("Temps(s)")
plt.ylabel("Vitesse V(cm/s)")
plt.legend()

plt.figure(2)
plt.plot(temps, w)
plt.plot(temps, [w_moy for i in range(len(temps))], label='moyenne', color='red')
plt.plot(temps, [w_moy+u_ecarttype for i in range(len(temps))], label='ecart-type',linestyle='--', color='orange')
plt.plot(temps, [w_moy-u_ecarttype for i in range(len(temps))],linestyle='--', color='orange')
plt.xlabel("Temps(s)")
plt.ylabel("Vitesse W(cm/s)")
plt.grid()


# Densité de probabilité
plt.figure(4)
plt.hist(u[1000:10000], bins=100, density=True, label='u')
plt.xlabel('vitesse '+ r'$(cm.s^{-1}$)')
plt.ylabel('pdf')
plt.legend()

# Densité de probabilité de la loi normale 
plt.figure(5)
plt.hist(u[1000:10000], bins=100, density=True, cumulative=True, label='u')
plt.xlabel('vitesse '+ r'$(cm.s^{-1}$)')
plt.ylabel('pdf')
plt.legend()


# Distribution de vitesse
def gaussian_density(x):
    return 1/(np.sqrt(2*np.pi))*np.exp(-0.5*x**2)
liste_x = np.linspace(min(u[1000:10000]), max(u[1000:10000]), 100)

#Densité de probabilité
def repartition_speed(u,precision=0.01):
    vitesse = np.sort(u)
    n = len(vitesse)
    v_max = vitesse[n-1]
    v_min = vitesse[0]
    repartition = np.linspace(v_min,v_max,int((v_max-v_min)/precision))
    return repartition


def repartition_vitesse(x,P):
    min_x = min(x)
    delta_x = (max(x) - min(x))/len(x)
    for i in range(len(x)):
        l = (x[i] - min_x)/delta_x
        l = min(float(l), len(x)-1) + 1
        P[l] = P[l] + 1
    return P

def pdf_speed(u,precision=0.01):
    vitesse = np.sort(u)
    n = len(vitesse)
    v_max = vitesse[n-1]
    v_min = vitesse[0]
    repartition = np.linspace(v_min,v_max,int((v_max-v_min)/precision))
    pdf = np.zeros(len(repartition))
    for vit in vitesse:
        for i in range(len(repartition)):
            if vit > repartition[i] and vit <= repartition[i+1]:
                pdf[i] += 1
                break
    proba = pdf / (np.sum(pdf)) 
    return proba

#calculer la moyenne et la variance à partir de la densité de probabilité
def moyenne(repartition, proba):
    moyenne = 0
    for i in range(len(repartition)):
        moyenne += repartition[i]*proba[i]
    return moyenne

def ecart_type(repartition, proba):
    u_moy_densite = moyenne(repartition, proba)
    variance = 0
    for i in range(len(repartition)):
        variance += (repartition[i]-u_moy_densite)**2*proba[i]
    return np.sqrt(variance)

repartition_u = repartition_speed(u)
#proba_u = pdf_speed(u)
delta_x_u= repartition_speed(u)[1]-repartition_speed(u)[0]
print("densité de probabilité de u = ", np.sum(pdf_speed(u))*(delta_x_u))
print("probabilite", np.sum(pdf_speed(u)))

u_moy_densite = moyenne(repartition_speed(u),  pdf_speed(u))
u_ecarttype_densite = ecart_type(repartition_speed(u),  pdf_speed(u))
print("u_moy_proba = ", u_moy_densite)
print("u_ecarttype_proba = ", ecart_type(repartition_speed(u),  pdf_speed(u)))


repartition_v = repartition_speed(v)
#proba_u = pdf_speed(u)
delta_x_v= repartition_speed(v)[1]-repartition_speed(v)[0]
print("densité de probabilité de v = ", np.sum(pdf_speed(v))*(delta_x_v))
print("probabilite", np.sum(pdf_speed(v)))

v_moy_densite = moyenne(repartition_speed(v),  pdf_speed(v))
v_ecarttype_densite = ecart_type(repartition_speed(v),  pdf_speed(v))
print("v_moy_proba = ", v_moy_densite)
print("v_ecarttype_proba = ", ecart_type(repartition_speed(v),  pdf_speed(v)))



repartition_w = repartition_speed(w)
#proba_u = pdf_speed(u)
delta_x_w= repartition_speed(w)[1]-repartition_speed(w)[0]
print("densité de probabilité de w = ", np.sum(pdf_speed(w))*(delta_x_w))
print("probabilite", np.sum(pdf_speed(w)))

w_moy_densite = moyenne(repartition_speed(w),  pdf_speed(w))
w_ecarttype_densite = ecart_type(repartition_speed(w),  pdf_speed(w))
print("w_moy_proba = ", w_moy_densite)
print("w_ecarttype_proba = ", ecart_type(repartition_speed(w),  pdf_speed(w)))


#Comparer la répartition de probabilité avec la loi normale
def loi_normale(x, mean, std):
    return 1/(std*np.sqrt(2*np.pi))*np.exp(-0.5*(x-mean)**2/std**2)

plt.figure(6)
plt.plot(repartition_u, pdf_speed(u), 'r', label='u '+ r'$(cm.s^{-1}$)')
plt.xlabel('vitesse '+ r'$(cm.s^{-1}$)')
plt.ylabel('pdf')
plt.legend()

#Dérivée temporelle du signal de vitesse
def derivation(u,temps) : 
    dt = np.diff(temps)
    du = np.diff(u)
    dt = np.concatenate(([dt[0]], dt))
    du = np.concatenate(([du[0]], du))
    return du/dt

u_derivate = derivation(u,temps)
u_derivate_moy = np.mean(u_derivate)
print("Dérivée u moyenne = ", u_derivate_moy)

v_derivate = derivation(v,temps)
v_derivate_moy = np.mean(v_derivate)
print("Dérivée v moyenne = ", v_derivate_moy)

w_derivate = derivation(w,temps)
w_derivate_moy = np.mean(w_derivate)
print("Dérivée w moyenne = ", w_derivate_moy)


plt.figure(7)
plt.plot(temps, u_derivate)
plt.plot(temps, [u_derivate_moy for i in range(len(temps))], label='moyenne', color='red')
plt.xlabel("Temps (s)")
plt.ylabel("du/dt (cm/s²)")
plt.legend()

plt.figure(8)
plt.plot(temps, v_derivate)
plt.plot(temps, [v_derivate_moy for i in range(len(temps))], label='moyenne', color='red')
plt.xlabel("Temps (s)")
plt.ylabel("dv/dt (cm/s²)")
plt.legend()

plt.figure(9)
plt.plot(temps, w_derivate)
plt.plot(temps, [w_derivate_moy for i in range(len(temps))], label='moyenne', color='red')
plt.xlabel("Temps (s)")
plt.ylabel("dw/dt (cm/s²)")
plt.legend()

plt.figure(10)
plt.hist(u_derivate[1000:10000], bins=100, density=True, label='du/dt')
plt.xlabel('dérivée_vitesse '+ r'$(cm.s^{-2}$)')
plt.ylabel('pdf')
plt.legend()
plt.show()
