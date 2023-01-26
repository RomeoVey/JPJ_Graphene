"""
Created on Wed Dec 14 18:02:04 2022

@author: NeoMagik
"""
from scipy import *
from cmath import *
import matplotlib 
import matplotlib.pyplot as plt
from numpy import *
from time import *
#Création des objets et paramètres
t = 1.0  #énergie de saut
e0=0.0 #énergie de site
eps = 0.00001  #élement infinitésimal
N = 1000
psimoins = zeros((2*N+1,2*N+1))  #tableau indicé de 0 à 2N compris
psiplus = zeros((2*N+1,2*N+1))  #tableau indicé de 0 à 2N compris
listea = zeros(N+1)  #tableau indicé de 0 à N compris
listeb = zeros(N+1)  #tableau indicé de 0 à N compris

#Initialisation
psimoins[N,N] = 1.0  #psi 0
psiplus[N,N] = 0.0
a = e0  #a0
b = 0  #b0
listea[0] = a #rentre le a0 dans la liste
listeb[0] = b #rentre le b0 dans la liste


########################################
ti = time()
#Calculs des coefs
for n in range(1,N):#parcours les nombres de 1 à N-1 compris
    aplus=0  #remet à zero la variable locale
    for i in range(N-n,N+n+1):#application du Hamiltonien
        for j in range(N-n,N+n+1):
            psiplus[i,j] = -b*psiplus[i,j]
    
    for i in range(N-n,N+n+1) :#application du Hamiltonien à psi_n
        for j in range(N-n,N+n+1):
            if n%2==0: #((abs(i)%2 == 0) and (abs(j)%2 == 0)) or ((abs(i)%2 == 1) and (abs(j)%2 == 1))
                psiplus[i,j] = psiplus[i,j]+t*(psimoins[i-1,j]+ psimoins[i,j-1]+psimoins[i,j+1])+e0*psimoins[i,j]
            else: psiplus[i,j] = psiplus[i,j]+t*(psimoins[i+1,j]+ psimoins[i,j-1]+psimoins[i,j+1])+e0*psimoins[i,j]
    for i in range(N-n,N+n+1) :  #parcours les nombres de N-n à N+n compris
        for j in range(N-n,N+n+1):
          aplus = aplus + psimoins[i,j]*psiplus[i,j]
    a = aplus
    bplus = 0  #remet à zero la variable locale
    for i in range(N-n,N+n+1) :#parcours les nombres de N-n à N+n compris
        for j in range(N-n,N+n+1):
            psiplus[i,j] = psiplus[i,j]-a*psimoins[i,j]
            bplus = bplus + (psiplus[i,j])**2
    bplus = bplus**(1/2)
    for i in range(N-n,N+n+1) :  #parcours les nombres de N-n à N+n compris
        for j in range(N-n,N+n+1):
            psiplus[i,j] = (1/bplus)*psiplus[i,j]
    temp = psiplus  #pour inverser psiplus et psimoins
    psiplus = psimoins  #pour inverser psiplus et psimoins
    psimoins = temp  #pour inverser psiplus et psimoins
    b = bplus

    listea[n] = a  #modifie la liste petit à petit
    listeb[n] = b

    print("Temps iteration :", round(time()-ti,2),"s\t"+str(n)+"/"+str(N))
    ti = time()
#Calculs des coefs infinis
ainf = 0.0
binf = 0.0        
for k in range(N-10,N):  #parcours les nombres de N-10 à N-1 compris
    ainf = ainf + listea[k]
    binf = binf + listeb[k]
ainf = ainf/10
binf = binf/10
listea[N] = ainf
listeb[N] = binf
#print("a=",listea)
#print("b=",listeb)

#Calcul de la densité sous forme de fonction        
def densite(E):
    G = (E - ainf -1j*sqrt(-(E - ainf)**2 + 4*binf**2))/(2*binf**2)  #calcul de Ginfini
    for i in range(N,0,-1): #parcours les nombres de N à 0 compris
        #print("in densité i=",i,"b=",listeb[i])
        G = 1/(E - listea[i] - G*listeb[i]**2)  #calcul de la fraction continue
    nmu = abs((1/pi)*G.imag)  #calcul de la densite en prenant la valeur absolue
    return nmu

energy = arange(ainf - 2*binf,ainf + 2*binf, 0.01)  #découpe des valeurs d'énergie     
dos = densite(energy)

pre_int = dos*energy

s=0
for i in range(len(energy)//2 + 1): s+=  0.01*pre_int[i]

#print(s)

plt.plot(energy,dos) 
plt.xlabel("Energie")
plt.ylabel('DOS')
plt.suptitle("DOS d'une couche de graphène")
plt.grid()
plt.show() # trace le graphique

