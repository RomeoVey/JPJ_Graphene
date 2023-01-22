from cmath import*
from numpy import*
import matplotlib.pyplot as plt
from scipy import *
from time import*

def plot_energy(L):
    t = 1.0
    s = 0.5
    e0 = 0.0
    eps = 0.00001
    psimoins = zeros((2*L+1, 2*L+1, 2))
    psiplus = zeros((2*L+1, 2*L+1, 2))
    listea = zeros(L+1)
    listeb = zeros(L+1)

    psimoins[L,L,0] = 1.0
    psiplus[L,L, 0] = 0.0
    a = e0
    b = 0.0 #b1
    listea[0] = a #rentre le a0 dans la liste
    listeb[0] = b #rentre le b0 dans la liste
    times = []
    ti = time()

    for n in range(1,L):#parcours les nombres de 1 à L-1 compris
        aplus=0  #remet à zero la variable locale
        for i in range(L-n,L+n+1):    
            for j in range(L-n,L+n+1):#application du Hamiltonien
                psiplus[i,j,0] = -b*psiplus[i,j,0]
                psiplus[i,j,1] = -b*psiplus[i,j,1]
        for i in range(L-n,L+n+1):    
            for j in range(L-n,L+n+1) :#application du Hamiltonien à psi_n
                if n%2==0: #((abs(i)%2 == 0) and (abs(j)%2 == 0)) or ((abs(i)%2 == 1) and (abs(j)%2 == 1))
                    psiplus[i,j,0] = psiplus[i,j,0]+t*(psimoins[i-1,j,0]+ psimoins[i,j-1,0]+psimoins[i,j+1,0])+s*psimoins[i,j,1]+e0*psimoins[i,j,0]
                    psiplus[i,j,1] = psiplus[i,j,1]+t*(psimoins[i-1,j,1]+ psimoins[i,j-1,1]+psimoins[i,j+1,1])+s*psimoins[i,j,0]+e0*psimoins[i,j,1]
            else: 
                psiplus[i,j,0] = psiplus[i,j,0]+t*(psimoins[i+1,j,0]+ psimoins[i,j-1,0]+psimoins[i,j+1,0])+s*psimoins[i,j,1]+e0*psimoins[i,j,0]
                psiplus[i,j,1] = psiplus[i,j,0]+t*(psimoins[i+1,j,1]+ psimoins[i,j-1,1]+psimoins[i,j+1,1])+s*psimoins[i,j,0]+e0*psimoins[i,j,1]
        for i in range(L-n,L+n+1):    
            for j in range(L-n,L+n+1) :  #parcours les nombres de L-n à L+n compris
                aplus = aplus + psimoins[i,j,0]*psiplus[i,j,0] + psimoins[i,j,1]*psiplus[i,j,1]
        a = aplus
        bplus = 0  #remet à zero la variable locale
        for i in range(L-n,L+n+1):    
            for j in range(L-n,L+n+1) :#parcours les nombres de L-n à L+n compris
                    psiplus[i,j,0] = psiplus[i,j,0]-a*psimoins[i,j,0]
                    psiplus[i,j,1] = psiplus[i,j,1]-a*psimoins[i,j,1]
                    bplus = bplus + (psiplus[i,j,0])**2 + (psiplus[i,j,1])**2 
        bplus = bplus**(1/2)
        for i in range(L-n,L+n+1):
            for j in range(L-n,L+n+1) :  #parcours les nombres de L-n à L+n compris
                psiplus[i,j,0] = (1/bplus)*psiplus[i,j,0]
                psiplus[i,j,1] = (1/bplus)*psiplus[i,j,1]
        temp = psiplus  #pour inverser psiplus et psimoins
        psiplus = psimoins  #pour inverser psiplus et psimoins
        psimoins = temp  #pour inverser psiplus et psimoins
        b = bplus
        
        listea[n] = a  #modifie la liste petit à petit
        listeb[n] = b
        print("Temps iteration :", round(time()-ti,2),"s\t"+str(n)+"/"+str(L))
        times.append(time()-ti)
        #with open('Cours_modelemicro/time3D.txt', 'a') as w:
        #    w.write('\n'+str(n)+'\t'+str(time()-ti))
        ti = time()

    #Calculs des coefs infinis
    ainf = 0.0
    binf = 0.0        
    for k in range(L-10,L):  #parcours les nombres de L-10 à L-1 compris
        ainf = ainf + listea[k]
        binf = binf + listeb[k]
    ainf = ainf/10
    binf = binf/10
    listea[L] = ainf
    listeb[L] = binf
    print("a=",listea)
    print("b=",listeb)

    #Calcul de la densité sous forme de fonction        
    def densite(E):
        G = (E - ainf -1j*sqrt(-(E - ainf)**2 + 4*binf**2))/(2*binf**2)  #calcul de Ginfini
        for i in range(L,0,-1): #parcours les nombres de L à 0 compris
            print("in densité i=",i,"b=",listeb[i])
            G = 1/(E - listea[i] - G*listeb[i]**2)  #calcul de la fraction continue
        nmu = abs((1/pi)*G.imag)  #calcul de la densite en prenant la valeur absolue
        return nmu

    energy = arange(ainf - 2*binf,ainf + 2*binf, 0.01)  #découpe des valeurs d'énergie
    dos = densite(energy)       
    plt.plot(energy,dos) 
    plt.grid()
    plt.show() # trace le graphique

'''

def read_time():
    with open('Cours_modelemicro/time3D.txt','r') as w:
        lines = w.read().splitlines()
        list = []
        for line in lines:
            list.append(line.split('\t'))
    iter = []
    time_iter = []
    for i in range(len(list)):
        iter.append(float(list[i][0]))
        time_iter.append(float(list[i][1]))
    fit_compl = polyfit(iter,time_iter,3)
    x = linspace(min(iter),max(iter),10*len(iter))
    fit = []
    for i in range(len(x)):
        fit.append(fit_compl[0]*x[i]**3+fit_compl[1]*x[i]**2+fit_compl[2]*x[i]+fit_compl[3])
    plt.plot(iter,time_iter)
    plt.plot(x,fit, color = 'green')
    plt.show()
'''
plot_energy(50)
#read_time()