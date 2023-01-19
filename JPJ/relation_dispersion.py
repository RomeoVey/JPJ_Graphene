"""
Created on Thu Jan 12 16:21:49 2023

@author: LeoMagik
"""
from scipy import *
from cmath import *
import matplotlib 
import matplotlib.pyplot as plt
from numpy import *
from mpl_toolkits import mplot3d

N = 1000

t = 0.1

a = 1

kx = linspace(-3,3,N)
ky = linspace(-3,3,N)

def energy(kx,ky,t,a):
    e_plus = t*sqrt(1 + 4*(cos(ky*a*sqrt(3)/2)**2 + cos(kx*3*a/2)*cos(ky*a*sqrt(3)/2)))
    e_moins = - e_plus
    return e_plus,e_moins

def energy_line(kx,ky,t,a):
    M_plus,M_moins = energy(0,ky,t,a)
    plt.plot(ky[N//2:],M_plus[N//2:])
    plt.plot(ky[N//2:],M_moins[N//2:]) 
    plt.grid()
    plt.show()
    
    
def plot_energy(kx,ky,t,a):
    ax = plt.axes(projection = '3d')
    M_plus = zeros((N,N))
    M_moins = zeros((N,N))
    for i in range(len(kx)):
        for j in range(len(ky)):
            M_plus[i,j] = energy(kx[i],ky[j],t,a)[0]
            M_moins[i,j] = energy(kx[i],ky[j],t,a)[1]

    step = 50
    l_kx = len(kx)//step
    l_ky = len(ky)//step
    kx_less = zeros(l_ky)
    ky_less = zeros(l_kx)
    MP = zeros((l_kx,l_ky))
    MM = zeros((l_kx,l_ky))
    for i in range(l_kx):
        ky_less[i] = ky[i*step]
    for j in range(l_ky):
        kx_less[j] = kx[j*step]

    for i in range(l_kx):
        for j in range(l_ky):
            MP[i,j] = M_plus[i*step,j*step]
            MM[i,j] = M_moins[i*step,j*step]
    for i in range(l_kx):
        kx_fix = kx_less[i]*ones(l_ky)
        ax.plot3D(kx_fix,ky_less,MP[i,:], color = 'blue')
        ax.plot3D(kx_fix,ky_less,MM[i,:], color = 'red')
    for j in range(l_ky):
        ky_fix = ky_less[j]*ones(l_kx)
        ax.plot3D(kx_less,ky_fix,MP[:,j], color = 'blue')
        ax.plot3D(kx_less,ky_fix,MM[:,j], color = 'red')
    plt.show()

# plot_energy(kx,ky,t,a)

energy_line(kx,ky,t,a)