# -*- coding: utf-8 -*-
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

N = 100

t = 0.1

a = 1

kx = linspace(-3,3,N)
ky = linspace(-3,3,N)

def energy(kx,ky,t,a):
    e_plus = t*sqrt(1 + 4*(cos(ky*a*sqrt(3)/2)**2 + cos(kx*3*a/2)*cos(ky*a*sqrt(3)/2)))
    e_moins = - e_plus
    return e_plus,e_moins

def plot_energy(kx,ky,t,a):
    ax = plt.axes(projection = '3d')
    M_plus = zeros((N,N))
    M_moins = zeros((N,N))
    for i in range(len(kx)):
        for j in range(len(ky)):
            M_plus[i,j] = energy(kx[i],ky[j],t,a)[0]
            M_moins[i,j] = energy(kx[i],ky[j],t,a)[1]

        kx_fix = kx[i]*ones(len(ky))
        ax.plot3D(kx_fix,ky,M_plus[i,:], color = 'blue')
        ax.plot3D(kx_fix,ky,M_moins[i,:], color = 'red')
    for j in range(len(ky)):
        ky_fix = ky[j]*ones(len(kx))
        ax.plot3D(kx,ky_fix,M_plus[:,j], color = 'blue')
        ax.plot3D(kx,ky_fix,M_moins[:,j], color = 'red')
    plt.show()

plot_energy(kx,ky,t,a)