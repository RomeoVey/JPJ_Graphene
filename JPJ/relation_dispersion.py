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

N = 100

t = 0.1

a = 1

kx = linspace(-3,3,N)
ky = linspace(-3,3,N)

def energy(kx,ky,t,a):
    e_plus = t*sqrt(1 + 4*(cos(ky*a*sqrt(3)/2)**2 + cos(kx*3*a/2)*cos(ky*a*sqrt(3)/2)))
    e_moins = - e_plus
    return e_plus,e_moins

