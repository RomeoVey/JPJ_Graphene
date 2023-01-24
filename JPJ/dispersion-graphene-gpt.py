import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

N = 1000

# Define the lattice constant
a = 1

# Define the hopping parameter
t = 2.8

# Define the energy dispersion function
def energy(kx, ky):
    return -t*np.sqrt(1 + 4*(np.cos(ky*a*np.sqrt(3)/2)**2 + np.cos(kx*3*a/2)*np.cos(ky*a*np.sqrt(3)/2)))

# Generate a grid of k-points in the Brillouin zone
kx, ky = np.mgrid[-np.pi/a:np.pi/a:2*np.pi/(a*N), -np.pi/a:np.pi/a:2*np.pi/(a*N)]

kx_n, ky_n = np.linspace(-np.pi/a,np.pi/a,N),np.linspace(-np.pi/a,np.pi/a,N)

# Calculate the energy at each k-point
E = energy(kx, ky)

E_min = np.min(E)
E_max = np.max(E)

def loren(En,E0,gamma): return gamma/(np.pi*((En-E0)**2 + gamma**2))

def DOS(En,E0,gamma):
    DOSi = 0
    for i in range(N):
        for j in range(N):
            fx,fy = kx_n[i], ky_n[j]
            DOSi += loren(En,E0(fx,fy),gamma)
    return DOSi

gamma = 0.1
Enj = np.linspace(-np.max([abs(E_min),abs(E_max)]),np.max([abs(E_min),abs(E_max)]),N)
dos1 = DOS(Enj,energy,gamma)
dos2 = dos1[::-1]


fig = plt.figure(1)
plt.plot(Enj,dos1+dos2)



# Plot the energy dispersion in 3D
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(kx, ky, E, cmap='rainbow')
ax.plot_surface(kx, ky, -E, cmap='rainbow_r',alpha = 0.75)
ax.set_xlabel('kx',fontsize = 'large')
ax.set_ylabel('ky',fontsize = 'large')
ax.set_zlabel('Energy',fontsize = 'large')


def rotate(angle):
        ax.view_init(azim=angle)   


print("Making animation")
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('JPJ/Dispersion_surfaces.gif', dpi=80)
