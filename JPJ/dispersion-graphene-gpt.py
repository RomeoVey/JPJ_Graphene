import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the lattice constant
a = 1

N = 100

# Define the hopping parameter
t = 2.8

# Define the energy dispersion function
def energy(kx, ky):
    return -t*np.sqrt(1 + 4*(np.cos(ky*a*np.sqrt(3)/2)**2 + np.cos(kx*3*a/2)*np.cos(ky*a*np.sqrt(3)/2)))

# Generate a grid of k-points in the Brillouin zone
kx, ky = np.mgrid[-np.pi/a:np.pi/a:1/N, -np.pi/a:np.pi/a:1/N]

# Calculate the energy at each k-point
E = energy(kx, ky)

# Plot the energy dispersion in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(kx, ky, E, cmap='rainbow')
ax.plot_surface(kx, ky, -E, cmap='rainbow_r',alpha = 0.75)
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('Energy')
plt.show()
