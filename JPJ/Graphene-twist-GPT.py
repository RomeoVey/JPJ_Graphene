import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters for bi-layer graphene lattice
a = 1.42 # lattice constant
d = 0.335 # interlayer distance
theta = np.pi/6 # rotation angle of top layer
N = 100 # number of k-points

# Create k-mesh for bi-layer graphene lattice
kx, ky = np.meshgrid(np.linspace(-np.pi/a, np.pi/a, N), np.linspace(-np.pi/a, np.pi/a, N))

# Calculate energy dispersion
E1 = np.sqrt(1 + 4*np.cos(np.sqrt(3)*kx*a/2)*np.cos(3*ky*a/2) + 4*np.cos(3*kx*a/2)*np.cos(np.sqrt(3)*ky*a/2))
E2 = np.sqrt(1 + 4*np.cos(np.sqrt(3)*kx*a/2)*np.cos(3*ky*a/2) + 4*np.cos(3*kx*a/2)*np.cos(np.sqrt(3)*ky*a/2) + 8*np.cos(3*kx*a/2)*np.cos(3*ky*a/2)*np.cos(d*np.sqrt(4/3) - 2*np.pi*theta))

# Create 3D figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the band structure
ax.plot_surface(kx, ky, E1, cmap='viridis', alpha=0.5)
ax.plot_surface(kx, ky, E2, cmap='viridis', alpha=0.5)
ax.set_xlabel('kx')
ax.set_ylabel('ky')
ax.set_zlabel('Energy (eV)')
plt.show()
