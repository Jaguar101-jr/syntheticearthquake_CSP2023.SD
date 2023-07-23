## A simple synthetic earthquake program simulation, Centre for Seismological Phenomena - Sudan July 2023.
# Hussein Muhammed; Abdelhafiz Gadelmula Taha; Ammar Awad Ali.
#Modeling a synthetic earthquake involves simulating the generation and propagation of seismic waves.
#One common approach is to use finite difference methods to solve the wave equations numerically.

##This code creates a synthetic earthquake with a single earthquake source located at the center of the model.
##The seismic wavefield is simulated using finite difference methods, and the resulting displacement is visualized using a color map.
#####################################################################################################################################

import numpy as np
import matplotlib.pyplot as plt

# Parameters
# controling factors for accurate propagation are: 'vp' and 'dt' must be fast and small enough to insure wave spreading
nx = 100  # Number of grid points in x-direction
nz = 100  # Number of grid points in z-direction
dx = 10.0  # Grid spacing in meters
dz = 10.0
vp = 3500.0  # P-wave velocity in m/s
dt = 0.0001  # Time step in seconds
nt = 1000  # Number of time steps

# Initialize the wavefield
u = np.zeros((nz, nx))

# Source parameters (location and shape)
source_x = int(nx / 2)  # Center of the source (x-coordinate)
source_z = int(nz / 4)  # Center of the source (z-coordinate)
source_amplitude = 1.0  # Amplitude of the source

# Apply the earthquake source
u[source_z, source_x] = source_amplitude

# Finite difference modeling
for it in range(1, nt):
    for iz in range(1, nz - 1):
        for ix in range(1, nx - 1):
            laplacian = (u[iz, ix + 1] - 2.0 * u[iz, ix] + u[iz, ix - 1]) + (u[iz + 1, ix] - 2.0 * u[iz, ix] + u[iz - 1, ix])
            u[iz, ix] = 2.0 * u[iz, ix] - u[iz, ix] + (vp * vp * dt * dt) * laplacian

# Plot the seismic wavefield
plt.imshow(u, cmap='seismic', extent=[0, nx * dx, 0, nz * dz], aspect='auto')
plt.colorbar(label='Displacement')
plt.xlabel('Distance (m)')
plt.ylabel('Depth (m)')
plt.title('Synthetic Earthquake Model')
plt.show()
