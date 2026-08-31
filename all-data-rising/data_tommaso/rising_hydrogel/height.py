import matplotlib.pyplot as plt
import numpy as np

file_name = '/home/tommaso/Projects/rheodata/rising hydrogel level/data.dat'

masses = []
height_hydro = []
height_water = []
normalized_height = []
volume = []
density = []

data = open(file_name, 'r')
lines = data.readlines()

lines.pop(0)

for count, line in enumerate(lines):
    obs = line.split()
    
    masses.append((float(obs[0])))
    height_hydro.append((float(obs[1])))
    height_water.append((float(obs[2])))
    volume.append((float(obs[3])))
    
    normalized_height.append(height_hydro[count]/height_water[count])
    density.append(masses[count]/volume[count])
    

fig = plt.figure()
ax = plt.gca()
plt.ylabel(r'$h_h / h_w$',fontsize=15)
plt.xlabel(r'$\rho$',fontsize=15)

x1 = np.linspace(0.8,3,100)
y1 = 0.33*x1

ax.plot(x1, y1, 'g--')
ax.plot(density, normalized_height, 'k+')
ax.axhline(y=1, color='r', linestyle='--')
plt.axvline(x = 2.95, color = 'b')
plt.show()