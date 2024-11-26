import matplotlib.pyplot as plt
import numpy as np


Ux = 0.45
Uy = 0.
a = 2.0


def stokes_solution_r(y):

    r = np.abs(y)
    
    #v_r = Ux * np.cos(theta) * (1 - 3*a/(2*r) + 0.5*(a/r)**3)
    v_theta = Ux * (1 - 3*a/(4*r) - 0.25*(a/r)**3)
    
    return v_theta


xlin = np.linspace(0, 6.8, 10000)
# The stoke solution is evaluated in the ball system of reference (with movid liquid)
# we now want to consider ourselves in the lab system of reference so we subtract the fluid velocity obtained to the ball velocity
v_r_stokes = [Ux - stokes_solution_r(x) if x**2 > a**2 else Ux for x in xlin]

v_r_PIV = []
positions = np.linspace(-225, 225, 91)
positions = [p*4/175 for p in positions]

fig, ax = plt.subplots()
plt.ylabel(r'$|v|$',fontsize=18)
plt.xlabel(r'$r$',fontsize=18)
ax.set_ylim([0.0, 0.52])
ax.set_xlim([0, 7.2])
plt.axvspan(0, 2.0, facecolor='orange', alpha=0.5)
plt.axhline(y=0.45, linewidth=2, color='#d62728', label='ball speed')
ax.tick_params(axis='both', which='major', labelsize=15)

#ax2.set_xlim([0,6.5])
plt.axhline(y=0.015, linewidth=2, color='#d62728', linestyle='--', label='backflow speed')

file_name = f'all-data-rising/data_tommaso/radial_velocity/radial_vel.txt'
data = open(file_name, 'r', encoding="utf-8")
lines = data.readlines()


for line in lines:
    obs = line.split()
    v_r_PIV.append(float(obs[0]))
    
mean_v_r_PIV = [(v_r_PIV[44 - i] + v_r_PIV[44 + i])/2 for i in range(45)]
error_v_r_PIV = [np.abs(v_r_PIV[44- i] - mean_v_r_PIV[i]) for i in range(45)]
    
ax.plot(xlin, v_r_stokes, color='k', label='Stokes solution')

ax.plot(positions[44:89], mean_v_r_PIV, 'bo--', label='PIV data')
ax.fill_between(positions[44:89], np.array(mean_v_r_PIV) - np.array(error_v_r_PIV), np.array(mean_v_r_PIV) + np.array(error_v_r_PIV), color='c', alpha=0.7)

plt.legend(fontsize=15)

plt.show()