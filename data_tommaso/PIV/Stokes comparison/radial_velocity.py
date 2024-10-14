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


xlin = np.linspace(-6.8, 6.8, 10000)
# The stoke solution is evaluated in the ball system of reference (with movid liquid)
# we now want to consider ourselves in the lab system of reference so we subtract the fluid velocity obtained to the ball velocity
v_r_stokes = [Ux - stokes_solution_r(x) if x**2 > a**2 else Ux for x in xlin]

v_r_PIV = []
positions = np.linspace(-225, 225, 91)
positions = [p*4/175 for p in positions]

ax1 = plt.subplot(1,2,1)
plt.ylabel(r'$|v^{\text{Stokes}}|$',fontsize=9)
plt.xlabel(r'$r$',fontsize=9)
ax1.set_ylim([0.0, 0.52])
plt.axvspan(-2.0, 2.0, facecolor='#2ca02c', alpha=0.5)
plt.axhline(y=0.45, linewidth=2, color='#d62728')
ax2 = plt.subplot(1,2,2)
plt.ylabel(r'$|v^{\text{hyd}}|$',fontsize=9)
plt.xlabel(r'$r$',fontsize=9)
ax2.set_ylim([0.0, 0.52])
plt.axvspan(-2.1, 2, facecolor='#2ca02c', alpha=0.5)
plt.axhline(y=0.45, linewidth=2, color='#d62728')
plt.axhline(y=0.015, linewidth=2, color='#d62728', linestyle='--')

file_name = f'/home/tommaso/Projects/PIV/stokes_sol/radial_vel_2.txt'
data = open(file_name, 'r', encoding="utf-8")
lines = data.readlines()


for line in lines:
    obs = line.split()
    v_r_PIV.append(float(obs[0]))
    
ax1.plot(xlin, v_r_stokes)

ax2.plot(positions, v_r_PIV)


plt.show()