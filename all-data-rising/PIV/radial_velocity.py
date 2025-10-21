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
#choose relevant interstep
interstep = 1
positions = np.linspace(-225, 225, int(450/interstep + 1))
positions = [p*4/175 for p in positions]

fig, ax = plt.subplots()
plt.ylabel(r'$|v|$',fontsize=18)
plt.xlabel(r'$r$',fontsize=18)
ax.set_ylim([0.0, 0.52])
ax.set_xlim([0, 7.2])
plt.axvspan(0, 2.0, facecolor='orange', alpha=0.5)
plt.axhline(y=0.45, linewidth=2, color='#d62728', label='ball speed')
ax.tick_params(axis='both', which='major', labelsize=15)

plt.axhline(y=0.015, linewidth=2, color='#d62728', linestyle='--', label='backflow speed')

file_name = f'all-data-rising/data_tommaso/radial_velocity/radial_vel_interstep{interstep}.txt'
data = open(file_name, 'r', encoding="utf-8")
lines = data.readlines()


for line in lines:
    obs = line.split()
    v_r_PIV.append(float(obs[0]))

mean = int((450 / interstep) / 2 - 1)
length = int((450 / interstep) / 2)
    
mean_v_r_PIV = [(v_r_PIV[mean - i] + v_r_PIV[mean + i])/2 for i in range(length)]
error_v_r_PIV = [np.abs(v_r_PIV[mean- i] - mean_v_r_PIV[i]) for i in range(length)]
    
ax.plot(xlin, v_r_stokes, color='k', label='Stokes solution')

ax.plot(positions[mean:(mean+length)], mean_v_r_PIV, 'bo--', label='PIV data')
ax.fill_between(positions[mean:(mean+length)], np.array(mean_v_r_PIV) - np.array(error_v_r_PIV), np.array(mean_v_r_PIV) + np.array(error_v_r_PIV), color='c', alpha=0.7)

plt.legend(fontsize=15)

plt.show()

###### vertical axis #######

v_vert_PIV = []
v_vert_component = []

#choose relevant interstep
interstep = 5
offset = 5
positions = np.linspace(-225 + offset, 225 + offset, int(450/interstep + 1))
positions = [p*4/175 for p in positions]

fig, ax = plt.subplots()
plt.ylabel(r'$velocity$',fontsize=28)
plt.xlabel(r'$r$',fontsize=28)
ax.set_ylim([0.0, 0.52])
#ax.set_xlim([-7.2, 7.2])
plt.axvspan(0, 2.0, facecolor='orange', alpha=0.5)
plt.axhline(y=0.45, linewidth=2, color='#d62728', label='ball speed')
ax.tick_params(axis='both', which='major', labelsize=25)

#load radial velocity
file_name = f'all-data-rising/data_tommaso/radial_velocity/radial_vel_vertical_interstep{interstep}.txt'
data = open(file_name, 'r', encoding="utf-8")
lines = data.readlines()

for line in lines:
    obs = line.split()
    v_vert_PIV.append(float(obs[0]))

#load vertical component only
file_name2 = f'all-data-rising/data_tommaso/radial_velocity/V_component_vertical_interstep{interstep}.txt'
data2 = open(file_name2, 'r', encoding="utf-8")
lines2 = data2.readlines()

for line in lines2:
    obs = line.split()
    v_vert_component.append(float(obs[0]))

#plotting

mean = int((450 / interstep) / 2 - 1)
length = int((450 / interstep) / 2)
    
#mean_v_r_PIV = [(v_r_PIV[mean - i] + v_r_PIV[mean + i])/2 for i in range(length)]
#error_v_r_PIV = [np.abs(v_r_PIV[mean- i] - mean_v_r_PIV[i]) for i in range(length)]

ax.plot(positions[mean:(mean+length)], v_vert_PIV[mean:(mean+length)], 'bo--', label='|v|, ahead')
ax.plot(positions[mean:(mean+length)], list(reversed(v_vert_PIV[:mean+1])), 'gs--', label='|v|, behind')

ax.plot(positions[mean:(mean+length)], v_vert_component[mean:(mean+length)], 'mo--', label='V component, ahead')
ax.plot(positions[mean:(mean+length)], list(reversed(v_vert_component[:mean+1])), 'cs--', label='V component, behind')

plt.legend(fontsize=25, loc='upper right')

plt.show()