import matplotlib.pyplot as plt
import numpy as np

files = [335, 336, 337, 338, 339, 340]
masses = [23.61, 25.96, 27.95, 28.90, 30.50, 31.13]
frame_rate = 50
pix_to_cm = 2/115

limits = [[90, 130, 150],
          [110, 170, 200],
          [110, 200, 280],
          [65, 330, 350],
          [100, 500, 800],
          [120, 600, 1000]]
speeds = []
etas = []

for i in range(6):
    file_name = f'/home/tommaso/Projects/rheodata/expo_velocity/center_position_{files[i]}.dat'
    
    space = []  #in pixels
    time = []   #in seconds
    
    data = open(file_name, 'r')
    lines = data.readlines()
    
    offset = float(lines[0].split()[0])
    
    for count, line in enumerate(lines):
        obs = line.split()
        
        time.append(count/frame_rate)
        space.append((offset - float(obs[0]))*pix_to_cm)
    
    m, b = np.polyfit(time[limits[i][0]:limits[i][1]], space[limits[i][0]:limits[i][1]], 1)
    speeds.append(np.abs(m))
    
    
    
    #fig = plt.figure()
    #ax = plt.gca()
    #
    #plt.ylabel(r'$\delta [px]$',fontsize=15)
    #plt.xlabel(r'$t [s]$',fontsize=15)
    #plt.plot(time, space)
    #plt.axvline(x=limits[i][0]/frame_rate, color='r', linestyle='--', alpha=0.5)
    #plt.axvline(x=limits[i][1]/frame_rate, color='r', linestyle='--', alpha=0.5)
    #
    ##linear fit
    #fitted_points = [m*t + b for t in time[limits[i][0]:limits[i][2]]]
    #plt.plot(time[limits[i][0]:limits[i][2]], fitted_points, 'b--')
    #plt.show()
    
fig2 = plt.figure()
ax2 = plt.gca()

y2 = [2.5e6*np.exp(-0.5*mass) for mass in masses[2:]]

plt.ylabel(r'$\dot{\delta}$',fontsize=15)
plt.xlabel(r'$mass [g]$',fontsize=15)
plt.yscale('log')
plt.plot(masses, speeds, 'k+-')
plt.plot(masses[2:], y2, 'r--', alpha = 0.8, label=r'$y=2.5e6 e^{-0.5 x}$')
plt.legend()
plt.show()