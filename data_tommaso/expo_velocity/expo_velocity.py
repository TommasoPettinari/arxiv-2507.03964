import matplotlib.pyplot as plt
import numpy as np

files = [335, 336, 337, 338, 339, 341, 340, 342, 343, 344]
masses = [23.61, 25.96, 27.95, 28.90, 30.50, 31.00, 31.13, 32.2, 32.8, 33.6]
frame_rate = [50, 50, 50, 50, 50, 50, 50, 50, 50, 5]
pix_to_cm = 2/115

limits = [[90, 130, 150],
          [110, 170, 200],
          [110, 200, 280],
          [65, 330, 350],
          [100, 500, 800],
          [80, 350, 350],
          [120, 600, 1000],
          [100, 500, 650],
          [100, 1250, 1250],
          [5, 750, 750]]
speeds = []
etas = []

for i in range(10):
    file_name = f'/home/tommaso/Projects/creep_rising_sinking/data_tommaso/expo_velocity/center_position_{files[i]}.dat'
    
    space = []  #in pixels
    time = []   #in seconds
    
    data = open(file_name, 'r')
    lines = data.readlines()
    
    offset = float(lines[0].split()[0])
    
    for count, line in enumerate(lines):
        obs = line.split()
        
        time.append(count/frame_rate[i])
        space.append((offset - float(obs[0]))*pix_to_cm)
    
    fit_x = [ t - time[limits[i][0]] for t in time[limits[i][0]:limits[i][1]]]
    fit_y = [ s - space[limits[i][0]] for s in space[limits[i][0]:limits[i][1]]]
    
    m, b = np.polyfit(fit_x, fit_y, 1)
    speeds.append(np.abs(m))
    
    
    
    fig = plt.figure()
    ax = plt.gca()
    
    plt.ylabel(r'$\delta [cm]$',fontsize=15)
    plt.xlabel(r'$t [s]$',fontsize=15)
    plt.plot(time, space)
    plt.axvline(x=limits[i][0]/frame_rate[i], color='r', linestyle='--', alpha=0.5)
    plt.axvline(x=limits[i][1]/frame_rate[i], color='r', linestyle='--', alpha=0.5)
    
    #linear fit
    fitted_points = [m*t - m*time[limits[i][0]] for t in time[limits[i][0]:limits[i][2]]]
    plt.plot(time[limits[i][0]:limits[i][2]], fitted_points, 'b--')
    plt.title(f"intruder mass {masses[i]}")
    plt.show()
    
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