import matplotlib.pyplot as plt
import numpy as np

#constants
radius = 0.02   #in metres
water_density = 1000  #in Kg
volume = 4/3*np.pi * radius**3
surface = np.pi * radius**2


files = [335, 336, 337, 338, 339, 341, 340, 342, 343, 344]
masses = [23.61, 25.96, 27.95, 28.90, 30.50, 31.00, 31.13, 32.2, 32.8, 33.51]    #in grams
stresses = []
frame_rate = [50, 50, 50, 50, 50, 50, 50, 50, 50, 5]
pix_to_cm = [2/115, 2/115, 2/115, 2/115, 2/115, 2/90, 2/115, 2/90, 2/90, 2/90]

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

for i in range(len(files)):
    file_name = f'/home/tommaso/Projects/creep_rising_sinking/data_tommaso/expo_velocity/center_position_{files[i]}.dat'
    
    space = []  #in cm
    time = []   #in seconds
    
    data = open(file_name, 'r')
    lines = data.readlines()
    
    offset = float(lines[0].split()[0])
    
    for count, line in enumerate(lines):
        obs = line.split()
        
        time.append(count/frame_rate[i])
        space.append((offset - float(obs[0]))*pix_to_cm[i])
    
    fit_x = [ t - time[limits[i][0]] for t in time[limits[i][0]:limits[i][1]]]
    fit_y = [ s - space[limits[i][0]] for s in space[limits[i][0]:limits[i][1]]]
    
    m, b = np.polyfit(fit_x, fit_y, 1)
    speeds.append(np.abs(m))
    
    stresses.append(9.81*(volume*water_density - masses[i]/1000)/surface)
    
    ball_density = masses[i]/1000/volume
    
    density_diff = water_density - ball_density
    
    etas.append((2*density_diff*(radius**2)*9.81)/(9*m/100))
    
    
    
    #fig = plt.figure()
    #ax = plt.gca()
    #
    #plt.ylabel(r'$\delta [cm]$',fontsize=15)
    #plt.xlabel(r'$t [s]$',fontsize=15)
    #plt.plot(time, space)
    #plt.axvline(x=limits[i][0]/frame_rate[i], color='r', linestyle='--', alpha=0.5)
    #plt.axvline(x=limits[i][1]/frame_rate[i], color='r', linestyle='--', alpha=0.5)
    #
    ##linear fit
    #fitted_points = [m*t - m*time[limits[i][0]] for t in time[limits[i][0]:limits[i][2]]]
    #plt.plot(time[limits[i][0]:limits[i][2]], fitted_points, 'b--')
    #plt.title(f"intruder mass {masses[i]}")
    #plt.show()
    
fig2 = plt.figure()
ax2 = plt.gca()

y2 = [2e4*np.exp(-0.3*mass) for mass in masses]

plt.ylabel(r'$\dot{\delta}$',fontsize=15)
plt.xlabel(r'$mass [g]$',fontsize=15)
plt.yscale('log')
plt.plot(masses, speeds, 'k+')
plt.plot(masses, y2, 'r--', alpha = 0.8, label=r'$y=2e4 e^{- 0.3x}$')
plt.legend()
plt.show()

fig3 = plt.figure()
ax3 = plt.gca()

y2 = [20*np.exp(-0.03*stress) for stress in stresses]

plt.ylabel(r'$\eta_{eff} [\text{Pa} \cdot s]$',fontsize=15)
plt.xlabel(r'$\sigma_S$ [N/$m^2$]',fontsize=15)
plt.yscale('log')
plt.plot(stresses[:-1], etas[:-1], 'k+')
plt.plot(stresses[:-1], y2[:-1], 'r--', alpha = 0.8, label=r'$y=2e4 e^{- 0.3x}$')
plt.legend()
plt.show()