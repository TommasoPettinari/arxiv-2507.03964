import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import os.path

fig = plt.figure()
ax = plt.gca()


time = []
displacement = []

#file_name = f'/home/tommaso/Projects/rheodata/creep_rising/2024-06-25.mqa'
file_name = '/home/tommaso/Projects/rheodata/sinking_with_liquid/sinking_with_liquid_3bis.dat'

data = open(file_name, 'r') #, encoding="utf-16"
lines = data.readlines()

#lines.pop(0)

for count in range(66):
    obs = lines[count].split()
    #time.append((float(obs[1])))
    displacement.append((float(obs[0])))
    
y_offset = displacement[0]

pix_to_cm = 4/733
displacement = [(i - y_offset)*pix_to_cm for i in displacement]
interval = 10
time = [n*interval for n in range(len(displacement))]

ax.scatter(time, displacement, marker='.', linewidths=1, label=f'data')

#ax.set_xscale('log')
#ax.set_yscale('log')

#p = np.polyfit(time[9000:70000], displacement[9000:70000], 1)

xlin1 = np.linspace(1, 800, 10^5)
ylin1 = 0.015*xlin1 + 0.5
ylin2 = 0.35*np.sqrt(xlin1)
plt.plot(xlin1, ylin2, c='red', label=f'y=sqrt(x)')
plt.plot(xlin1, ylin1, c='black', linestyle='dashed', label=f'y= 0.015 x')


#x1 = np.linspace(0,10**3,10**6)
#y1 = 20*x1**(0.5)
#plt.plot(x1,y1, c='red', label=f'y=t^{1/2}')
ax.set_xlabel('time (s)', fontsize = 15.0)
ax.set_ylabel('displacement (cm)', fontsize = 15.0)
ax.legend(loc='upper left', fontsize=12.0)            
plt.show()
fig.savefig(f'sinking_with_liquid/sinking_1bis.png', format='png', bbox_inches='tight')          