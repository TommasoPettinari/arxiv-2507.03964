import matplotlib.pyplot as plt
import numpy as np


####SIMULATION PARAMETERS#######

#relative velocity at infinity
Ux = 0.5
Uy = 0.

##if the ball is accelerating is useful to consider the lab frame
##there we define the ball velocity at inf and the fluid velocity at inf
#uBall_inf_x = 0.
#uBall_inf_y = 0.
#
#uFluid_inf_x = Ux - uBall_inf_x
#uFluid_inf_y = Uy - uBall_inf_y

a = 3.      #ball radius


########FUNCTIONS########


def numerically_find_flow_line_2d(starting_point, step, vf_x, vf_y, x_lim, y_lim):
    func_x = []
    func_y = []
    
    # Find the points to the "right" of the point
    x = starting_point[0]
    y = starting_point[1]
    while True:
        func_x.append(x)
        func_y.append(y)
        
        delta_x = vf_x(x, y)
        delta_y = vf_y(x, y)
        
        # Normalize to unit vector and scale by 1/step
        magnitude = np.sqrt(delta_x ** 2 + delta_y ** 2)
        delta_x = delta_x / magnitude * step
        delta_y = delta_y / magnitude * step
        
        # Get the new point
        x = x + delta_x
        y = y + delta_y
        
        # Break if either x or y is outside of bounds
        if x > x_lim[1] or x < x_lim[0] or y > y_lim[1] or y < y_lim[0]:
            break
        
        # Break if it's taking too long for the curve to leave the bounds
        if len(func_x) > 100000:
            break
    
    # Find the points to the "left" of the point- nearly identical 
    x = starting_point[0]
    y = starting_point[1]
    while True:
        func_x.insert(0, x)
        func_y.insert(0, y)
        
        delta_x = vf_x(x, y)
        delta_y = vf_y(x, y)
        
        magnitude = np.sqrt(delta_x ** 2 + delta_y ** 2)
        delta_x = delta_x / magnitude * step
        delta_y = delta_y / magnitude * step
        
        # Here, we subtract delta instead of add
        x = x - delta_x
        y = y - delta_y
        
        if x > x_lim[1] or x < x_lim[0] or y > y_lim[1] or y < y_lim[0]:
            break
        if len(func_x) > 200000:
            break 
    
    return func_x, func_y


def stokes_solution_x(x,y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan(y/x)
    
    v_r = Ux * np.cos(theta) * (1 - 3*a/(2*r) + 0.5*(a/r)**3)
    v_theta = - Ux * np.sin(theta) * (1 - 3*a/(4*r) - 0.25*(a/r)**3)
    
    v_x = v_r * np.cos(theta) - v_theta * np.sin(theta)
    
    return v_x


def stokes_solution_y(x,y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan(y/x)
    
    v_r = Ux * np.cos(theta) * (1 - 3*a/(2*r) + 0.5*(a/r)**3)
    v_theta = - Ux * np.sin(theta) * (1 - 3*a/(4*r) - 0.25*(a/r)**3)
    
    v_y = v_r * np.sin(theta) + v_theta * np.cos(theta)
    
    return v_y


########VISUALISATION########

x_lim = (-10, 10)
y_lim = (-10, 10)

grid_step = 0.4
scale = 1

X, Y = np.meshgrid(np.arange(x_lim[0], x_lim[1], grid_step), np.arange(y_lim[0], y_lim[1], grid_step))
U = np.zeros(X.shape)
V = np.zeros(Y.shape)

for i in range(X.shape[0]):
    for j in range(Y.shape[0]):
        if (X[i,j]**2 + Y[i,j]**2 > a*2):
            U[i,j] = stokes_solution_x(X[i, j], Y[i, j])
            V[i,j] = stokes_solution_y(X[i, j], Y[i, j])
        else:
            U[i,j] = 0.
            V[i,j] = 0.
          

##QUIVER PLOT

fig, ax = plt.subplots()
Q = ax.quiver(X, Y, U, V, units='width', color='blue', alpha=0.7) # ,units='xy', scale=scale
qk = ax.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{[L]}{[T]}$', labelpos='E', coordinates='figure')

ball = plt.Circle((0, 0), a, color='r')          
ax.add_patch(ball)

##FLOW LINES
int_step = 0.02
func_x = []
func_y = []


func_x, func_y = numerically_find_flow_line_2d([-5,0.05], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
plt.plot(func_x, func_y, 'm')
func_x, func_y = numerically_find_flow_line_2d([-5,-0.05], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
plt.plot(func_x, func_y, 'm')


for i in range(1,10):
    func_x, func_y = numerically_find_flow_line_2d([-5,float(i)/5.], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
    plt.plot(func_x, func_y, 'g')
    func_x, func_y = numerically_find_flow_line_2d([-5,-float(i)/5.], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
    plt.plot(func_x, func_y, 'g')
    
for i in range(2,9):
    func_x, func_y = numerically_find_flow_line_2d([-5,i], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
    plt.plot(func_x, func_y, 'k')
    func_x, func_y = numerically_find_flow_line_2d([-5,-i], int_step, stokes_solution_x, stokes_solution_y, x_lim, y_lim)
    plt.plot(func_x, func_y, 'k')

plt.xlim(x_lim)
plt.ylim(y_lim)
ax.set_title('Velocity flow field')

##SCALAR HEATMAP PLOT

fig_vx, ax_vx = plt.subplots()
fig_vy, ax_vy = plt.subplots()

half = int(Y.shape[0] / 2)

horizontal = ax_vx.imshow(np.transpose(U[:half,:]), origin='upper', extent=(10,0,-10,10))
fig_vx.colorbar(horizontal)
ax_vx.set_title('y velocity')
ball2 = plt.Circle((0, 0), (a), color='k') 
ax_vx.add_patch(ball2)

vertical = ax_vy.imshow(np.transpose(V[:half,:]), origin='upper', extent=(10,0,-10,10))
fig_vy.colorbar(vertical)
ax_vy.set_title('x velocity')
ball3 = plt.Circle((0, 0), (a), color='k') 
ax_vy.add_patch(ball3)

plt.show()
