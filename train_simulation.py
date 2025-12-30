import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as animation

l = 0.4
gamma = 0.9 # spring damping
r0 = 1
w = 10
b = 0.7
d = 0.5
atan = np.arctan(b/d)


t_span = (0, 50)
t_eval = np.linspace(*t_span, 5000)
frame_frac = 5

k = 2000000
m = 15000
l0 = 10

n = 8 #no bogies
# s0 = [l0*i for i in range(n)]
# # s0[3] += 1
# # s0 = [0,l0,2*l0+0.2,3*l0]
# y0 = [0] *n
# # y0 = [0,0.2,0,0]
# # phi0 = [0.00,0.1,0,0,0,0,0,0]
# phi0 = [0]*n
# u0 = s0
# v0 = [0]*n
# vsx0 = [0]*n
# vsy0 = [0]*n
# state0 = s0+y0+phi0+u0+v0+vsx0+vsy0


# Finding the track
#Geen bocht
# p = [0,10000,200000,600000]
# sharp = 0.04
# h = 0.1
# off = 4.6
# s0 = [l0*i for i in range(n)]
# y0 = [0] *n
# phi0 = [0]*n
# u0 = s0
# v0 = [0]*n
# vsx0 = [0]*n
# vsy0 = [0]*n
# state0 = s0+y0+phi0+u0+v0+vsx0+vsy0

# bocht 1
p = [0,10,200,600]
sharp = 0.04
h = 0.1
off = 4.6
s0 = [l0*i for i in range(n)]
y0 = [0] *n
phi0 = [0]*n
u0 = s0
v0 = [0]*n
vsx0 = [0]*n
vsy0 = [0]*n
state0 = s0+y0+phi0+u0+v0+vsx0+vsy0

#bocht 2
# p = [0,100,200,600]
# sharp = 0.04
# h = 2
# off = 4.6
# s0 = [l0*i for i in range(n)]
# y0 = [0] *n
# phi0 = [0]*n
# u0 = s0
# v0 = [0]*n
# vsx0 = [0]*n
# vsy0 = [0]*n
# state0 = s0+y0+phi0+u0+v0+vsx0+vsy0

# #veer
# p = [0,10,200,600]
# sharp = 0.04
# h = 0.1
# off = 4.6
# s0 = [l0*i for i in range(n)]
# s0[3]+=1/10 * l0
# y0 = [0] *n
# phi0 = [0]*n
# u0 = s0
# v0 = [0]*n
# vsx0 = [0]*n
# vsy0 = [0]*n
# state0 = s0+y0+phi0+u0+v0+vsx0+vsy0

#crazy
# p = [0,100,200,600]
# sharp = 0.04
# h = 5
# off = 4.6
# s0 = [l0*i for i in range(n)]
# y0 = [0] *n
# phi0 = [0]*n
# u0 = s0
# v0 = [0]*n
# vsx0 = [0]*n
# vsy0 = [0]*n
# state0 = s0+y0+phi0+u0+v0+vsx0+vsy0

def find_theta(s):
    if p[0]<s<p[1]:
        return 0
    elif p[1]<s<p[2]:
        return h*(np.arctan(sharp*(s-p[1])-off)-np.arctan(-off))
    elif p[2]<s<p[3]:
        return -h*(np.arctan(sharp*(s-p[2])-off)-np.arctan(sharp*(p[2]-p[1])-off))
    else:
        return -h*(np.arctan(sharp*(p[3]-p[2])-off)-np.arctan(sharp*(p[2]-p[1])-off))
    

def find_dtheta(s):
    if p[0]<s<p[1]:
        return 0
    elif p[1]<s<p[2]:
        return h*sharp/(1+(sharp*(s-p[1])-off)**2)
    elif p[2]<s<p[3]:
        return -h*sharp/(1+(sharp*(s-p[2])-off)**2)
    else: return 0


def system(t, state):
    s = np.array(state[0:n])
    z = np.array(state[n:2*n])
    phi = np.array(state[2*n:3*n])
    x = np.array(state[3*n:4*n])
    y = np.array(state[4*n:5*n])
    vx_slip = np.array(state[5*n:6*n]) 
    vy_slip =np.array(state[6*n:7*n]) 
    theta = np.array([find_theta(i) for i in s])
    ds_dt, dz_dt, dphi_dt, dx_dt, dy_dt, dvx_slip_dt, dvy_slip_dt = np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n),np.zeros(n)

    ds_dt = w * r0 * np.cos(phi)
    dz_dt = w * r0 * np.sin(phi)


    gx = x - z*np.sin(theta)
    gy = y + z*np.cos(theta)


    Fx,Fy = np.zeros(n-1),np.zeros(n-1)

    for i in range(n-1):
        dx = gx[i+1] - gx[i]
        dy = gy[i+1] - gy[i]
        dis = np.sqrt(dx**2 + dy**2)
        ux, uy = dx/dis, dy/dis  
        Fx[i] = k * (dis - l0) * ux 
        Fy[i] = k * (dis - l0) * uy



    dvx_slip_dt[1:] -= Fx/m 
    dvx_slip_dt[:n-1] += Fx/m
    dvy_slip_dt[1:] -= Fy/m 
    dvy_slip_dt[:n-1] += Fy/m

    dvx_slip_dt -= gamma*dvx_slip_dt
    dvy_slip_dt -= gamma*dvy_slip_dt

    for i in range(n):
        dphi_dt[i] = -z[i]*w*l*np.sin(atan)/np.sqrt(b**2+d**2) - find_dtheta(s[i]) * ds_dt[i]
    
    ds_dt += vx_slip * np.cos(theta) + vy_slip * np.sin(theta)
    dz_dt += vx_slip * -np.sin(theta) + vy_slip * np.cos(theta)
    dx_dt = np.cos(theta)*ds_dt
    dy_dt = np.sin(theta)*ds_dt 
    return np.concatenate((ds_dt, dz_dt, dphi_dt,dx_dt,dy_dt,dvx_slip_dt,dvy_slip_dt))

sol = solve_ivp(system, t_span, state0, t_eval=t_eval, method='RK45')  

s_values = [sol.y[i] for i in range(n)]

x_center = np.array([sol.y[3*n+i] for i in range(n)])
y_center = np.array([sol.y[4*n+i] for i in range(n)])
theta_values = [[find_theta(s) for s in s_values[j]] for j in range(n)]

# Compute left and right boundaries
x_left, y_left, x_right, y_right = n*[[]],n*[[]],n*[[]],n*[[]]



for i in range(n):
    x_left[i] = [x - b * np.sin(theta) for x, theta in zip(x_center[i], theta_values[i])]
    y_left[i] = [y + b * np.cos(theta) for y, theta in zip(y_center[i], theta_values[i])]
    x_right[i] = [x + b * np.sin(theta) for x, theta in zip(x_center[i], theta_values[i])]
    y_right[i] = [y - b * np.cos(theta) for y, theta in zip(y_center[i], theta_values[i])]


def update(frame):
    frame *= frame_frac
    r = np.sqrt(b**2 + d**2)
    
    for i in range(n):
        phi = sol.y[2*n+ i][frame]
        z = sol.y[n + i][frame]
        s = sol.y[i][frame]
        cost = np.cos(find_theta(s))
        sint = np.sin(find_theta(s))


        flwx = r * np.cos(phi + atan)
        flwy = z + r * np.sin(phi + atan)
        frwx = r * np.cos(phi - atan)
        frwy = z + r * np.sin(phi - atan)
        blwx = r * np.cos(np.pi + phi - atan)
        blwy = z + r * np.sin(np.pi + phi - atan)
        brwx = r * np.cos(np.pi + phi + atan)
        brwy = z + r * np.sin(np.pi + phi + atan)

        animated_flw[i].set_data([flwx * cost - flwy * sint + x_center[i][frame]], [flwx * sint + flwy * cost + y_center[i][frame]])
        animated_frw[i].set_data([frwx * cost - frwy * sint + x_center[i][frame]], [frwx * sint + frwy * cost + y_center[i][frame]])
        animated_blw[i].set_data([blwx * cost - blwy * sint + x_center[i][frame]], [blwx * sint + blwy * cost + y_center[i][frame]])
        animated_brw[i].set_data([brwx * cost - brwy * sint + x_center[i][frame]], [brwx * sint + brwy * cost + y_center[i][frame]])
        animated_center[i].set_data([x_center[i][frame]], [y_center[i][frame]])

    wid = (l0-1.5*d)/2 
    hig = b
    bl,br,fl,fr = (-wid,hig),(-wid,-hig),(wid,hig),(wid,-hig)

    x_mid = [x_center[i][frame] + sol.y[n+i][frame] * -np.sin(theta_values[i][frame]) for i in range (n)]
    y_mid = [y_center[i][frame] + sol.y[n+i][frame] * np.cos(theta_values[i][frame]) for i in range (n)]

    for i in range(n-1):
        dx = x_mid[i+1] - x_mid[i]
        dy = y_mid[i+1] - y_mid[i]
        mid_x = (x_mid[i] + x_mid[i+1]) / 2 
        mid_y = (y_mid[i] + y_mid[i+1]) / 2
        a = np.arctan2(dy,dx)
        bl1 = [bl[0]*np.cos(a)-bl[1]*np.sin(a)+mid_x,bl[0]*np.sin(a)+bl[1]*np.cos(a)+mid_y]
        br1 = [br[0]*np.cos(a)-br[1]*np.sin(a)+mid_x,br[0]*np.sin(a)+br[1]*np.cos(a)+mid_y]
        fl1 = [fl[0]*np.cos(a)-fl[1]*np.sin(a)+mid_x,fl[0]*np.sin(a)+fl[1]*np.cos(a)+mid_y]
        fr1 = [fr[0]*np.cos(a)-fr[1]*np.sin(a)+mid_x,fr[0]*np.sin(a)+fr[1]*np.cos(a)+mid_y]
        animated_tl[i].set_data([bl1[0],fl1[0]],[bl1[1],fl1[1]])
        animated_bl[i].set_data([br1[0],fr1[0]],[br1[1],fr1[1]])
        animated_rl[i].set_data([fl1[0],fr1[0]],[fl1[1],fr1[1]])
        animated_ll[i].set_data([bl1[0],br1[0]],[bl1[1],br1[1]])

    ax1.set_xlim([np.mean(x_center[:,frame]) - (1 + n*l0)/1.5, np.mean(x_center[:,frame]) + (1 + n*l0)/1.5])
    ax1.set_ylim([np.mean(y_center[:,frame]) - (1 + n*l0)/1.5, np.mean(y_center[:,frame]) + (1 + n*l0)/1.5])



    return [*animated_flw, *animated_frw, *animated_blw, *animated_brw, *animated_center, *animated_tl, *animated_bl, *animated_rl, *animated_ll]


fig,(ax1,ax2) = plt.subplots(1,2,figsize = (15,15))


ax1.plot(x_left[0], y_left[0], color='black')
ax1.plot(x_right[0], y_right[0],color='black')
ax2.plot(x_left[0], y_left[0], color='black')
ax2.plot(x_right[0], y_right[0],color='black')



animated_flw, animated_frw, animated_blw, animated_brw, animated_center = [], [], [], [], []

animated_tl, animated_bl, animated_rl, animated_ll = [],[],[],[]

for _ in range(n):
    animated_flw.append(ax1.plot([], [], 'o', markersize=3, color='blue')[0])
    animated_frw.append(ax1.plot([], [], 'o', markersize=3, color='blue')[0])
    animated_blw.append(ax1.plot([], [], 'o', markersize=3, color='blue')[0])
    animated_brw.append(ax1.plot([], [], 'o', markersize=3, color='blue')[0])
    animated_center.append(ax2.plot([], [], 's', markersize=5, color='blue')[0])
    
for _ in range(n-1):
    animated_tl.append(ax1.plot([], [], color = 'orange',lw=2)[0])
    animated_bl.append(ax1.plot([], [], color = 'orange',lw=2)[0])
    animated_rl.append(ax1.plot([], [], color = 'orange',lw=2)[0])
    animated_ll.append(ax1.plot([], [], color = 'orange',lw=2)[0])


ani = animation.FuncAnimation(fig=fig, func= update, frames = int(len(t_eval)/frame_frac), interval = (t_eval[1]-t_eval[0])*1000*frame_frac)
ax1.set_ylabel('y (m)')
ax1.set_xlabel('x (m)')
ax2.set_ylabel('y (m)')
ax2.set_xlabel('x (m)')

ax2.set_aspect('equal')
ax1.set_aspect('equal')
plt.show()


#_______________________________________________________________________________________________________________


