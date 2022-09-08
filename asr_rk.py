#!/usr/bin/python3

# Simulates adiabatic spin rotation by solving the Bloch equations
# using Runge-Kutta integration

# by Jeff
# Edited by Emma

#B shift from +x to +z

from math import *
from scipy.constants import *
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


gamma=physical_constants['neutron gyromag. ratio'][0]
gamma=gamma/1e6 # rad/s/uT
# gamma is positive in the physical constants library
# to make it negative, you have to add a minus sign

# Set the field and rate at which it should change, here.
B1=1 # uT
T1=1 # s, time to go around once
omega1=2*pi/T1 # Hz
# At 1 uT, gamma*B1 = 2*pi*(30 Hz)
# So pick an omega1 that's a bit slower than that so we're fairly adiabatic
# I picked omega1=2*pi*(1 Hz) here.

s0=[0,1,0]      # starting polarization vector
t0=0.               # starting time
t1=T1               # ending time (seconds)
dt=.001           # time step (seconds)

i=0
sparse=1 # sparsification factor

def b(t):
    bx=0
    by=B1*cos(omega1*t)
    bz=B1*sin(omega1*t)
    return [bx,by,bz]  # uT

def dsdt(t,s):
    # The Bloch equations in the non-rotating frame
    # f = dS/dt = gamma S x B + relaxation
    sx,sy,sz=s
    bx,by,bz=b(t)
    dsxdt=-gamma*(sy*bz-sz*by)
    dsydt=-gamma*(sz*bx-sx*bz)
    dszdt=-gamma*(sx*by-sy*bx)
    return [dsxdt,dsydt,dszdt]

r=ode(dsdt).set_integrator('dop853',atol=1e-8)
r.set_initial_value(s0,t0)

t_save=[]
sx_save=[]
sy_save=[]
sz_save=[]
sx_perfect=[]
sy_perfect=[]
sz_perfect=[]
b_save=[]
while r.successful() and r.t < t1:
    r.integrate(r.t+dt) # do an integration step
    t=r.t               # time after the step
    s=r.y               # polarization vector after the step
    sx,sy,sz=s          # polarization vector after the step
    i+=1
    if(i%sparse==0):
        t_save.append(t)
        sx_save.append(sx)
        sy_save.append(sy)
        sz_save.append(sz)
        b_t=b(t)
        b_save.append(b_t)
        sx_perfect.append(b_t[0]/B1)
        sy_perfect.append(b_t[1]/B1)
        sz_perfect.append(b_t[2]/B1)

#had to add this to not throw error for newer version of matplotlib
t_save = np.array(t_save)
sx_save = np.array(sx_save)
sy_save = np.array(sy_save)
sz_save = np.array(sz_save)
b_save = np.array(b_save)
sx_perfect = np.array(sx_perfect)
sy_perfect = np.array(sy_perfect)
sz_perfect = np.array(sz_perfect)

# print(b_save[0][2])

#plt.plot(t_save,sx_save)
#plt.plot(t_save,sy_save)
#plt.plot(t_save,sz_save)
#plt.xlabel('Time (s)')
#plt.ylabel('Spin components')
#plt.show()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

line,=ax.plot(sx_save[:1],sy_save[:1],sz_save[:1], color='tomato')
line2,=ax.plot(sx_perfect[:1],sy_perfect[:1],sz_perfect[:1], color='royalblue')
a=Arrow3D([0,sx_save[0]],[0,sy_save[0]],[0,sz_save[0]],mutation_scale=20,arrowstyle="-|>",color="r")
a2=Arrow3D([0,sx_perfect[0]],[0,sy_perfect[0]],[0,sz_perfect[0]],mutation_scale=20,arrowstyle="-|>",color="b")
ax.add_artist(a)
ax.add_artist(a2)

B_shift = [0.75, 0, 0]

B1=Arrow3D([B_shift[0],b_save[0][0]+B_shift[0]],[0,b_save[0][1]],[0,b_save[0][2]],mutation_scale=15,arrowstyle="Fancy",color="orange")
ax.add_artist(B1)

ax.legend([a, a2, B1], ['P(t)', 'P(t) perfect', 'B(t)'])

def update(num):
    line.set_xdata(sx_save[:num])
    line.set_ydata(sy_save[:num])
    line.set_3d_properties(sz_save[:num])
    line2.set_xdata(sx_perfect[:num])
    line2.set_ydata(sy_perfect[:num])
    line2.set_3d_properties(sz_perfect[:num])

    a._verts3d=[0,sx_save[num]],[0,sy_save[num]],[0,sz_save[num]]
    a2._verts3d=[0,sx_perfect[num]],[0,sy_perfect[num]],[0,sz_perfect[num]]
    B1._verts3d=[B_shift[0],b_save[num][0]+B_shift[0]],[0,b_save[num][1]],[0,b_save[num][2]]
    return [line,line2,a2,a,B1]

ani=animation.FuncAnimation(fig,update,frames=i,interval=1,blit=False,repeat=True)
ani.save('asr_rk.mp4')
plt.show()
