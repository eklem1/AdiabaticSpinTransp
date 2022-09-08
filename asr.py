#!/usr/bin/python3

# The exact solution is graphed.  For the derivation of the solution,
# please see the handwritten notes.

from scipy.constants import *
from math import *

from matplotlib import pyplot as plt
import numpy as np
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
# I think I've corrected signs below so that either sign can be selected
# After cheating by doing RK, I fixed the signs again.

# Set the field and rate at which it should change, here.
B1=1 # uT
T1=1 # s, time to go around once

omegaL = gamma*B1
Omega = 2*pi/T1 # Hz

# Omega=2*pi/T1 # Hz
# At 1 uT, gamma*B1 = 2*pi*(30 Hz) = omegaL
# So pick an Omega that's a bit slower than that so we're fairly adiabatic
# I picked Omega=2*pi*(1 Hz) here.

a = np.sqrt(omegaL**2 + Omega**2)
cosbeta = -omegaL/a #the negative here is needed if gamma > 0
sinbeta = Omega/a

print(f"B1: {B1} uT, wL: {omegaL}, Omega: {Omega} Hz, a:{a}")

"""
conversion to Jeff's version:
    sintheta=sinbeta
    costheta=cosbeta
    omega=a
    omega1=Omega
"""

# The 100 in the line below should make sure we get about 100 points
# in each of the revolutions in the rotating frame.
# N=abs(int(omega/omega1))
N=abs(int(a/Omega*10))

print(N)

# def spin_rot(t): #not sure what this is for
#     sx=sintheta**2*np.cos(omega*t)+costheta**2
#     sy=sintheta*np.sin(omega*t)
#     sz=sintheta*costheta*(np.cos(omega*t)-1)
#     return sx,sy,sz

def spin(t): # P(t) in the non-rotating frame (F_UCN), as a function of a, beta and Omega

    sx=sinbeta*cosbeta*(1 - np.cos(a*t))

    sy=(sinbeta**2*np.cos(a*t)+cosbeta**2)*np.cos(Omega*t)+sinbeta*np.sin(a*t)*np.sin(Omega*t)

    sz=(sinbeta**2*np.cos(a*t)+cosbeta**2)*np.sin(Omega*t)-sinbeta*np.sin(a*t)*np.cos(Omega*t)
    return sx,sy,sz


def spin_perfect(t): # non-rotating frame, but assuming perfect adiabatic
    sx=0*t
    sy=np.cos(Omega*t)
    sz=np.sin(Omega*t)
    return sx,sy,sz


def b(t):
    bx=0*t
    by=B1*np.cos(Omega*t)
    bz=B1*np.sin(Omega*t)
    return [bx,by,bz]  # uT


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

t=np.linspace(0,T1,N)
sx,sy,sz=spin(t)
sx_perfect,sy_perfect,sz_perfect=spin_perfect(t)
b_save=np.array(b(t)).T

ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

fig.suptitle("Precession in a rotating magnetic field\n"+
    f"$B$: {B1} uT, $\omega_L$: {wL:.3}, $\Omega$: {Omega:.3} Hz, a:{a:.3}")

line,=ax.plot(sx[:1],sy[:1],sz[:1], color='tomato')
line2,=ax.plot(sx_perfect[:1],sy_perfect[:1],sz_perfect[:1], color='royalblue')
a=Arrow3D([0,sx[0]],[0,sy[0]],[0,sz[0]],mutation_scale=20,arrowstyle="-|>",color="r")
a2=Arrow3D([0,sx_perfect[0]],[0,sy_perfect[0]],[0,sz_perfect[0]],mutation_scale=20,arrowstyle="-|>",color="b")
ax.add_artist(a)
ax.add_artist(a2)

B_shift = [0.75, 0, 0]

B1_ar=Arrow3D([B_shift[0],b_save[0][0]+B_shift[0]],[0,b_save[0][1]],[0,b_save[0][2]],mutation_scale=15,arrowstyle="Fancy",color="orange")
ax.add_artist(B1_ar)

ax.legend([a, a2, B1_ar], ['P(t)', 'P(t) perfect', 'B(t)'])

def update(num):
    line.set_xdata(sx[:num])
    line.set_ydata(sy[:num])
    line.set_3d_properties(sz[:num])
    line2.set_xdata(sx_perfect[:num])
    line2.set_ydata(sy_perfect[:num])
    line2.set_3d_properties(sz_perfect[:num])
    a._verts3d=[0,sx[num]],[0,sy[num]],[0,sz[num]]
    a2._verts3d=[0,sx_perfect[num]],[0,sy_perfect[num]],[0,sz_perfect[num]]
    B1_ar._verts3d=[B_shift[0],b_save[num][0]+B_shift[0]],[0,b_save[num][1]],[0,b_save[num][2]]
    return [line,line2,a2,a,B1_ar]

ani=animation.FuncAnimation(fig,update,frames=N,interval=1,blit=False,repeat=True)
ani.save('asr.mp4')
plt.show()

# writergif = animation.PillowWriter(fps=30)
# ani.save('filename.gif',writer=writergif)