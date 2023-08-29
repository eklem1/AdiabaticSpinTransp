#!/usr/bin/python3

# The exact solution is graphed.  For the derivation of the solution,
# please see the handwritten notes.

# by Jeff
# Edited by Emma

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
B1=1.5 # uT
T1=.5# s, time to go around once

omegaL = gamma*B1
Omega = 2*pi/T1 # Hz

# Omega=2*pi/T1 # Hz
# At 1 uT, gamma*B1 = 2*pi*(30 Hz) = omegaL
# So pick an Omega that's a bit slower than that so we're fairly adiabatic
# I picked Omega=2*pi*(1 Hz) here.

omega_eff = np.sqrt(omegaL**2 + Omega**2)
cosbeta = -omegaL/omega_eff #the negative here is needed if gamma > 0
sinbeta = Omega/omega_eff

k = omegaL/Omega

print(f"B1: {B1} uT, wL: {omegaL}, Omega: {Omega} Hz, omega_eff:{omega_eff}, k: {k}")

"""
conversion to Jeff's version:
    sintheta=sinbeta
    costheta=cosbeta
    omega=omega_eff
    omega1=Omega
"""

# The 100 in the line below should make sure we get about 100 points
# in each of the revolutions in the rotating frame.
# N=abs(int(omega/omega1))
N=abs(int(omega_eff/Omega*10)) #40

print("N: ",N)

def spinR2(t): # P(t) in the rotating frame (F_rot2), as a function of a, beta and Omega

    sx=-sinbeta*np.cos(omega_eff*t)
    sy=cosbeta +0*t
    sz=-sinbeta*np.sin(omega_eff*t)
    return sx,sy,sz

def spinR(t): # P(t) in the single rotating frame (F_rot2), as a function of a, beta and Omega

    sx=sinbeta*cosbeta*(1 - np.cos(omega_eff*t))
    sy=sinbeta**2*np.cos(omega_eff*t)+cosbeta**2
    sz=-sinbeta*np.sin(omega_eff*t)
    return sx,sy,sz

def spin(t): # P(t) in the non-rotating frame (F_UCN), as a function of a, beta and Omega

    sx=sinbeta*cosbeta*(1 - np.cos(omega_eff*t))
    sy=(sinbeta**2*np.cos(omega_eff*t)+cosbeta**2)*np.cos(Omega*t)+sinbeta*np.sin(omega_eff*t)*np.sin(Omega*t)
    sz=(sinbeta**2*np.cos(omega_eff*t)+cosbeta**2)*np.sin(Omega*t)-sinbeta*np.sin(omega_eff*t)*np.cos(Omega*t)
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

#update of plot in time for only lab frame
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

    #saves a still of the video for a number frame number (magic_value)
    # magic_value=73
    # if num == magic_value: fig.savefig(f'fig{num}.pdf')

    return [line,line2,a2,a,B1_ar]

#update of plot in time for all frames
def update_multi(num):
    line.set_xdata(sx[:num])
    line.set_ydata(sy[:num])
    line.set_3d_properties(sz[:num])
    line2.set_xdata(sx_perfect[:num])
    line2.set_ydata(sy_perfect[:num])
    line2.set_3d_properties(sz_perfect[:num])
    a._verts3d=[0,sx[num]],[0,sy[num]],[0,sz[num]]
    a2._verts3d=[0,sx_perfect[num]],[0,sy_perfect[num]],[0,sz_perfect[num]]
    B1_ar._verts3d=[B_shift[0],b_save[num][0]+B_shift[0]],[0,b_save[num][1]],[0,b_save[num][2]]

    line_r.set_xdata(sx_r[:num])
    line_r.set_ydata(sy_r[:num])
    line_r.set_3d_properties(sz_r[:num])
    line_r2.set_xdata(sx_r2[:num])
    line_r2.set_ydata(sy_r2[:num])
    line_r2.set_3d_properties(sz_r2[:num])
    a_r._verts3d=[0,sx_r[num]],[0,sy_r[num]],[0,sz_r[num]]
    a_r2._verts3d=[0,sx_r2[num]],[0,sy_r2[num]],[0,sz_r2[num]]

    return [line,line2,a2,a,B1_ar, line_r, line_r2, a_r, a_r2]

#plotting
fig = plt.figure()

fig.suptitle("Precession in a rotating magnetic field \n"+
    f"$B$: {B1} $\mu$T, $\omega_L$: {omegaL:.3}, $\Omega$: {Omega:.3} Hz, "+"$\omega_{eff}$ "+f":{omega_eff:.3}, $k$: {k:.3}")

#change this to false if you only want to plot in the lab frame
multiplot = False

if multiplot:

    ax1 = fig.add_subplot(2, 2, 1, projection='3d')
    ax2 = fig.add_subplot(2, 2, 3, projection='3d')
    ax3 = fig.add_subplot(2, 2, 4, projection='3d')

    t=np.linspace(0,T1,N)
    sx_r2,sy_r2,sz_r2=spinR2(t)
    sx_r,sy_r,sz_r=spinR(t)
    sx,sy,sz=spin(t)
    sx_perfect,sy_perfect,sz_perfect=spin_perfect(t)
    b_save=np.array(b(t)).T

    ax1.set_xlim(-1,1)
    ax1.set_ylim(-1,1)
    ax1.set_zlim(-1,1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.view_init(elev=30., azim=60)

    ax2.set_xlim(-0.5,0.5)
    ax2.set_ylim(-0.5,0.5)
    ax2.set_zlim(-0.5,0.5)
    ax2.set_xlabel("x'")
    ax2.set_ylabel("y'")
    ax2.set_zlabel("z'")
    ax2.view_init(elev=30., azim=60)

    ax3.set_xlim(-0.5,0.5)
    ax3.set_ylim(-0.5,0.5)
    ax3.set_zlim(-0.5,0.5)
    ax3.set_xlabel('x"')
    ax3.set_ylabel('y"')
    ax3.set_zlabel('z"')
    ax3.view_init(elev=30., azim=60)

    line,=ax1.plot(sx[:1],sy[:1],sz[:1], color='tomato')
    line2,=ax1.plot(sx_perfect[:1],sy_perfect[:1],sz_perfect[:1], color='royalblue')
    a=Arrow3D([0,sx[0]],[0,sy[0]],[0,sz[0]],mutation_scale=20,arrowstyle="-|>",color="r")
    a2=Arrow3D([0,sx_perfect[0]],[0,sy_perfect[0]],[0,sz_perfect[0]],mutation_scale=20,arrowstyle="-|>",color="b")
    ax1.add_artist(a)
    ax1.add_artist(a2)

    line_r,=ax2.plot(sx_r[:1],sy_r[:1],sz_r[:1], color='tomato')
    a_r=Arrow3D([0,sx_r[0]],[0,sy_r[0]],[0,sz_r[0]],mutation_scale=20,arrowstyle="-|>",color="r")
    ax2.add_artist(a_r)

    line_r2,=ax3.plot(sx_r2[:1],sy_r2[:1],sz_r2[:1], color='tomato')
    a_r2=Arrow3D([0,sx_r2[0]],[0,sy_r2[0]],[0,sz_r2[0]],mutation_scale=20,arrowstyle="-|>",color="r")
    ax3.add_artist(a_r2)

    #a shift for the arrow showing the B field direction so it can be seen better
    B_shift = [0.75, 0, 0]

    B1_ar=Arrow3D([B_shift[0],b_save[0][0]+B_shift[0]],[0,b_save[0][1]],[0,b_save[0][2]],mutation_scale=15,arrowstyle="Fancy",color="orange")
    ax1.add_artist(B1_ar)

    ax1.legend([a, a2, B1_ar], ['P(t)', 'P(t) perfect', 'B(t)'], title="$F_{lab}$", bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax2.legend([a_r], ['P(t)'], title="$F_{rot}$")
    ax3.legend([a_r2], ['P(t)'], title="$F_{rot2}$")

    ani=animation.FuncAnimation(fig,update_multi,frames=N,interval=1,blit=False,repeat=True)
    ani.save('asr_multiFrame.mp4')

else:
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
    ax.view_init(elev=30., azim=60)

    line,=ax.plot(sx[:1],sy[:1],sz[:1], color='tomato')
    line2,=ax.plot(sx_perfect[:1],sy_perfect[:1],sz_perfect[:1], color='royalblue')
    a=Arrow3D([0,sx[0]],[0,sy[0]],[0,sz[0]],mutation_scale=20,arrowstyle="-|>",color="r")
    a2=Arrow3D([0,sx_perfect[0]],[0,sy_perfect[0]],[0,sz_perfect[0]],mutation_scale=20,arrowstyle="-|>",color="b")
    ax.add_artist(a)
    ax.add_artist(a2)

    #a shift for the arrow showing the B field direction so it can be seen better
    B_shift = [0.75, 0, 0]

    B1_ar=Arrow3D([B_shift[0],b_save[0][0]+B_shift[0]],[0,b_save[0][1]],[0,b_save[0][2]],mutation_scale=15,arrowstyle="Fancy",color="orange")
    ax.add_artist(B1_ar)

    ax.legend([a, a2, B1_ar], ['P(t)', '$P_{perfect}(t) $', 'B(t)'], title="$F_{lab}$", zzbbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

    ani=animation.FuncAnimation(fig,update,frames=N,interval=1,blit=False,repeat=True)

    ani.save('./animation2.gif', writer='imagemagick', fps=60)
    # ani.save('asr.mp4')

plt.show()