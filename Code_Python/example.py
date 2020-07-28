"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

import elements as element
import kinematics as kin
import dynamics as dyn
import utils as utils
import user as user


def plot_res(results):
    """
    """
    for result in results:

        result = result.flatten()
        plt.plot(result)

    plt.grid(b=True)

pi = np.pi

# Links numbering starts from zero (the base)
# --------------
# 0     Foot/base
# 1     Leg
# 2     Trunk
# 3     Upper arm
# 4     Lower arm

# 0 - Foot (base)
name = 'foot/base'
mass = 2.35
inertia = mass * utils.inertia('none')
cc = {1: [-0.15, 0.05, 0.0]}      # Leg connection
cc = {}
foot = element.base(name=name, mass=mass, inertia=inertia, cc=cc)

# 1 - Leg
name = 'leg'
mass = 0.43
inertia = mass * utils.inertia('cylinder', 0.2, 0.1, 0.4)
j_type = 'P'
Qi = [0.0, 0.0, pi/4]
cc = {1: [-0.25, 0.0, 0.0],      # Foot/base connection
      2: [ 0.20, 0.0, 0.0],      # Trunk connection
      4: [ 0.20, 0.0, 0.0]}      # Lower arm connection
leg = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                   Qi=Qi, cc=cc)

# 2 - Trunk
name = 'trunk'
mass = 1.41
inertia = mass * utils.inertia('sphere', 0.3, 0.15)
j_type = 'R'
Qi = [0.0, -pi/6, 0.0]
cc = {2: [-0.2, 0.0, 0.0],      # Leg connection
      3: [ 0.3, 0.0, 0.0]}      # Upper arm connection
trunk = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                     Qi=Qi, cc=cc)

# 3 - Upper arm
name = 'upper arm'
mass = 0.72
inertia = mass * utils.inertia('bar', 0.5)
j_type = 'R'
cc = {3: [-0.15, 0.0, 0.0]}      # Trunk connection
upperArm = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                        cc=cc)

# 4 - Lower arm
name = 'lower arm'
mass = 0.88
inertia = mass * utils.inertia('prism', 0.1, 0.2, 0.6)
j_type = 'P'
Qi = [0.0, 0.0, pi/3]
cc = {4: [-0.23, 0.0, 0.0]}      # Leg connection
lowerArm = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                        Qi=Qi, cc=cc)

# Endpoints
# ee = {
#     0: (3, [0.4, -0.05, 0.5], [0.0,  0.0, -pi/4]),
#     1: (4, [0.15, 0.0,  1.0], [0.0, -pi/3, pi/2]),
#     }

ee = {
    0: (0, [0.0, +1.0, 0.0], [0.0,  pi/2, 0.0]),
    1: (0, [0.0, -1.0, 0.0], [0.0,  pi/2, 0.0]),
    }

# System
name = 'simple robot'
# bodies = [foot, leg, trunk, upperArm, lowerArm]
# robot = element.model(name=name, bodies=bodies, ee=ee, load='example')
bodies = [foot]
robot = element.model(name=name, bodies=bodies, ee=ee, load='example')

# Initial conditions
t0 = 0.0
R0 = np.array([0.0, 0.0, 0.0])
Q0 = np.array([0.0, 0.0, 0.0])
A0 = utils.rpy2dc(Q0).T
v0 = np.array([0.0, 0.0, 0.0])
w0 = np.array([0.0, 0.0, 0.0])
# q = np.array([0.1, 0.5, 1.0, 1.5])
# qd = np.array([-0.1, 0.2, -0.3, 0.4])
q = np.array([])
qd = np.array([])
robot.set_init(t0=t0, R0=R0, A0=A0, v0=v0, w0=w0, q=q, qd=qd)

dt = 0.01
n_steps = 250
pos = np.zeros((3, n_steps))
apos = np.zeros((3, n_steps))
vel = np.zeros((3, n_steps))
avel = np.zeros((3, n_steps))
acc = np.zeros((3, n_steps))
aacc = np.zeros((3, n_steps))
pos[:, 0] = R0
apos[:, 0] = Q0
vel[:, 0] = v0
avel[:, 0] = w0
acc[:, 0] = robot.vd[:, 0]
aacc[:, 0] = robot.wd[:, 0]

for i in range(n_steps):
    t = t0 + float(i) * dt
    # print('time = ', t)
    Fe, Te, tau = user.calc_forces(t, robot.num_j, robot.num_e, load='example')
    R0, A0, v0, w0, q, qd = \
        dyn.f_dyn_nb(dt, R0, A0, v0, w0, q, qd, Fe, Te, tau, robot.SS, robot.SE,
                     robot.BB, robot.j_type, robot.cc, robot.ce, robot.mass,
                     robot.inertia, robot.Qi, robot.Qe)
        
    # Save results
    vd0, wd0, qdd = dyn.f_dyn(R0, A0, v0, w0, q, qd, Fe, Te, tau, robot.SS,
                              robot.SE, robot.BB, robot.j_type, robot.cc,
                              robot.ce, robot.mass, robot.inertia, robot.Qi,
                              robot.Qe)
    pos[:, i] = R0
    apos[:, i] = utils.dc2rpy(A0.T)
    vel[:, i] = v0
    avel[:, i] = w0
    acc[:, i] = vd0
    aacc[:, i] = wd0


print('R0 = ', R0)
print('Q0 = ', utils.dc2rpy(A0.T))
print('v0 = ', v0)
print('w0 = ', w0)
print('vd0 = ', vd0)
print('wd0 = ', wd0)

plt.subplot(231)
plot_res(pos)
plt.subplot(232)
plot_res(vel)
plt.subplot(233)
plot_res(acc)

plt.subplot(234)
plot_res(apos)
plt.subplot(235)
plot_res(avel)
plt.subplot(236)
plot_res(aacc)
plt.show()

# plt.subplot(221)
# plt.plot(pos[0, :], pos[2, :])
# plt.xlabel('X')
# plt.ylabel('Z')
# plt.grid(b=True)
# plt.subplot(222)
# plt.plot(pos[0, :], pos[1, :])
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.grid(b=True)
# plt.subplot(223)
# plt.plot(pos[1, :], pos[2, :])
# plt.xlabel('Y')
# plt.ylabel('Z')
# plt.grid(b=True)
# plt.subplot(224, projection='3d')
# plt.plot(pos[0, :], pos[1, :], pos[2, :])

# plt.show()
