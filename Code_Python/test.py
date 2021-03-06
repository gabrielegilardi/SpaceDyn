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
import utils as utils

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
foot = element.base(name=name, mass=mass, inertia=inertia, cc=cc)
# foot = element.base(name=name, mass=mass, inertia=inertia)

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

# System
name = 'simple robot'

ee = {
    0: (3, [0.4, -0.05, 0.5], [0.0,  0.0, -pi/4], 'I'),
    1: (4, [0.15, 0.0,  1.0], [0.0, -pi/3, pi/2], 'I'),
    2: (0, [0.0,  0.0,  0.0], [0.0,  0.0,  0.0], 'I')
    }
param = {
    'solver': 'nb',
    'reps': 10,
    }
load = 'full_system'
bodies = [foot, leg, trunk, upperArm, lowerArm]
robot = element.model(name=name, bodies=bodies, ee=ee, load=load, param=param)
R0 = [1.0, 2.0, 3.0]
Q0 = [0.1, 0.2, 0.3]
v0 = [1.3, -2.2, 3.7]
w0 = [-0.2, 0.5, -0.4]
q = [0.1, 0.5, 1.0, 1.5]
qd = [-0.1, 0.2, -0.3, 0.4]
robot.set_init(R0=R0, Q0=Q0, v0=v0, w0=w0, q=q, qd=qd)

# ee = {
#     0: (0, [0.0, 0.0, 0.0], [0.0,  0.0, 0.0], 'I'),
#     }
# bodies = [foot]
# load = 'base_only'
# param = {
#     'solver': 'nb',
#     'reps': 10,
#     }
# robot = element.model(name=name, bodies=bodies, ee=ee, load=load, param=param)
# robot.set_init()


# Simulate
ts = 0.0
tf = 1.0
dt = 0.01
rec = 0.1
robot.simulate(ts=ts, tf=tf, dt=dt, rec=rec)

aa = robot.res
plt.subplot(221)
plt.plot(aa[:, 0], aa[:, 1])
plt.grid(b=True)
plt.subplot(222)
plt.plot(aa[:, 0], aa[:, 2])
plt.grid(b=True)
plt.subplot(223)
plt.plot(aa[:, 0], aa[:, 3])
plt.grid(b=True)
plt.subplot(224)
plt.plot(aa[:, 0], aa[:, 4])
plt.grid(b=True)
plt.show()

print('time = ', aa[-1, 0])
print('R0[2] = ', aa[-1, 1])
print('q[2] = ', aa[-1, 2])
print('Q[0] = ', aa[-1, 3])
print('qd[0] = ', aa[-1, 4])
print('CoM acc = ', robot.calc_CoM()[2])
print('Kin Ene = ', robot.calc_kin_ener()[0])
print('Pot Ene = ', robot.calc_pot_ener()[0])
print('Lin Mom = ', robot.calc_lin_mom()[0])
print('Ang Mom = ', robot.calc_ang_mom1()[0])

tot = aa[-1, 0] + aa[-1, 1] + aa[-1, 2] + aa[-1, 3] + aa[-1, 4] \
      + (robot.calc_CoM()[2]).sum() + robot.calc_kin_ener()[0] \
      + robot.calc_pot_ener()[0] + (robot.calc_lin_mom()[0]).sum() \
      + (robot.calc_ang_mom1()[0]).sum()
print('tot = ', tot)

# robot.info()
