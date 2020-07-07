"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np

import elements as element


d2r = 180.0 / np.pi

# Links numbering starts from zero (the base)
# --------------
# 0     Foot/base
# 1     Leg
# 2     Trunk
# 3     Upper arm
# 4     Lower arm

# 0 - Foot (base)
name = 'foot/base'
mass = 2.3
inertia = 1.5 * np.eye(3)
cc = {1: [-0.15, 0.05, 0.0]}      # Leg connection
foot = element.base(name=name, mass=mass, inertia=inertia, cc=cc)

# 1 - Leg
name = 'leg'
mass = 0.4
inertia = 0.75 * np.eye(3)
j_type = 'R'
Qi = [0.0, 0.0, 30.0 / d2r]
cc = {1: [-0.25, 0.0, 0.0],      # Foot/base connection
      2: [ 0.20, 0.0, 0.0],      # Trunk connection
      4: [ 0.20, 0.0, 0.0]}      # Lower arm connection
leg = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                   Qi=Qi, cc=cc)

# 2 - Trunk
name = 'trunk'
mass = 1.4
inertia = 0.5 * np.eye(3)
j_type = 'R'
Qi = [0.0, 0.0, -30.0 / d2r]
cc = {2: [-0.2, 0.0, 0.0],      # Leg connection
      3: [ 0.3, 0.0, 0.0]}      # Upper arm connection
trunk = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                     Qi=Qi, cc=cc)

# 3 - Upper arm
name = 'upper arm'
mass = 0.7
inertia = 2.0 * np.eye(3)
j_type = 'R'
cc = {3: [-0.15, 0.0, 0.0]}      # Trunk connection
upperArm = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                        cc=cc)

# 4 - Lower arm
name = 'lower arm'
mass = 0.8
inertia = 4.0 * np.eye(3)
j_type = 'R'
cc = {4: [-0.23, 0.0, 0.0]}      # Leg connection
lowerArm = element.link(name=name, mass=mass, inertia=inertia, j_type=j_type,
                        cc=cc)

name = 'simple robot'
ee = {3: [0.13, 0.0, 0.0, 0.0, 90.0 / d2r, 90.0 / d2r],
      4: [0.17, 0.0, 0.0, -90.0 / d2r, 90.0 / d2r, 0.0],
      0: [0.4, -0.05, 0.0, 0.0, 0.0, -90.0 / d2r]}
bodies = [foot, leg, trunk, upperArm, lowerArm]
robot = element.model(name=name, bodies=bodies, ee=ee)

robot.info()

# Initial conditions
# R0 = np.array([1.0, 2.0, 3.0])
# Q0 = np.array([0.1, 0.2, 0.3])
# v0 = np.array([1.0, 2.0, 3.0])
# w0 = np.array([0.1, 0.2, 0.3])
# q = np.array([0.1, 0.5, 1.0, 1.5])
# qd = np.array([ 0.1, 0.2, 0.3, 0.4])


# robot.set_init(R0=R0, Q0=Q0, v0=v0, w0=w0, q=q, qd=qd)

# TK, TKi = robot.calc_kin_ener()
# VG, VGi = robot.calc_pot_ener()

# LM, LMi = robot.calc_lin_mom()
# LM1, LM1i = robot.calc_lin_mom1()

# HM, HMi = robot.calc_ang_mom(np.array([1.0, 2.0, 3.0]))
# HM1, HM1i = robot.calc_ang_mom1(np.array([1.0, 2.0, 3.0]), np.array([0.1, 0.2, 0.3]))
