"""
SpaceDyn - A Toolbox for Space and Mobile Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.

- check arrays shape
- add initial data checks and cross-checks    
"""

import numpy as np

import elements as element
import model as model

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
cc = {1: [0.05, 0.0, 0.0]}      # Leg connection
foot = element.link(name=name, mass=mass, inertia=inertia, cc=cc)

# 1 - Leg
name = 'leg'
mass = 0.4
inertia = 0.75 * np.eye(3)
cc = {1: [-0.2, 0.0, 0.0],      # Foot/base connection
      2: [ 0.2, 0.0, 0.0],      # Trunk connection
      4: [ 0.2, 0.0, 0.0]}      # Lower arm connection
leg = element.link(name=name, mass=mass, inertia=inertia, cc=cc)

# 2 - Trunk
name = 'trunk'
mass = 1.4
inertia = 0.5 * np.eye(3)
cc = {2: [-0.2, 0.0, 0.0],      # Leg connection
      3: [ 0.2, 0.0, 0.0]}      # Upper arm connection
trunk = element.link(name=name, mass=mass, inertia=inertia, cc=cc)

# 3 - Upper arm
name = 'upper arm'
mass = 0.7
inertia = 2.0 * np.eye(3)
cc = {3: [-0.2, 0.0, 0.0]}      # Trunk connection
upperArm = element.link(name=name, mass=mass, inertia=inertia, cc=cc)

# 4 - Lower arm
name = 'lower arm'
mass = 0.8
inertia = 4.0 * np.eye(3)
cc = {4: [-0.2, 0.0, 0.0]}      # Leg connection
lowerArm = element.link(name=name, mass=mass, inertia=inertia, cc=cc)

# Joints numbering starts from one
# ------
# 1     Foot - Leg (0-1)
# 2     Leg - Trunk (1-2)
# 3     Leg - Lower arm (1-4)
# 4     Trunk - Upper arm (2-3)

# Foot - Leg
name = 'foot - leg'
j_type = 'R'
Qi = [0.1, 0.2, 0.3]
j_foot_leg = element.joint(name=name, j_type=j_type, Qi=Qi)

# Leg - Trunk
name = 'leg - trunk'
j_type = 'R'
j_leg_trunk = element.joint(name=name, j_type=j_type)

# Leg - Lower arm
name = 'leg - lower arm'
j_type = 'R'
j_leg_lowerArm = element.joint(name=name, j_type=j_type)

# Trunk - Upper arm
name = 'trunk - upper arm'
j_type = 'R'
j_trunk_upperArm = element.joint(name=name, j_type=j_type)


# Model
# ----

name = 'humanoid robot'
ee = {3: [0.3, 0.0, 0.0, -0.1, 0.0, 0.0],
      4: [0.2, 0.0, 0.0, -0.2, 0.0, 0.0],
      0: [0.0, 0.0, 0.0,  0.0, 0.0, 0.0]}
bodies = [foot, leg, trunk, upperArm, lowerArm]
joints = [j_foot_leg, j_leg_trunk, j_leg_lowerArm, j_trunk_upperArm]
human_robot = model.model(name=name, bodies=bodies, joints=joints, ee=ee)
human_robot.info()

R0 = np.array([1.0, 2.0, 3.0])
Q0 = np.array([0.1, 0.2, 0.3])
v0 = np.array([1.0, 2.0, 3.0])
w0 = np.array([0.1, 0.2, 0.3])
q = np.array([0.1, 0.5, 1.0, 1.5])
qd = np.array([ 0.1, 0.2, 0.3, 0.4])
human_robot.set_init(R0=R0, Q0=Q0, v0=v0, w0=w0, q=q, qd=qd)

