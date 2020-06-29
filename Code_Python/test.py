"""
SpaceDyn - A Toolbox for Space and Mobile Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
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
foot = element.base(name=name, mass=mass, inertia=inertia, cc=cc)

# 1 - Leg
name = 'leg'
mass = 0.4
inertia = 0.75 * np.eye(3)
cc = {1: [-0.2, 0.0, 0.0],      # Foot/base connection
      2: [ 0.2, 0.0, 0.0],      # Trunk connection
      4: [ 0.2, 0.0, 0.0]}      # Lower arm connection
Qe = [0.1, 0.2, 0.3]            # End-point relative orientation
leg = element.link(name=name, mass=mass, inertia=inertia, cc=cc, Qe=Qe)

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
ce = [0.2, 0.0, 0.0]            # End-point position
upperArm = element.link(name=name, mass=mass, inertia=inertia, cc=cc, ce=ce)

# 4 - Lower arm
name = 'lower arm'
mass = 0.8
inertia = 4.0 * np.eye(3)
cc = {4: [-0.2, 0.0, 0.0]}      # Leg connection
ce = [0.2, 0.0, 0.0]            # End-point position
lowerArm = element.link(name=name, mass=mass, inertia=inertia, cc=cc, ce=ce)

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

# Links and joints must be listed in order
single_body = False
if (single_body):
    pass
    # foot = element.base(2.3, 1.5*np.eye(3), name='foot')
    # robot = mc.model([foot])
    # robot.set_init(R0=[1.0, 2.0, 3.0], Q0=[0.1, 0.2, 0.3],
    #                    v0=[1.0, 2.0, 3.0], w0=[0.1, 0.2, 0.3])

else:
    name = 'humanoid robot'
    bodies = [foot, leg, trunk, upperArm, lowerArm]
    joints = [j_foot_leg, j_leg_trunk, j_leg_lowerArm, j_trunk_upperArm]
    human_robot = model.model(name=name, bodies=bodies, joints=joints)
    # robot.set_init(R0=[1.0, 2.0, 3.0], Q0=[0.1, 0.2, 0.3],
    #                v0=[1.0, 2.0, 3.0], w0=[0.1, 0.2, 0.3],
    #                q=[0.1, 0.5, 1.0, 1.5], qd=[ 0.01, 0.02, 0.03, 0.04 ])

# robot.update_elements()
# robot.set_param(Gravity=[-9.81, 0.0, 0.7], Ez=[0.0, 1.0, 1.0])

# R0=np.array([[ 1.7,  1.7], [-2.0, -2.0 ], [3.4,  3.4 ]])
# Q0=np.array([[-0.1, -0.1], [ 0.23, 0.23], [0.76, 0.76]])
# print(R0, type(R0))
# print(Q0, type(Q0))
# dd = R0 - np.reshape(Q0[:, 0], (3, 1))
# print(dd, dd.ndim)
