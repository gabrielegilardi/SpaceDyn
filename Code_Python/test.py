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

import elements as element
import kinematics as kin
import utils as utils

pi = np.pi
d2r = 180.0 / pi

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
Qi = [0.0, 0.0, pi/6]
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
Qi = [0.0, 0.0, -pi/6]
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
ee = {
      0: (0, [0.4, -0.05, 0.0], [ 0.0,  0.0, -pi/2]),
      2: (3, [0.13, 0.0,  0.0], [ 0.0,  pi/2, pi/2]),
      1: (3, [0.15, 0.0,  0.0], [ 0.0, -pi/2, pi/2]),
      3: (4, [0.17, 0.0,  0.0], [-pi/2, pi/2, 0.0])
     }
bodies = [foot, leg, trunk, upperArm, lowerArm]
robot = element.model(name=name, bodies=bodies, ee=ee)

# Initial conditions
R0 = np.array([1.0, 2.0, 3.0])
Q0 = np.array([0.1, 0.2, 0.3])
v0 = np.array([1.0, 2.0, 3.0])
w0 = np.array([0.1, 0.2, 0.3])
q = np.array([0.1, 0.5, 1.0, 1.5])
qd = np.array([ 0.1, 0.2, 0.3, 0.4])
robot.set_init(R0=R0, Q0=Q0, v0=v0, w0=w0, q=q, qd=qd)

# ------- Test ------- #

if len(sys.argv) != 2:
    print("Usage: python test.py <block-to-test>")
    sys.exit(1)
test = sys.argv[1]

if (test == 'utils'):

    print('Data:')
    roll = pi/6     # 30 deg
    pitch = pi/3    # 60 deg
    yaw = pi/4      # 45 deg
    phi = pi/6      # 30 deg
    theta = pi/3    # 60 deg
    psi = pi/4      # 45 deg
    print('\nroll= ', roll, 'pitch=', pitch, 'yaw= ', yaw)
    print('\nphi= ', phi, 'theta=', theta, 'psi= ', psi)

    print('\nSimple rotation in X')
    # [[ 1.         0.         0.       ]
    #  [ 0.         0.8660254  0.5      ]
    #  [ 0.        -0.5        0.8660254]]
    rotX = utils.rotX(roll)
    print(rotX)

    print('\nSimple rotation in Y')
    # [[ 0.8660254  0.       -0.5      ]
    # [ 0.         1.         0.       ]
    # [ 0.5        0.         0.8660254]]
    rotY = utils.rotY(roll)
    print(rotY)

    print('\nSimple rotation in Z')
    # [[ 0.8660254  0.5        0.       ]
    #  [-0.5        0.8660254  0.       ]
    #  [ 0.         0.         1.       ]]    
    rotZ = utils.rotZ(roll)
    print(rotZ)

    print('\nRotation matrix from RPY angles (array)')
    # [[ 0.35355339  0.91855865 -0.1767767 ]
    #  [-0.35355339  0.30618622  0.88388348]
    #  [ 0.8660254  -0.25        0.4330127 ]]    
    Cmat = utils.rpy2dc(np.array([roll, pitch, yaw]))
    print(Cmat)

    print('\nRotation matrix from RPY angles (scalars)')
    # [[ 0.35355339  0.91855865 -0.1767767 ]
    #  [-0.35355339  0.30618622  0.88388348]
    #  [ 0.8660254  -0.25        0.4330127 ]]    
    Cmat = utils.rpy2dc(roll, pitch, yaw)
    print(Cmat)

    print('\nRPY angles from rotation matrix (general case)')
    # [0.52359878 1.04719755 0.78539816]
    rpy = utils.dc2rpy(Cmat)
    print(rpy)

    print('\nRPY angles from rotation matrix (singular case)')
    # [1.30899694 1.57079633 0.        ]
    Cmat = utils.rpy2dc(roll, pi/2, yaw)
    rpy = utils.dc2rpy(Cmat)
    print(rpy)

    print('\nRotation matrix from Euler angles (array)')
    # [[ 0.43559574  0.65973961  0.61237244]
    #  [-0.78914913 -0.04736717  0.61237244]
    #  [ 0.4330127  -0.75        0.5       ]]
    Cmat = utils.eul2dc(np.array([phi, theta, psi]))
    print(Cmat)

    print('\nRotation matrix from Euler angles (scalars)')
    # [[ 0.43559574  0.65973961  0.61237244]
    #  [-0.78914913 -0.04736717  0.61237244]
    #  [ 0.4330127  -0.75        0.5       ]]
    Cmat = utils.eul2dc(phi, theta, psi)
    print(Cmat)

    print('\nEuler angles from rotation matrix (general case)')
    # [0.52359878 1.04719755 0.78539816]
    eul = utils.dc2eul(Cmat)
    print(eul)

    print('\nEuler angles from rotation matrix (singular case)')
    # [1.30899694 0.         0.        ]
    Cmat = utils.eul2dc(phi, 0.0, psi)
    eul = utils.dc2eul(Cmat)
    print(eul)

    print('\nSkew-symmetric matrix from vactor')
    # [[ 0.         -0.78539816  1.04719755]
    #  [ 0.78539816  0.         -0.52359878]
    #  [-1.04719755  0.52359878  0.        ]]
    Bmat = utils.tilde(np.array([roll, pitch, yaw]))
    print(Bmat)

