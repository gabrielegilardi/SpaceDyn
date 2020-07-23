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
import dynamics as dyn
import utils as utils
import user as user

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
mass = 2.35
inertia = mass * utils.inertia('none')
cc = {1: [-0.15, 0.05, 0.0]}      # Leg connection
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
ee = {
      0: (3, [0.4, -0.05, 0.5], [0.0,  0.0, -pi/4]),
      1: (4, [0.15, 0.0,  1.0], [0.0, -pi/3, pi/2]),
     }

# System
name = 'simple robot'
bodies = [foot, leg, trunk, upperArm, lowerArm]
robot = element.model(name=name, bodies=bodies, ee=ee)

# Initial conditions
R0 = np.array([1.0, 2.0, 3.0])
Q0 = np.array([0.1, 0.2, 0.3])
A0 = utils.rpy2dc(Q0).T
v0 = np.array([1.3, -2.2, 3.7])
w0 = np.array([-0.2, 0.5, -0.4])
q = np.array([0.1, 0.5, 1.0, 1.5])
qd = np.array([-0.1, 0.2, -0.3, 0.4])
robot.set_init(R0=R0, A0=A0, v0=v0, w0=w0, q=q, qd=qd)

# ------- Test ------- #

if (len(sys.argv) != 2):
    print("Usage: python test.py <block-to-test>")
    sys.exit(1)
test = sys.argv[1]

if (test == 'Jac_endpoint'):

    print('\nJacobian due to the end point')
    # ie =  0 body =  3
    # [[ 0.19866933 -0.84606681 -0.34695908  0.        ]
    #  [-0.0978434  -0.33038082 -0.35348895  0.        ]
    #  [ 0.97517033 -0.26317159 -0.24426412  0.        ]
    #  [ 0.         -0.05657752 -0.05657752  0.        ]
    #  [ 0.         -0.52939695 -0.52939695  0.        ]
    #  [ 0.          0.84648559  0.84648559  0.        ]]
    # ie =  1 body =  4
    # [[ 0.19866933  0.          0.          0.19866933]
    #  [-0.0978434   0.          0.         -0.0978434 ]
    #  [ 0.97517033  0.          0.          0.97517033]
    #  [ 0.          0.          0.          0.        ]
    #  [ 0.          0.          0.          0.        ]
    #  [ 0.          0.          0.          0.        ]]
    for ie in range(0, robot.num_e):
        seq_link = kin.j_num(robot.SE[ie], robot.BB)
        Jacobian = kin.calc_je(robot.RR, robot.AA, q, seq_link, robot.j_type,
                               robot.cc, robot.ce[:, ie], robot.Qe[:, ie])
        print('ie = ', ie, ', body = ', robot.SE[ie])
        print(Jacobian)

elif (test == 'rne'):

    print('\nData:')
    vd0 = np.array([-1.7, 2.4, -4.5])
    wd0 = np.array([0.3, -0.2, 0.13])
    qdd = np.array([0.1, -0.3, 0.6, -1.1])
    print('vd0 = ', vd0)
    print('wd0 = ', wd0)
    print('qdd = ', qdd)

    print('\nExternal forces')
    # [[-10.3  -12.36]
    #  [ 11.4   13.68]
    #  [ 20.4   24.48]]
    F0, T0, tau, Fe, Te = user.calc_forces(robot.num_j, robot.num_e)
    print(Fe)

    print('\nExternal moments')
    # [[ 2.2  -1.54]
    #  [-4.4   3.08]
    #  [ 1.6  -1.12]]
    print(Te)

    # <Force> has shape (6+num_j, )
    print('\nRNE result F0')
    Force = dyn.r_ne(robot.RR, robot.AA, v0, w0, q, qd, vd0, wd0, qdd, Fe, Te,
                     robot.SS, robot.SE, robot.BB, robot.j_type, robot.cc,
                     robot.ce,robot.mass, robot.inertia)
    # [ 12.07580754 -12.51535814 -14.9985643 ]
    print(Force[0:3])

    print('\nRNE result T0')
    # [ 24.47189078  39.44989961 -11.31933556]
    print(Force[3:6])

    print('\nRNE result tau')
    # [-21.82565933  -2.32047012   1.56963594 -17.20710116]
    print(Force[6:])

elif (test == 'f_dyn'):

    print('\nControl terms')
    # F0 =  [-1.7  2.4 -4.5], shape (3, )
    # T0 =  [ 0.3  -0.2   0.13], shape (3, )
    # tau =  [ 0.1 -0.3  0.6 -1.1], shape (num_j, )
    F0, T0, tau, Fe, Te = user.calc_forces(robot.num_j, robot.num_e)
    print('F0 = ', F0)
    print('T0 = ', T0)
    print('tau = ', tau)

    print('\nExternal forces')
    # [[-10.3  -12.36]
    #  [ 11.4   13.68]
    #  [ 20.4   24.48]]
    print(Fe)

    print('\nExternal moments')
    # [[ 2.2  -1.54]
    #  [-4.4   3.08]
    #  [ 1.6  -1.12]]
    print(Te)

    print('\nBase linear acceleration')
    # [ -2.00047845   1.68790646 -11.44146866], shape (3, )
    vd0, wd0, qdd = dyn.f_dyn(R0, A0, v0, w0, q, qd, F0, T0, tau, Fe, Te,
                              robot.SS, robot.SE, robot.BB, robot.j_type,
                              robot.cc, robot.ce, robot.mass, robot.inertia,
                              robot.Qi, robot.Qe)
    print(vd0)

    print('\nBase angular acceleration')
    # [-7.80976388 -8.68344374  0.92761215], shape (3, )
    print(wd0)

    print('\nJoint angular/linear accelerations')
    # [13.36313061   15.57472009 -109.54087316   16.22245975], shape (num_j, )
    print(qdd)

elif (test == 'matrix_HH'):

    print('\nMatrix HH shape')
    # (10, 10)
    HH = kin.calc_hh(robot.RR, robot.AA, robot.mass, robot.inertia, robot.BB,
                     robot.j_type, robot.cc)
    print(HH.shape)

    print('\nMatrix HH_bb')
    # [[ 5.79        0.          0.          0.          2.13252875 -1.90414099]
    #  [ 0.          5.79        0.         -2.13252875  0.          0.28314725]
    #  [ 0.          0.          5.79        1.90414099 -0.28314725  0.        ]
    #  [ 0.         -2.13252875  1.90414099  6.2071716  -0.13503105 -0.3638451 ]
    #  [ 2.13252875  0.         -0.28314725 -0.13503105  5.07870753 -1.08971231]
    #  [-1.90414099  0.28314725  0.         -0.3638451  -1.08971231  3.71116261]]
    print(HH[0:6, 0:6])

    print('\nMatrix HH_bq')
    # [[ 0.6834225  -0.70083883 -0.05998451  0.17482901]
    #  [-0.33658128 -0.04461712 -0.07428796 -0.08610219]
    #  [ 3.35458593 -0.07474657 -0.05046937  0.85814989]
    #  [ 2.06551565 -0.04429423 -0.01282562  0.51411802]
    #  [ 0.14755127 -0.27832932 -0.02809198  0.07226354]
    #  [-0.4059985   0.58144965  0.05659342 -0.09748961]]
    print(HH[0:6, 6:])

    print('\nMatrix HH_qq')
    # [[ 3.44       -0.20776033 -0.05386473  0.88      ]
    #  [-0.20776033  0.36715122  0.04537632  0.        ]
    #  [-0.05386473  0.04537632  0.0162      0.        ]
    #  [ 0.88        0.          0.          0.88      ]]
    print(HH[6:, 6:])

elif (test == 'acceleration'):

    print('\nData:')
    vd0 = np.array([-1.7, 2.4, -4.5])
    wd0 = np.array([0.3, -0.2, 0.13])
    qdd = np.array([0.1, -0.3, 0.6, -1.1])
    print('vd0 = ', vd0)
    print('wd0 = ', wd0)
    print('qdd = ', qdd)

    print('\nLinear acceleration')
    # [[-1.7        -1.82564463 -1.89752749 -1.86707295 -2.02768068]
    #  [ 2.4         2.25985011  2.16232466  1.99136753  1.67073172]
    #  [-4.5        -4.40821138 -4.34017092 -4.37590531 -5.88368248]]
    vd, wd = kin.calc_acc(robot.AA, robot.ww, vd0, wd0, q, qd, qdd, robot.BB,
                          robot.j_type, robot.cc)
    print(vd)

    print('\nAngular acceleration')
    # [[ 0.3         0.3         0.35927006  0.26187834  0.3       ]
    #  [-0.2        -0.2        -0.00279529 -0.3780119  -0.2       ]
    #  [ 0.13        0.13       -0.09711205  0.37052886  0.13      ]]
    print(wd)

elif (test == 'connectivity'):
    
    print('\nBodies upper connection(s) matrix')
    # [[-1  1  0  0  0]
    #  [ 0 -1  1  0  1]
    #  [ 0  0 -1  1  0]
    #  [ 0  0  0 -1  0]
    #  [ 0  0  0  0 -1]]
    SS = robot.SS
    print(SS)

    print('\nBodies endpoint(s) vector')
    # [0 3 3 4], shape (num_e, )
    SE = robot.SE
    print(SE)

    print('\nLinks lower connection vector')
    # [0 1 2 1], shape (num_j, )
    BB = robot.BB
    print(BB)

elif (test == 'properties'):

    print('\nMasses')
    # [2.35 0.43 1.41 0.72 0.88], shape (num_b, )
    print(robot.mass)

    print('\nInertia (transpose)')
    # [[2.35       0.         0.        ]
    #  [0.         2.35       0.        ]
    #  [0.         0.         2.35      ]
    #  [0.01110833 0.         0.        ]
    #  [0.         0.01110833 0.        ]
    #  [0.         0.         0.01075   ]
    #  [0.05619857 0.         0.        ]
    #  [0.         0.05619857 0.        ]
    #  [0.         0.         0.05619857]
    #  [0.015      0.         0.        ]
    #  [0.         0.015      0.        ]
    #  [0.         0.         0.        ]
    #  [0.02713333 0.         0.        ]
    #  [0.         0.00366667 0.        ]
    #  [0.         0.         0.02933333]]
    print(robot.inertia.T)

    print('\nJoint type')
    # ['P', 'R', 'R', 'P'], shape (num_j, )
    print(robot.j_type)

    print('\nBody-to-body connections')
    # cc[0, :, :] = X coordinate
    # [[[ 0.   -0.15  0.    0.    0.  ]
    #   [ 0.   -0.25  0.2   0.    0.2 ]
    #   [ 0.    0.   -0.2   0.3   0.  ]
    #   [ 0.    0.    0.   -0.15  0.  ]
    #   [ 0.    0.    0.    0.   -0.23]]
    #
    # cc[1, :, :] = Y coordinate
    #  [[ 0.    0.05  0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]]
    #
    # cc[2, :, :] = Z coordinate
    #  [[ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]
    #   [ 0.    0.    0.    0.    0.  ]]]
    print(robot.cc)

    print('\nJoint/centroid frames')
    # [[ 0.          0.          0.          0.        ]
    #  [ 0.         -0.52359878  0.          0.        ]
    #  [ 0.78539816  0.          0.          1.04719755]]
    print(robot.Qi)

    print('\nBody-to-endpoint connections')
    # [[ 0.4   0.15]
    #  [-0.05  0.  ]
    #  [ 0.5   1.  ]]
    print(robot.ce)

    print('\nEndpoint frames')
    # [[ 0.          0.        ]
    #  [ 0.         -1.04719755]
    #  [-0.78539816  1.57079633]]
    print(robot.Qe)

elif (test == 'state'):

    print('\nSystem rotation matrices (transpose)')
    # [[ 0.93629336  0.31299183 -0.15934508]
    #  [-0.28962948  0.94470249  0.153792  ]
    #  [ 0.19866933 -0.0978434   0.97517033]
    #  [ 0.45726042  0.88932418 -0.00392662]
    #  [-0.86685835  0.44668689  0.22142135]
    #  [ 0.19866933 -0.0978434   0.97517033]
    #  [ 0.01910228  0.84711447  0.53106702]
    #  [-0.99821545  0.04621626 -0.03781495]
    #  [-0.05657752 -0.52939695  0.84648559]
    #  [-0.82964834  0.49658754  0.25511655]
    #  [-0.55541212 -0.6878515  -0.46730899]
    #  [-0.05657752 -0.52939695  0.84648559]
    #  [-0.52209115  0.83150428  0.1897932 ]
    #  [-0.82942832 -0.54683388  0.11411123]
    #  [ 0.19866933 -0.0978434   0.97517033]]
    AA = kin.calc_aa(A0, q, robot.BB, robot.j_type, robot.Qi)
    print(AA.T)

    print('\nSystem positions')
    # [[1.         0.97925656 1.0745291  0.95581253 1.24863168]
    #  [2.         2.21283306 2.56012078 2.88874326 2.43517878]
    #  [3.         3.12812674 3.23355482 3.43114241 4.63374934]]
    RR = kin.calc_pos(R0, AA, q, robot.BB, robot.j_type, robot.cc)
    print(RR)

    print('\nSystem linear velocities')
    # [[ 1.3         1.42932966  1.58103017  1.75971121  2.35054698]
    #  [-2.2        -2.15629294 -2.17146769 -2.07137279 -2.00205582]
    #  [ 3.7         3.57028808  3.45168166  3.45005619  3.7811995 ]]
    vv, ww = kin.calc_vel(AA, v0, w0, q, qd, robot.BB, robot.j_type, robot.cc)
    print(vv)

    print('\nSystem angular velocities')
    # [[-0.2        -0.2        -0.2113155  -0.19434225 -0.2       ]
    #  [ 0.5         0.5         0.39412061  0.5529397   0.5       ]
    #  [-0.4        -0.4        -0.23070288 -0.48464856 -0.4       ]]
    print(ww)

elif (test == 'endpoint'):

    print('\nLink sequence')
    # body =  0 sequence =  []
    # body =  1 sequence =  [1]
    # body =  2 sequence =  [1, 2]
    # body =  3 sequence =  [1, 2, 3]
    # body =  4 sequence =  [1, 4]    
    # all with shape shape (num_links_in_sequence, )
    for i in range(robot.num_b):
        seq_link = kin.j_num(i, robot.BB)
        print('body = ', i, ', sequence = ', seq_link)

    print('\nPosition joints in link sequence')
    # [[0.84507452 1.07070864 1.08025978]
    #  [2.00028635 2.39069789 2.81425512]
    #  [3.03159136 3.12734141 3.39287493]]    
    seq_link = kin.j_num(3, robot.BB)
    POS_jnt, ORI_jnt = kin.f_kin_j(robot.RR, robot.AA, q, seq_link,
                                   robot.j_type, robot.cc)
    print(POS_jnt)

    print('\nOrientation joints in link sequence (transpose)')
    # [[ 0.45726042  0.88932418 -0.00392662]
    #  [-0.86685835  0.44668689  0.22142135]
    #  [ 0.19866933 -0.0978434   0.97517033]
    #  [ 0.01910228  0.84711447  0.53106702]
    #  [-0.99821545  0.04621626 -0.03781495]
    #  [-0.05657752 -0.52939695  0.84648559]
    #  [-0.82964834  0.49658754  0.25511655]
    #  [-0.55541212 -0.6878515  -0.46730899]
    #  [-0.05657752 -0.52939695  0.84648559]]    
    print(ORI_jnt.T)

    print('\nPosition/orientation endpoint in link sequence')
    # ie =  0 body =  3
    # [0.62343504 2.85707237 3.97979727], shape (3, )
    # [[-0.19391429 -0.97938564 -0.05657752]
    #  [ 0.83752487 -0.13524404 -0.52939695]
    #  [ 0.510832   -0.15004271  0.84648559]]
    # ie =  1 body =  4
    # [1.36898733 2.46206103 5.63738865], shape (3, )
    # [[-0.82942832  0.08899289  0.55147886]
    #  [-0.54683388 -0.33101728 -0.76902553]
    #  [ 0.11411123 -0.93941888  0.32321943]]    
    for ie in range(robot.num_e):
        seq_link = kin.j_num(robot.SE[ie], robot.BB)
        POS_e, ORI_e = kin.f_kin_e(robot.RR, robot.AA, seq_link, robot.Qe[:, ie],
                                robot.ce[:, ie])
        print('ie = ', ie, ', body = ', robot.SE[ie])
        print(POS_e)
        print(ORI_e)

elif (test == 'utils'):

    print('\nData:')
    roll = pi/6
    pitch = pi/3
    yaw = pi/4
    V1 = np.array([1.1, 2.2, 3.3])
    V2 = np.array([0.7, -0.2, -2.4])
    M1 = np.tile(V1.reshape(3, 1), (1, 5))
    M2 = np.tile(V2.reshape(3, 1), (1, 5))
    w = np.array([0.4, -0.7, 0.3])
    dt = 0.2
    print('roll= ', roll)
    print('pitch=', pitch)
    print('yaw= ', yaw)
    print('V1= ', V1)
    print('V2= ', V2)
    print('M1=')
    print(M1)
    print('M2=')
    print(M2)
    print('dt= ', dt)

    print('\nRotation about X')
    # [[ 1.         0.         0.       ]
    #  [ 0.         0.8660254  0.5      ]
    #  [ 0.        -0.5        0.8660254]]
    Rx = utils.rotX(roll)
    print(Rx)

    print('\nRotation about Y')
    # [[ 0.8660254  0.       -0.5      ]
    # [ 0.         1.         0.       ]
    # [ 0.5        0.         0.8660254]]
    Ry = utils.rotY(roll)
    print(Ry)

    print('\nRotation about Z')
    # [[ 0.8660254  0.5        0.       ]
    #  [-0.5        0.8660254  0.       ]
    #  [ 0.         0.         1.       ]]
    Rz = utils.rotZ(roll)
    print(Rz)

    print('\nRotation matrix from RPY angles')
    # [[ 0.35355339  0.91855865 -0.1767767 ]
    #  [-0.35355339  0.30618622  0.88388348]
    #  [ 0.8660254  -0.25        0.4330127 ]]
    R = utils.rpy2dc(np.array([roll, pitch, yaw]))
    print(R)

    print('\nRPY angles from rotation matrix (general case)')
    # [0.52359878 1.04719755 0.78539816], shape (3, )
    rpy = utils.dc2rpy(R)
    print(rpy)

    print('\nRPY angles from rotation matrix (singular case)')
    # [1.30899694 1.57079633 0.        ], shape (3, )
    R = utils.rpy2dc(np.array([roll, pi/2, yaw]))
    rpy = utils.dc2rpy(R)
    print(rpy)

    print('\nRotation matrix from Euler angles')
    # [[ 0.43559574  0.65973961  0.61237244]
    #  [-0.78914913 -0.04736717  0.61237244]
    #  [ 0.4330127  -0.75        0.5       ]]
    R = utils.eul2dc(np.array([roll, pitch, yaw]))
    print(R)

    print('\nEuler angles from rotation matrix (general case)')
    # [0.52359878 1.04719755 0.78539816], shape (3, )
    eul = utils.dc2eul(R)
    print(eul)

    print('\nEuler angles from rotation matrix (singular case)')
    # [1.30899694 0.         0.        ], shape (3, )
    R = utils.eul2dc(np.array([roll, 0.0, yaw]))
    eul = utils.dc2eul(R)
    print(eul)

    print('\nSkew-symmetric matrix from vector')
    # [[ 0.         -0.78539816  1.04719755]
    #  [ 0.78539816  0.         -0.52359878]
    #  [-1.04719755  0.52359878  0.        ]]
    B_skew = utils.tilde(np.array([roll, pitch, yaw]))
    print(B_skew)

    print('\nCross product (two vectors)')
    # [-4.62  4.95 -1.76], shape (3, )
    V12 = utils.cross(V1, V2)
    print(V12)

    print('\nCross product (two matrices)')
    # [[-4.62 -4.62 -4.62 -4.62 -4.62]
    #  [ 4.95  4.95  4.95  4.95  4.95]
    #  [-1.76 -1.76 -1.76 -1.76 -1.76]]
    M12 = utils.cross(M1, M2)
    print(M12)

    print('\nRotation matrix for a rotation about a vector')
    # [[ 0.98842859 -0.06529064 -0.13691627]
    #  [ 0.05411824  0.99501232 -0.08379557]
    #  [ 0.14170444  0.07541627  0.98703204]]
    R = utils.rotW(w, dt)
    print(R)
