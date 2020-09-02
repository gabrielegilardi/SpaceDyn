"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np
from utils import rpy2dc, cross, tilde

# Joint axis direction
Ez = np.array([0.0, 0.0, 1.0])


def calc_je(ie, RR, AA, q, Conn, Prop):
    """
    Returns the Jacobian (6 x num_j) of the endpoint <ie> (eqs 3.25 and 3.26).
    """
    SE, BB, j_type = Conn[1], Conn[2], Conn[3]
    cc, ce = Prop[2], Prop[3]
    
    # Link sequence from the base (excluded) to the endpoint
    seq_link = []
    i = SE[0, ie]
    while (i > 0):
        seq_link.insert(0, i)
        i = BB[i-1]

    # Position of the endpoint
    j = seq_link[-1]
    A_I_j = AA[:, 3*j:3*(j+1)]
    POS_e = RR[:, j] + A_I_j @ ce[:, ie]

    # Position for all joints in the link sequence
    n_links = len(seq_link)
    POS_jnt = np.zeros((3, n_links))
    Jacobian = np.zeros((6, len(q)))
    for i in range(n_links):

        j = seq_link[i]             # Joint/link number in the sequence
        idxj = j - 1                # Index link/joint <j> in q, j_type
        A_I_j = AA[:, 3*j:3*(j+1)]
        Ez_I_j = A_I_j @ Ez

        # Position for a rotational joint
        if (j_type[idxj] == 'R'):
            POS_jnt[:, i] = RR[:, j] + A_I_j @ cc[:, j, j]
            JJ_te = cross(Ez_I_j, POS_e - POS_jnt[:, i])
            JJ_re = Ez_I_j

        # Prismatic joint (rotational component already set to zero)
        elif (j_type[idxj] == 'P'):
            POS_jnt[:, i] = RR[:, j] + A_I_j @ (cc[:, j, j] - Ez * q[idxj])
            JJ_te = Ez_I_j
            JJ_re = np.zeros(3)

        # Assemble the endpoint Jacobian in the global Jacobian
        Jacobian[0:3, j-1] = JJ_te
        Jacobian[3:6, j-1] = JJ_re

    return Jacobian


def calc_aa(Q0, q, Conn, Prop):
    """
    Returns the rotation matrices of all bodies with respect to the inertial
    frame.

    AA = [AA_0, AA_1, ... ]         (3, 3*num_b)

    Note: for the links the centroid frame and the corresponding joint frame
          are always parallel.
    """
    BB, j_type = Conn[2], Conn[3]
    Qi = Prop[4]

    num_j = len(q)                      # Number of joints/links
    num_b = num_j + 1                   # Number of bodies
    AA = np.zeros((3, 3*num_b))

    # Base
    AA[:, 0:3] = rpy2dc(Q0).T

    # i = current link, k = lower link connection
    for i in range(1, num_b):

        idxi = i - 1                # Index link/joint <i> in BB, j_type, Qi
        k = BB[idxi]                # Index lower link connection
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Rotational joint
        if (j_type[idxi] == 'R'):
            A_rel = rpy2dc(Qi[:, idxi] + Ez * q[idxi]).T

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            A_rel = rpy2dc(Qi[:, idxi]).T

        # Update rotation matrix (integral of eqs. 3.8 and 3.12)
        AA[:, 3*i:3*(i+1)] = A_I_k @ A_rel

    return AA


def calc_pos(R0, AA, q, Conn, Prop):
    """
    Returns the position of all body centroids with respect to the inertial
    frame (figure 3.5).

    RR = [RR_0, RR_1, ... ]         (3, num_b)
    """
    BB, j_type = Conn[2], Conn[3]
    cc = Prop[2]

    num_j = len(q)                      # Number of joints/links
    num_b = num_j + 1                   # Number of bodies
    RR = np.zeros((3, num_b))

    # Base
    RR[:, 0] = R0

    # i = current link, k = lower link connection
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, q
        k = BB[idxi]             # Index lower link connection

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Rotational joint (integral of eq. 3.9)
        if (j_type[idxi] == 'R'):
            RR[:, i] = RR[:, k] + A_I_k @ cc[:, k, i] - A_I_i @ cc[:, i, i]

        # Prismatic joint (integral of eq. 3.13)
        elif (j_type[idxi] == 'P'):
            RR[:, i] = RR[:, k] + A_I_k @ cc[:, k, i] \
                                + A_I_i @ (Ez * q[idxi] - cc[:, i, i])

    return RR


def calc_vel(AA, q, Yd, Conn, Prop):
    """
    Returns the linear velocity of all body centroids and the angular velocity
    of all body (eqs. 3.8-3.9 and 3.12-3.13, figure 3.5).

    vv = [vv_0, vv_1, ... ]         (3, num_b)
    ww = [ww_0, ww_1, ... ]         (3, num_b)
    """
    BB, j_type = Conn[2], Conn[3]
    cc = Prop[2]
    v0, w0, qd = Yd[0:3], Yd[3:6], Yd[6:]

    num_b = len(q) + 1                   # Number of bodies
    vv = np.zeros((3, num_b))
    ww = np.zeros((3, num_b))

    # Base
    vv[:, 0] = v0
    ww[:, 0] = w0

    # i = current link, k = lower link connection
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, q, qd
        k = BB[idxi]             # Index lower link connection

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Relative positions
        cc_I_i = A_I_i @ cc[:, i, i]
        cc_I_k = A_I_k @ cc[:, k, i]

        # Joint values
        q_I_i = A_I_i @ (Ez * q[idxi])
        qd_I_i = A_I_i @ (Ez * qd[idxi])

        # Rotational joint (eqs. 3.8-3.9)
        if (j_type[idxi] == 'R'):
            ww[:, i] = ww[:, k] + qd_I_i
            vv[:, i] = vv[:, k] + cross(ww[:, k], cc_I_k) \
                                - cross(ww[:, i], cc_I_i)

        # Prismatic joint (eqs. 3.12-3.13)
        elif (j_type[idxi] == 'P'):
            ww[:, i] = ww[:, k]
            vv[:, i] = vv[:, k] + cross(ww[:, k], cc_I_k) \
                                - cross(ww[:, i], cc_I_i) \
                                + cross(ww[:, i], q_I_i) \
                                + qd_I_i

    return vv, ww


def calc_acc(AA, ww, q, qd, Ydd, Conn, Prop):
    """
    Returns the linear acceleration of all body centroids and the angular
    acceleration of all body (eqs. 3.10-3.11 and 3.14-3.15, figure 3.5).

    vd = [vd_0, vd_1, ... ]         (3, num_b)
    wd = [wd_0, wd_1, ... ]         (3, num_b)
    """
    BB, j_type = Conn[2], Conn[3]
    cc = Prop[2]
    vd0, wd0, qdd = Ydd[0:3], Ydd[3:6], Ydd[6:]

    num_b = len(q) + 1                   # Number of bodies
    vd = np.zeros((3, num_b))
    wd = np.zeros((3, num_b))

    # Base
    vd[:, 0] = vd0
    wd[:, 0] = wd0

    # i = current link, k = lower link connection
    for i in range(1, num_b):

        idxi = i - 1        # Index link/joint <i> in BB, j_type, q, qd, qdd
        k = BB[idxi]        # Index lower link connection

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Relative positions
        cc_I_i = A_I_i @ cc[:, i, i]
        cc_I_k = A_I_k @ cc[:, k, i]

        # Joint values
        q_I_i = A_I_i @ (Ez * q[idxi])
        qd_I_i = A_I_i @ (Ez * qd[idxi])
        qdd_I_i = A_I_i @ (Ez * qdd[idxi])

        # Rotational joint (eqs. 3.10-3.11)
        if (j_type[idxi] == 'R'):
            wd[:, i] = wd[:, k] + cross(ww[:, i], qd_I_i) + qdd_I_i
            vd[:, i] = vd[:, k] + cross(wd[:, k], cc_I_k) \
                                + cross(ww[:, k], cross(ww[:, k], cc_I_k)) \
                                - cross(wd[:, i], cc_I_i) \
                                - cross(ww[:, i], cross(ww[:, i], cc_I_i))

        # Prismatic joint (eqs. 3.14-3.15)
        elif (j_type[idxi] == 'P'):
            wd[:, i] = wd[:, k]
            vd[:, i] = vd[:, k] + cross(wd[:, k], cc_I_k) \
                                + cross(ww[:, k], cross(ww[:, k], cc_I_k)) \
                                - cross(wd[:, i], cc_I_i) \
                                - cross(ww[:, i], cross(ww[:, i], cc_I_i)) \
                                + cross(wd[:, i], q_I_i) \
                                + cross(ww[:, i], cross(ww[:, i], q_I_i)) \
                                + 2.0 * cross(ww[:, i], qd_I_i) \
                                + qdd_I_i

    return vd, wd


def calc_jac(RR, AA, Conn, Prop):
    """
    Returns the translational Jacobians (3 x num_j x num_j) of all link
    centroids (equation 3.25-3.26).
    """
    BB, j_type = Conn[2], Conn[3]
    cc = Prop[2]

    num_j = len(j_type)                     # Number of joints/links
    JJ_t = np.zeros((3, num_j*num_j))
    JJ_r = np.zeros((3, num_j*num_j))

    # Loop over all links
    for i in range(1, num_j+1):

        j = i       # Initial index from link/joint i to base (itself)

        # Follow the branch until the base (j = 0) is reached
        while(j > 0):

            idxj = j - 1         # Index link/joint <j> in BB, j_type
            A_I_j = AA[:, 3*j:3*(j+1)]
            Ez_I_j = A_I_j @ Ez
            cc_I_j = A_I_j @ cc[:, j, j]
            col = (i-1)*num_j+(j-1)

            # Rotational joint
            if (j_type[idxj] == 'R'):
                JJ_t[:, col] = \
                    cross(Ez_I_j, RR[:, i] - (RR[:, j] + cc_I_j))
                JJ_r[:, col] = Ez_I_j

            # Prismatic joint
            elif (j_type[idxj] == 'P'):
                JJ_t[:, col] = Ez_I_j
                JJ_r[:, col] = np.zeros(3)

            j = BB[idxj]     # Previous link/joint along the branch

    return JJ_t, JJ_r


def calc_hh(RR, AA, Conn, Prop):
    """
    Returns the system inertia matrix HH with shape (6 + num_j) x (6 + num_j)
    (equations 3.19-3.24).

    Subscripts: b = base, q = joints, w = base orientation, r = base position
    """
    mass, inertia = Prop[0], Prop[1]

    num_b = len(mass)               # Number of bodies
    num_j = num_b - 1               # Number of joints/links

    # System CoM with respect to the base (eq. 3.27)
    mt = mass.sum()
    Rg = (mass * RR[:, 0:num_b]).sum(axis=1) / mt
    r0g = Rg - RR[:, 0]

    # Translational and rotational jacobian of all centroids - shape of each
    # is (3 x num_j x num_j)
    JJ_t, JJ_r = calc_jac(RR, AA, Conn, Prop)       # Eq. 3.25-3.26

    # Initialize
    HH_ww = AA[:, 0:3] @ inertia[:, 0:3] @ AA[:, 0:3].T         # Eq. 3.21
    HH_wq = np.zeros((3, num_j))                                # Eq. 3.22
    HH_qq = np.zeros((num_j, num_j))                            # Eq. 3.23
    HH_rq = np.zeros((3, num_j))                                # Eq. 3.24

    # Assemble links
    for i in range(1, num_j+1):

        # Start/end index for link <i> in JJ_t and JJ_r
        id1 = (i - 1) * num_j
        id2 = i * num_j

        # Position link centroid <i> wrt the base centroid (eq. 3.28)
        r0i = RR[:, i] - RR[:, 0]

        # Link rotation and inertia matrix
        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I_i = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T

        # Eq. 3.21 - shape (3 x 3)
        HH_ww += In_I_i + mass[i] * tilde(r0i).T @ tilde(r0i)

        # Eq. 3.22 - shape (3 x num_j)
        HH_wq += In_I_i @ JJ_r[:, id1:id2] \
                 + mass[i] * tilde(r0i) @ JJ_t[:, id1:id2]

        # Eq. 3.23 - shape (num_j x num_j)
        HH_qq += JJ_r[:, id1:id2].T @ In_I_i @ JJ_r[:, id1:id2] \
                + mass[i] * JJ_t[:, id1:id2].T @ JJ_t[:, id1:id2]

        # Eq. 3.24 - shape (3 x num_j)
        HH_rq += mass[i] * JJ_t[:, id1:id2]

    # Matrix HH_bb (eq. 3.19) - shape (6 x 6)
    HH_bb = np.block([[mt * np.eye(3),  mt * tilde(r0g).T],
                      [mt * tilde(r0g), HH_ww            ]])

    # Matrix HH_bm (eq. 3.20) - shape (6 x num_j)
    HH_bq = np.block([[HH_rq],
                      [HH_wq]])

    # Matrix HH (eq. 3.18) - shape (6+num_j x 6+num_j)
    HH = np.block([[HH_bb,   HH_bq],
                   [HH_bq.T, HH_qq ]])

    return HH
