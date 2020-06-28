"""
SpaceDyn - A Toolbox for Space and Mobile Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np

from utils import rpy2dc, cross, tilde


def j_num(num_link, BB):
    """
    Returns the link sequence from the base (excluded) to the endpoint
    associated with the link defined by number <num_link>.
    """
    seq_link = [num_link]           # Add current joint

    # Walk back until the base (j_number = 0) is reached
    j_number = BB[num_link-1]          # Previous link number
    while (j_number != 0):
        seq_link.insert(0, j_number)
        j_number = BB[j_number-1]      # Previous link number

    return seq_link


def f_kin_e(RR, AA, seq_link, Qe, ce):
    """
    Returns position and orientation of the endpoint defined by sequence
     <seq_link> (the link number of this endpoint is the last value in
     <seq_link>).
     """
    # Endpoint link number
    i = seq_link[-1]

    # Orientation
    A_I_i = AA[:, 3*i:3*(i+1)]           # Rotation matrix link wrt the inertial frame
    A_I_e = rpy2dc(Qe[:, i]).T         # Rotation matrix from joint/centroid to endpoint
    ORI_e = A_I_i @ A_I_e               # Rotation matrix wrt the inertial frame

    # Position
    POS_e = RR[:, i] + A_I_i @ ce[:, i]

    return POS_e, ORI_e


def f_kin_j(RR, AA, q, seq_link, j_type, cc, Ez):
    """
    Returns position and orientation of all joints in sequence <seq_link>
    (the link number of this endpoint is the last value in <seq_link>).
    """
    n_links = len(seq_link)
    POS_jnt = np.zeros((3, n_links))
    ORI_jnt = np.zeros((3, 3*n_links))

    # Loop over the joint sequence
    for i in range(n_links):

        j = seq_link[i]         # Link/joint index 
        idxj = j - 1             # Index link/joint <j> in q, j_type

        # Orientation
        A_I_j = AA[:, 3*j:3*(j+1)]
        ORI_jnt[:, 3*i:3*(i+1)] = A_I_j

        # Position (rotational joint)
        if (j_type[idxj] == 'R'):
            POS_jnt[:, i] = RR[:, j] + A_I_j @ cc[:, j, j]

        # Position (prismatic joint)
        elif (j_type[idxj] == 'P'):
            POS_jnt[:, i] = RR[:, j] + A_I_j @ (cc[:, j, j] - Ez*q[idxj])

    return POS_jnt, ORI_jnt


def calc_jte(RR, AA, q, seq_link, j_type, cc, ce, Qe, Ez):
    """
    Returns the translational Jacobian (3 x n_links) of the endpoint defined
    by sequence <seq_link> (the link number of this endpoint is the last value
    in <seq_link>).
    """
    n_links = len(seq_link)
    JJ_te = np.zeros((3, n_links))

    # Calculation of position and orientation for all joints in the sequence
    POS_jnt, ORI_jnt = f_kin_j(RR, AA, q, seq_link, j_type, cc, Ez)

    # Calculation of position and orientation for the endpoint
    POS_e, ORI_e = f_kin_e(RR, AA, seq_link, Qe, ce)

    # Loop over the sequence
    for i in range(n_links):

        j = seq_link[i]             # Link number in the sequence
        A_I_j = AA[:, 3*j:3*(j+1)]

        # Rotational joint
        if (j_type[i] == 'R'):
            JJ_te[:, i] = cross((A_I_j @ Ez), (POS_e - POS_jnt[:, i]))

        # Prismatic joint
        elif (j_type[i] == 'P'):
            JJ_te[:, i] = A_I_j @ Ez

    return JJ_te


def calc_jre(AA, seq_link, j_type, Ez):
    """
    Returns the rotational Jacobian (3 x n_links) of the endpoint defined by
    sequence <seq_link> (the link number of this endpoint is the last value
    in <seq_link>).
    """
    n_links = len(seq_link)
    JJ_re = np.zeros((3, n_links))

    # Loop over the sequence
    for i in range(n_links):

        j = seq_link[i]             # Link number in the sequence
        A_I_j = AA[:, 3*j:3*(j+1)]

        # Rotational joint
        if (j_type[i] == 'R'):
            JJ_re[:, i] = A_I_j @ Ez

        # Prismatic joint
        elif (j_type[i] == 'P'):
            pass                    # Already set to zero

    return JJ_re


def calc_je(RR, AA, q, seq_link, j_type, cc, ce, Qe, Ez):
    """
    Calculation of the Jacobian (6 x num_j) of the endpoint defined by
    sequence <seq_link>.
    """
    n_links = len(seq_link)
    num_j = len(q)
    Jacobian = np.zeros((6, num_j))

    # Calculation of Jacobian
    JJ_te = calc_jte(RR, AA, q, seq_link, j_type, cc, ce, Qe, Ez)   # Translation
    JJ_re = calc_jre(AA, seq_link, j_type, Ez)                      # Rotation
    JJ = np.block([[ JJ_te ],
                   [ JJ_re ]])

    # Assemble the endpoint Jacobian in the correct indeces (corresponding to the joints)
    for i in range(n_links):
        j = seq_link[i]             # Joint number in the sequence
        Jacobian[:, j-1] = JJ[:, i]

    return Jacobian


def calc_aa(Q0, q, BB, j_type, Qi):
    """
    Calculates the rotation matrices of all bodies (link centroid frame is set
    equal to the corresponding joint frame).
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    AA = np.zeros((3, 3*num_b))

    # Base
    AA[:, 0:3] = rpy2dc(Q0).T

    # Links
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, Qi
        k = BB[idxi]             # Connected previous link/joint

        # Rotational joint
        if (j_type[idxi] == 'R'):
            A_rel = rpy2dc(Qi[0, idxi], Qi[1, idxi], Qi[2, idxi] + q[idxi]).T

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            A_rel = rpy2dc(Qi[:, idxi]).T

        # Update rotation matrix
        AA[:, 3*i:3*(i+1)] = AA[:, 3*k:3*(k+1)] @ A_rel

    return AA


def calc_pos(R0, AA, q, BB, j_type, cc, Ez):
    """
    Calculates the position vector of all body centroids.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    RR = np.zeros((3, num_b))

    # Base
    RR[:, 0] = R0

    # Links
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, q
        k = BB[idxi]             # Connected previous link/joint

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Rotational joint
        if (j_type[idxi] == 'R'):
            RR[:, i] = RR[:, k] + A_I_k @ cc[:, k, i] - A_I_i @ cc[:, i, i]

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            RR[:, i] = RR[:, k] + A_I_k @ cc[:, k, i] - A_I_i @ (Ez*q[idxi] - cc[:, i, i])

    return RR


def calc_vel(AA, v0, w0, q, qd, BB, j_type, cc, Ez):
    """
    Calculates the linear velocity vector of all body centroids and the angular
    velocity vector of all bodies.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    vv = np.zeros((3, num_b))
    ww = np.zeros((3, num_b))

    # Base
    vv[:, 0] = v0
    ww[:, 0] = w0

    # Links
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, q, qd
        k = BB[idxi]             # Connected previous link/joint

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Rotational joint
        if (j_type[idxi] == 'R'):
            ww[:, i] = ww[:, k] + A_I_i @ (Ez*qd[idxi])
            vv[:, i] = vv[:, k] + cross(ww[:, k], (A_I_k @ cc[:, k, i])) \
                                - cross(ww[:, i], (A_I_i @ cc[:, i, i]))

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            ww[:, i] = ww[:, k]
            vv[:, i] = vv[:, k] + cross(ww[:, k], (A_I_k @ cc[:, k, i])) \
                                - cross(ww[:, i], (A_I_i @ cc[:, i, i])) \
                                + A_I_i @ (Ez*qd[idxi]) \
                                + cross(ww[:, i], (A_I_i @ (Ez*q[idxi])))

    return vv, ww


def calc_acc(AA, ww, vd0, wd0, q, qd, qdd, BB, j_type, cc, Ez):
    """
    Calculates the linear acceleration vector of all body centroids and the
    angular acceleration vector of all bodies.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    vd = np.zeros((3, num_b))
    wd = np.zeros((3, num_b))

    # Base
    vd[:, 0] = vd0
    wd[:, 0] = wd0

    # Links
    for i in range(1, num_b):

        idxi = i - 1             # Index link/joint <i> in BB, j_type, q, qd, qdd
        k = BB[idxi]             # Connected previous link/joint

        # Rotation matrices
        A_I_i = AA[:, 3*i:3*(i+1)]
        A_I_k = AA[:, 3*k:3*(k+1)]

        # Relative positions
        cc_I_i = A_I_i @ cc[:, i, i]
        cc_I_k = A_I_k @ cc[:, k, i]

        # Joint values
        q_I_i = A_I_i @ (Ez*q[idxi])
        qd_I_i = A_I_i @ (Ez*qd[idxi])
        qdd_I_i = A_I_i @ (Ez*qdd[idxi])

        # Rotational joint
        if (j_type[idxi] == 'R'):
            wd[:, i] = wd[:, k] + cross(ww[:, i], qd_I_i) + qdd_I_i
            vd[:, i] = vd[:, k] + cross(wd[:, k], cc_I_k) \
                                + cross(ww[:, k], cross(ww[:, k], cc_I_k)) \
                                - cross(ww[:, i], cross(ww[:, i], cc_I_i)) \
                                - cross(wd[:, i], cc_I_i)

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            wd[:, i] = wd[:, k]
            vd[:, i] = vd[:, k] + cross(wd[:, k], cc_I_k) \
                                + cross(ww[:, k], cross(ww[:, k], cc_I_k)) \
                                + cross(wd[:, i], q_I_i) \
                                + cross(ww[:, i], cross(ww[:, i], q_I_i)) \
                                + 2.0 * cross(ww[:, i], qd_I_i) + qdd_I_i \
                                - cross(ww[:, i], cross(ww[:, i], cc_I_i)) \
                                - cross(wd[:, i], cc_I_i)

    return vd, wd


def calc_jt(RR, AA, BB, j_type, cc, Ez):
    """
    Calculates the translational jacobians wrt the link centroids (base not
    included).
    """
    num_j = len(j_type)                     # Number of joints/links
    JJ_t = np.zeros((3, num_j*num_j))

    # Loop over all links
    for i in range(1, num_j+1):

        j = i               # Initial index from link/joint i to base

        # Follow the branch until the base is reached
        while(j > 0):

            idxj = j - 1         # Index link/joint <j> in BB, j_type
            A_I_j = AA[:, 3*j:3*(j+1)]

            # Rotational joint
            if (j_type[idxj] == 'R'):
                JJ_t[:, (i-1)*num_j+(j-1)] = \
                    cross((A_I_j @ Ez), (RR[:, i] - RR[:, j] - A_I_j @ cc[:, j, j]))

            # Prismatic joint
            elif (j_type[idxj] == 'P'):
                JJ_t[:, (i-1)*num_j+(j-1)] = A_I_j @ Ez

            j = BB[idxj]     # Previous link/joint along the branch

    return JJ_t


def calc_jr(AA, BB, j_type, Ez):
    """
    Calculates the rotational jacobians wrt the link centroids (base not
    included).
    """
    num_j = len(j_type)                     # Number of joints/links
    JJ_r = np.zeros((3, num_j*num_j))

    # Loop over all links
    for i in range(1, num_j+1):

        idxi = i - 1         # Index link/joint <i> in BB, j_type
        A_I_i = AA[:, 3*i:3*(i+1)]

        # Rotational joint
        if (j_type[idxi] == 'R'):
            JJ_r[:, (i-1)*num_j+(i-1)] = A_I_i @ Ez

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            JJ_r[:, (i-1)*num_j+(i-1)] = np.zeros((3, 1))

        j = BB[idxi]         # Previous link/joint along the branch
        if (j > 0):
            JJ_r[:, (i-1)*num_j:(i-1)*num_j+i-1] = JJ_r[:, (j-1)*num_j:(j-1)*num_j+i-1]

    return JJ_r


def calc_hh(RR, AA, mass, inertia, BB, j_type, cc, Ez):
    """
    Calculates the inertia matrices HH (6 + num_j) x (6 + num_j).
    """
    num_j = len(j_type)         # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    JJ_tg = np.zeros((3, num_j))
    HH_w = np.zeros((3, 3))
    HH_wq = np.zeros((3, num_j))
    HH_q = np.zeros((num_j, num_j))

    # Calculation of partial translational & rotational jacobians
    JJ_t = calc_jt(RR, AA, BB, j_type, cc, Ez)
    JJ_r = calc_jr(AA, BB, j_type, Ez)

    # Position of the gravity centroid (relative to the base)
    Rm = np.zeros(3)
    for i in range(num_b):
        Rm = Rm + mass[i] * RR[:, i]
    wr0g = Rm - RR[:, 0] * mass.sum()
    wr0g_tilde = tilde(wr0g)

    # Loop over all links
    for i in range(1, num_j+1):

        id1 = (i - 1) * num_j       # Start index joint/link i in matrices JJ_t and JJ_r
        id2 = i*num_j           # End index joint/link i in matrices JJ_t and JJ_r

        r0i = RR[:, i] - RR[:, 0]     # Position link centroid i wrt the base
        r0i_tilde = tilde(r0i)      # Matrix form of previous vector

        A_I_i = AA[:, 3*i:3*(i+1)]                               # Rotation matrix link i
        In_I_i = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T       # Inertia matrix link i

        # (3 x num_j)
        JJ_tg = JJ_tg + mass[i] * JJ_t[:, id1:id2]

        # (3 x 3)
        HH_w = HH_w + In_I_i + mass[i] * r0i_tilde.T @ r0i_tilde

        # (3 x num_j)
        HH_wq = HH_wq + In_I_i @ JJ_r[:, id1:id2] + mass[i] * r0i_tilde @ JJ_t[:, id1:id2]

        # (num_j x num_j)
        HH_q = HH_q + JJ_r[:, id1:id2].T @ In_I_i.T @ JJ_r[:, id1:id2] \
               + mass[i] * JJ_t[:, id1:id2].T @ JJ_t[:, id1:id2]

    # Base contribution to linear and rotational unknowns
    wE = mass.sum() * np.eye(3)
    HH_w = HH_w + AA[:, 0:3] @ inertia[:, 0:3] @ AA[:, 0:3].T

    # Assemble matrix HH
    HH = np.block([[         wE,  wr0g_tilde.T,  JJ_tg ],
                    [ wr0g_tilde,          HH_w,  HH_wq ],
                    [    JJ_tg.T,       HH_wq.T,  HH_q  ]])

    return HH
