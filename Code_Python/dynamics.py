"""
SpaceDyn - A Toolbox for Space and Mobile Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np

import kinematics as kin
from utils import cross


def calc_Forces(num_j=0):
    """
    Returns the external forces/moments/torques applied to base, links, and
    joints.
    """
    num_b = num_j + 1
    Fe = np.zeros((3, num_b))           # Forces on base and link end-points
    Te = np.zeros((3, num_b))           # Moments on base and link end-points
    tau = np.zeros(num_j)               # Torques/forces on joints

    return Fe, Te, tau


def r_ne(RR, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te, SS, SE, j_type, cc, ce,
         mass, inertia, Ez, Gravity, BB):
    """
    Inverse Dynamics computation by the Recursive Newton-Euler method.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies

    # Base and links velocities
    vv, ww = kin.calc_vel(AA, v0, w0, q, qd, BB, j_type, cc, Ez)

    # Base and links accelerations
    vd, wd = kin.calc_acc(AA, ww, vd0, wd0, q, qd, qdd, BB, j_type, cc, Ez)

    # Inertial forces and moments on the base and link centroids (includes
    # also the gravity effects)
    FF = np.zeros((3, num_b))
    TT = np.zeros((3, num_b))
    for i in range(0, num_b):

        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T

        # Eq. 3.30 and 3.31
        FF[:, i] = mass[i] * (vd[:, i] - Gravity)
        TT[:, i] = In_I @ wd[:, i] + cross(ww[:, i], (In_I @ ww[:, i]))

    # Forces and moments on joints (Eqs. 3.32 and 3.33)
    Fjnt = np.zeros((3, num_j))
    Tjnt = np.zeros((3, num_j))

    # Loop over the joints from the last
    for i in range(num_j, 0, -1):

        idxi = i - 1            # Index link/joint <i> in Fjnt, Tjnt, j_type, q
        F_tmp = np.zeros(3)
        T_tmp = np.zeros(3)

        # Forces - Eq. 3.32 (each link may have more than one upper connection)
        for j in range(i+1, num_j+1):
            idxj = j - 1        # Index link/joint <j> in Fjnt
            F_tmp = F_tmp + float(SS[i, j]) * Fjnt[:, idxj]
        Fjnt[:, idxi] = FF[:, i] + F_tmp + float(SE[i]) * Fe[:, i]

        # Moments (Eq. 3.33)

        # Add attached links contribution (each link may have more than one
        # upper connection)
        A_I_i = AA[:, 3*i:3*(i+1)]
        for j in range(i+1, num_j+1):
            idxj = j - 1        # Index link/joint <j> in j_type, q, Fjnt, Tjnt
            d = cc[:, i, j] - cc[:, i, i] + \
                float(j_type[idxj] == 'P') * Ez * q[idxj]
            T_tmp = T_tmp + float(SS[i, j]) * (cross((A_I_i @ d), Fjnt[:, idxj])
                                               + Tjnt[:, idxj])

        # Add inertial terms contribution (rotational joint)
        if (j_type[idxi] == 'R'):
            Tjnt[:, idxi] = TT[:, i] + T_tmp \
                                     - cross((A_I_i @ cc[:, i, i]), FF[:, i])

        # Add inertial terms contribution (prismatic joint)
        elif (j_type[idxi] == 'P'):
            Tjnt[:, idxi] = TT[:, i] + T_tmp \
                - cross((A_I_i @ (Ez * q[idxi] - cc[:, i, i])), FF[:, i])

        # Add endpoint contribution
        d = ce[:, i] - cc[:, i, i] + float(j_type[idxi] == 'P') * Ez * q[idxi]
        Tjnt[:, idxi] = Tjnt[:, idxi] - float(SE[i]) \
                        * (cross((A_I_i @ d), Fe[:, i]) + Te[:, i])

    # Reaction forces on the base (Eqs. 3.38 and 3.39)
    F_tmp = np.zeros(3)
    T_tmp = np.zeros(3)

    # Forces/moments from the links connected to the base
    for i in range(1, num_j+1):

        idxi = i - 1             # Index link/joint <i> in Fjnt, Tjnt

        # Add if the link is connected (Eq. 3.38 and 3.39)
        if (SS[0, i] > 0):
            F_tmp = F_tmp + float(SS[0, i]) * Fjnt[:, idxi]
            T_tmp = T_tmp + float(SS[0, i]) \
                    * cross((AA[:, 0:3] @ cc[:, 0, i]), Fjnt[:, idxi]) + \
                    Tjnt[:, idxi]

    # Add inertial and external terms
    FF0 = FF[:, 0] + F_tmp + float(SE[i]) * Fe[:, 0]
    TT0 = TT[:, 0] + T_tmp + float(SE[i]) * Te[:, 0]

    # Calculation of torque at each joint (Eq. 3.36 and 3.37)
    tau = np.zeros(num_j)
    for i in range(1, num_j+1):

        idxi = i - 1        # Index link/joint <i> in Fjnt, Tjnt, j_type, tau
        Ez_I_i = AA[:, 3*i:3*(i+1)] @ Ez

        if (j_type[idxi] == 'R'):
            tau[idxi] = Tjnt[:, idxi] @ Ez_I_i      # Eq. 3.36

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            tau[idxi] = Fjnt[:, idxi] @ Ez_I_i       # Eq. 3.37

    # Compose generalized forces (return [FF0, TT0] if single-body)
    Force = np.block([FF0, TT0, tau])

    return Force
