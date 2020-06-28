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
    Fe = np.zeros((3, num_b))          # Forces on base and link endpoints
    Te = np.zeros((3, num_b))          # Moments on base and link endpoints
    tau = np.zeros(num_j)             # Torques/forces on joints

    return Fe, Te, tau


def r_ne(RR, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te, SS, SE, j_type, cc, ce,
         mass, inertia, Ez, Gravity, BB):
    """
    Inverse Dynamics computation by the Recursive Newton-Euler method.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies

    # Calculation of body velocity vectors (base + links)
    vv, ww = kin.calc_vel(AA, v0, w0, q, qd, BB, j_type, cc, Ez)

    # Calculation of body acceleration vectors (base + links)
    vd, wd = kin.calc_acc(AA, ww, vd0, wd0, q, qd, qdd, BB, j_type, cc, Ez)

    # Calculation of inertial forces and moments on body centroids (base + links)
    # Includes also the gravity effect
    FF = np.zeros((3, num_b))
    TT = np.zeros((3, num_b))
    for i in range(0, num_b):
        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T
        FF[:, i] = mass[i] * (vd[:, i] - Gravity)                             # Eq. 3.30
        TT[:, i] = In_I @ wd[:, i] + cross(ww[:, i], (In_I @ ww[:, i]))       # Eq. 3.31

    # Calculation of forces and moments on joints (Eqs. 3.32 and 3.33)
    Fjnt = np.zeros((3, num_j))
    Tjnt = np.zeros((3, num_j))

    # Loop over the joints from the last (index 0 is the fake joint on the base)
    for i in range(num_j, 0, -1):

        idxi = i - 1             # Index link/joint <i> in Fjnt, Tjnt, j_type, q
        F_tmp = np.zeros(3)
        T_tmp = np.zeros(3)

        # Forces (one link may have more than one upper connection, thus the loop)
        for j in range(i+1, num_j+1):
            idxj = j - 1        # Index link/joint <j> in Fjnt
            F_tmp = F_tmp + float(SS[i, j]) * Fjnt[:, idxj]        # Add previous joint/joints
        Fjnt[:, idxi] = FF[:, i] + F_tmp + float(SE[i]) * Fe[:, i]  # Eq. 3.32

        # Moments (Eq. 3.33)

        # Add attached links contribution (may have more than one upper connection)
        A_I_i = AA[:, 3*i:3*(i+1)]
        for j in range(i+1, num_j+1):
            idxj = j - 1        # Index link/joint <j> in j_type, q, Fjnt, Tjnt
            d = cc[:, i, j] - cc[:, i, i] + float(j_type[idxj] == 'P') * Ez * q[idxj]
            T_tmp = T_tmp + float(SS[i, j]) * (cross((A_I_i @ d), Fjnt[:, idxj]) + Tjnt[:, idxj])

        # Add inertial terms contribution (Rotational joint)
        if (j_type[idxi] == 'R'):
            Tjnt[:, idxi] = TT[:, i] + T_tmp - cross((A_I_i @ cc[:, i, i]), FF[:, i])

        # Add inertial terms contribution (Prismatic joint)
        elif (j_type[idxi] == 'P'):
            Tjnt[:, idxi] = TT[:, i] + T_tmp - cross((A_I_i @ (Ez*q[idxi] - cc[:, i, i])), FF[:, i])

        # Add endpoint contribution
        d = ce[:, i] - cc[:, i, i] + float(j_type[idxi] == 'P') * Ez * q[idxi]
        Tjnt[:, idxi] = Tjnt[:, idxi] - float(SE[i]) * (cross((A_I_i @ d), Fe[:, i]) + Te[:, i])

    # Reaction forces on the base (Eqs. 3.38 and 3.39)
    F_tmp = np.zeros(3)
    T_tmp = np.zeros(3)

    # Forces/moments from the links connected to the base
    for i in range(1, num_j+1):

        idxi = i - 1             # Index link/joint <i> in Fjnt, Tjnt

        # Add if the link is connected
        if (SS[0, i] > 0):
            F_tmp = F_tmp + float(SS[0, i]) * Fjnt[:, idxi]        # Eq. 3.38
            T_tmp = T_tmp + float(SS[0, i]) * cross((AA[:, 0:3] @ cc[:, 0, i]), Fjnt[:, idxi]) + \
                Tjnt[:, idxi]     # Eq. 3.39

    # Add inertial and external terms
    FF0 = FF[:, 0] + F_tmp + float(SE[i]) * Fe[:, 0]
    TT0 = TT[:, 0] + T_tmp + float(SE[i]) * Te[:, 0]

    # Calculation of torque at each joint (Eq. 3.36 and 3.37)
    tau = np.zeros(num_j)
    for i in range(1, num_j+1):

        idxi = i - 1                 # Index link/joint <i> in Fjnt, Tjnt, j_type, tau
        A_I_i = AA[:, 3*i:3*(i+1)]

        if (j_type[idxi] == 'R'):
            tau[idxi] = Tjnt[:, idxi] @ (A_I_i @ Ez)       # Eq. 3.36

        # Prismatic joint
        elif (j_type[idxi] == 'P'):
            tau[idxi] = Fjnt[:, idxi] @ (A_I_i @ Ez)       # Eq. 3.37

    # Compose generalized forces (return [FF0, TT0] if single-body)
    Force = np.block([FF0, TT0, tau])

    return Force


# def f_dyn(R0, Q0, v0, w0, q, qd, Fe, Te, tau, cc, ce, SS, SE, BB, j_type, Qi, Ez, mass, \
#     inertia, Gravity):
#     """Forward dynamics computation, i.e. return the accelerations given the forces,
#        moments, and torques on the system"""


# !!!!! check for calls to calc_je when no links

#     num_j = len(q) - 1      # Number of joints/links
#     num_b = num_j + 1       # Number of bodies (base + links)

#     # Calculation of rotation matrices
#     AA = kin.calc_aa(Q0, q, BB, j_type, Qi)

#     # Calculation of position vectors
#     RR = kin.calc_pos(R0, AA, q, BB, j_type, cc, Ez)

#     # Calculation of inertia matrix
#     HH = kin.calc_hh(RR, AA, mass, inertia, BB, j_type, cc, Ez)

#     # Calculation of velocty dependent term, Force0. This is obtained by the RNE inverse
#     # dynamics computation with accerelations, forces/moments, and torques set to zero.
#     vd0 = np.zeros(3)
#     wd0 = np.zeros(3)
#     qdd = np.zeros(num_j+1)
#     Fe = np.zeros((3,num_b))
#     Te = np.zeros((3,num_b))
#     tau = np.zeros(num_j+1)
#     Force0 = r_ne(RR, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te, SS, SE, j_type, cc, ce, mass, \
#     inertia, Ez, Gravity, BB)


# % Force = forces on the generalized coordinate.
# % Force_ex = forces on the end points.
# Force = zeros(6+num_q,1);
# Force_ex = zeros(6+num_q,1);
# % F0, T0 are forces on the centroid of the 0-th body.
# Force(1:3) = F0;
# Force(4:6) = T0;
# % If Multi body system, tau is a joint torque.
# if ( num_q ~= 0 )
#    Force(7:num_q+6) = tau;
# end
# % Calculate external forces
# % If single body system, no external forces.
# if num_q == 0
#    % Note that the body 0 cannot have an endpoint.
#    Fx   = zeros(3,1);
#    Tx   = zeros(3,1);
#    taux = [];
# % Multi body system
# else
#    Fx    = zeros(3,1);
#    Tx    = zeros(3,1);
#    taux  = zeros(num_q,1);
#    E_3 = eye(3,3);
#    O_3 = zeros(3,3);
#    num_e = 1;
#    for i = 1 : num_q
#       if SE(i)==1
#          joints = j_num(num_e);
#          tmp = calc_je(RR, AA, q, joints);
#          JJ_tx_i = tmp(1:3,:);
#          JJ_rx_i = tmp(4:6,:);
#          num_e = num_e + 1;
#          A_I_i = AA(:,i*3-2:i*3);
#          Re0i = RR(:,i) - R0 + A_I_i*ce(:,i);
#          Me_i = [         E_3      O_3;
#                   tilde(Re0i)      E_3;
#                      JJ_tx_i'  JJ_rx_i' ];
#          F_ex(:,i) = Me_i * [ Fe(:,i) ; Te(:,i) ];
#       end
#    end
#    for i = 1 : num_q
#       Fx   = Fx   + F_ex(1:3,i);
#       Tx   = Tx   + F_ex(4:6,i);
#       taux = taux + F_ex(7:6+num_q,i);
#    end
# end
# Force_ex(1:3) = Fx;
# Force_ex(4:6) = Tx;
# Force_ex(7:6+num_q) = taux;
# % Calculation of the acclelation
# a_Force = Force - Force0 + Force_ex;
# Acc = HH\a_Force;
# %Acc = inv(HH)*a_Force;
# vd0 = Acc(1:3);
# wd0 = Acc(4:6);
# qdd = Acc(7:6+num_q);
# if num_q == 0
#    qdd=[];
# end
