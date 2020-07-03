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
from utils import cross, rotW


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
    Inverse Dynamics computation by the recursive Newton-Euler method.
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


def f_dyn_nb2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce):
    """
    Forward dynamics using the Newmark-beta method and the Rodrigues formula
    to update the rotations.
    """
    # Newmark parameters
    n_reps = 1
    beta = 1.0/6.0
    k1 = dt
    k2 = dt * dt / 3.0
    k3 = dt * dt * beta
    k4 = dt / 2.0

    # 1st step: prediction

    # Get the acceleration terms for current state
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te,
                                      tau, SE, ce)

    # Predict R0 and v0
    Rdd0 = vd0_tmp
    Rdd0_pred = Rdd0
    R0_pred = R0 + k1 * v0 + k2 * Rdd0 + k3 * Rdd0_pred
    Rd0_pred = v0 + k4 * (Rdd0 + Rdd0_pred)
    v0_pred = Rd0_pred

    # Predict q and qd
    qdd = qdd_tmp
    qdd_pred = qdd
    q_pred = q + k1 * qd + k2 * qdd + k3 * qdd_pred
    qd_pred = qd + k4 * (qdd + qdd_pred)

    # Predict w0 and A0
    wd0 = wd0_tmp
    wd0_pred = wd0
    w0_pred = w0 + k4 *(wd0 + wd0_pred)
    A0_pred = rotW(w0_pred) @ A0

    # 2nd step: correction
    for i in range(n_reps):

        # Get the acceleration terms for current state
        vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0_pred, A0_pred, v0_pred, w0_pred,
                                          q_pred, qd_pred, F0, T0, Fe, Te,
                                          tau, SE, ce)

        # Correct R0 and v0
        Rdd0 = vd0_tmp
        Rdd0_corr = Rdd0
        R0_corr = R0 + k1 * v0 + k2 * Rdd0 + k3 * Rdd0_corr
        Rd0_corr = v0 + k4 * (Rdd0 + Rdd0_corr)
        v0_corr = Rd0_corr
   
        # Correct q and qd
        qdd = qdd_tmp
        qdd_corr = qdd
        q_corr = q + k1 * qd + k2 * qdd + k3 * qdd_corr
        qd_corr = qd + k4 * (qdd + qdd_corr)
   
        # Correct w0 and A0
        wd0 = wd0_tmp
        wd0_corr = wd0
        w0_corr = w0 + k4 *(wd0 + wd0_corr)
        A0_corr = rotW(w0_corr) @ A0
   
        # Next step
        R0_pred = R0_corr
        A0_pred = A0_corr
        v0_pred = v0_corr
        w0_pred = w0_corr
        q_pred  = q_corr
        qd_pred = qd_corr

    return R0_pred, A0_pred, v0_pred, w0_pred, q_pred, qd_pred 


def f_dyn_rk2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce):
    """
    Forward dynamics using the Runge-Kutta method and the Rodrigues formula
    to update the rotations.
    """
    # 1st step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te,
                                      tau, SE, ce)

    k1_R0 = v0 * dt
    k1_A0 = rotW(w0) @ A0 - A0
    k1_q  = qd * dt
    k1_v0 = vd0_tmp * dt
    k1_w0 = wd0_tmp * dt
    k1_qd = qdd_tmp * dt

    # 2nd Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k1_R0 / 2.0, A0 + k1_A0 / 2.0,
                                      v0 + k1_v0 / 2.0, w0 + k1_w0 / 2.0,
                                      q + k1_q / 2.0, qd + k1_qd / 2.0,
                                      F0, T0, Fe, Te, tau, SE, ce)

    k2_R0 = (v0 + k1_v0 / 2.0) * dt
    k2_A0 = rotW(w0 + k1_w0 / 2.0) @ A0 - A0
    k2_q  = (qd + k1_qd / 2.0) * dt
    k2_v0 = vd0_tmp * dt
    k2_w0 = wd0_tmp * dt
    k2_qd = qdd_tmp * dt

    # 3rd Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k2_R0 / 2.0, A0 + k2_A0 / 2.0,
                                      v0 + k2_v0 / 2.0, w0 + k2_w0 / 2.0,
                                      q + k2_q / 2.0, qd + k2_qd / 2.0,
                                      F0, T0, Fe, Te, tau, SE, ce)

    k3_R0 = (v0 + k2_v0 / 2.0) * dt
    k3_A0 = rotW(w0 + k2_w0 / 2.0) @ A0 - A0
    k3_q  = (qd + k2_qd / 2.0) * dt
    k3_v0 = vd0_tmp * dt
    k3_w0 = wd0_tmp * dt
    k3_qd = qdd_tmp * dt

    # 4th Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k3_R0, A0 + k3_A0, v0 + k3_v0,
                                      w0 + k3_w0, q + k3_q, qd + k3_qd,
                                      F0, T0, Fe, Te, tau, SE, ce)

    k4_R0 = (v0 + k3_v0) * dt
    k4_A0 = rotW(w0 + k3_w0) @ A0 - A0
    k4_q  = (qd + k3_qd) * dt
    k4_v0 = vd0_tmp * dt
    k4_w0 = wd0_tmp * dt
    k4_qd = qdd_tmp * dt

    # Compute Values at the Next Time Step
    R0_pred = R0 + (k1_R0 + 2.0 * k2_R0 + 2.0 * k3_R0 + k4_R0) / 6.0
    A0_pred = A0 + (k1_A0 + 2.0 * k2_A0 + 2.0 * k3_A0 + k4_A0) / 6.0
    q_pred  = q  + (k1_q  + 2.0 * k2_q  + 2.0 * k3_q  + k4_q) / 6.0
    v0_pred = v0 + (k1_v0 + 2.0 * k2_v0 + 2.0 * k3_v0 + k4_v0) / 6.0
    w0_pred = w0 + (k1_w0 + 2.0 * k2_w0 + 2.0 * k3_w0 + k4_w0) / 6.0
    qd_pred = qd + (k1_qd + 2.0 * k2_qd + 2.0 * k3_qd + k4_qd) / 6.0

    return R0_pred, A0_pred, v0_pred, w0_pred, q_pred, qd_pred 


def f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, SE, ce):
    """
    Forward dynamics: returns the accelerations given the state and any
    external input.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies

    # Rotation matrices
    AA = calc_aa(A0, q)

    # Position vectors
    RR = calc_pos(R0, A0, AA, q)

    # Inertia matrice
    HH = calc_hh(R0, RR, A0, AA)

    # # Calculation of velocty dependent term, accomplished by recursive Newton
    # # Eulero inverse dynamics with accelerations and external forces set to 0.
    # qdd0 = np.zeros((num_j, 1))
    # acc0 = np.zeros((3, 1))
    # fe0  = np.zeros((3, num_j))
    # Force0 = r_ne(R0, RR, A0, AA, v0, w0, acc0, acc0, q, qd, qdd0, fe0, fe0 )

    # # % Force = forces on the generalized coordinate.
    # # % Force_ex = forces on the end points.
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

    return vd0, wd0, qdd
