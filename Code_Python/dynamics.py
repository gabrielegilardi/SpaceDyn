"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

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
    Inverse dynamics computation by the recursive Newton-Euler method
    (eqs. 3.30-3.39).

    FF, TT          Inertial forces and moments on the centroids
    Fj, Tj          Forces and moments on the joints
    Fe, Te          Forces and moments on the endpoints

    F0, T0          Reaction force and moment of the base
    tau             Torques/forces on the revolute/prismatic joints

    !!!! check signs !!!!!

    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = len(SE)             # Number of endpoints

    # Linear and angular velocities (all bodies)
    vv, ww = kin.calc_vel(AA, v0, w0, q, qd, BB, j_type, cc, Ez)

    # Linear and angular accelerations (all bodies)
    vd, wd = kin.calc_acc(AA, ww, vd0, wd0, q, qd, qdd, BB, j_type, cc, Ez)

    # Inertial forces and moments on the centroids, including the gravitational
    # force (eqs. 3.30-3.31) - base is included
    FF = np.zeros((3, num_b))
    TT = np.zeros((3, num_b))
    for i in range(num_b):

        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I_i = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T

        # Eq. 3.30 and 3.31
        FF[:, i] = mass[i] * (vd[:, i] - Gravity)
        TT[:, i] = In_I_i @ wd[:, i] + cross(ww[:, i], (In_I_i @ ww[:, i]))

    # Forces and moments on the joints (eqs. 3.32-3.35)
    Fj = np.zeros((3, num_j))
    Tj = np.zeros((3, num_j))

    # Start from the last link
    for i in range(num_j, 0, -1):

        idxi = i - 1            # Index joint <i> in Fjnt, Tjnt, j_type, q
        A_I_i = AA[:, 3*i:3*(i+1)]
        is_P = float(j_type[idxi] == 'P')       # Joint <i> is prismatic

        # Vector centroid <i> to joint <i>
        L_ii = cc[:, i, i] - is_P * Ez * q[idxi]

        # Add inertial force and moment (!!!! why moment is negative !!!!)
        Fj[:, idxi] = FF[:, i]
        Tj[:, idxi] = TT[:, i] - cross((A_I_i @ L_ii), FF[:, i])

        # Add connected links force and moments (may be more than one)
        for j in range(i+1, num_j+1):
            
            idxj = j - 1        # Index joint <j> in j_type, q, Fjnt, Tjnt

            # Add contribution if link <j> is connected to link <i>
            if (SS[i, j]):

                # Vector joint <i> to joint <j>
                L_ij = cc[:, i, j] - cc[:, i, i] + is_P * Ez * q[idxj]

                # Add joint force and moment
                Fj[:, idxi] += Fj[:, idxj]
                Tj[:, idxi] += cross(A_I_i @ L_ij, Fj[:, idxj]) + Tj[:, idxj]
        
        # Add external forces and moments
        for ie in range(num_e):

            # Add contribution if endpoint <ie> is connected to link <i>
            if (SE[ie] == i):

                # Vector centroid <i> to endpoint <ie>
                L_ie = ce[:, ie] - cc[:, i, i] + is_P * Ez * q[idxi]

                # Add external force and moment (!!!!! check sign both !!!!)
                Fj[:, idxi] += - Fe[:, ie]
                Tj[:, idxi] += - cross(A_I_i @ L_ie, Fe[:, ie]) - Te[:, ie]

    # Reaction force and moment on the base centroid (Eqs. 3.38 and 3.39)
    FF0 = np.zeros(3)
    TT0 = np.zeros(3)

    # Add inertial force and moment of the base
    FF0 += FF[:, 0]
    TT0 += TT[:, 0]
    
    # Add forces and moments from the links connected to the base
    for i in range(1, num_j+1):

        # Add if link <i> is connected
        if (SS[0, i] > 0):
            idxi = i - 1             # Index link/joint <i> in Fjnt, Tjnt
            FF0 += Fj[:, idxi]
            TT0 += cross(AA[:, 0:3] @ cc[:, 0, i], Fj[:, idxi]) + Tj[:, idxi]

    # Add external forces and moments !!! check signs !!!!!
    for ie in range(num_e):

        # Add contribution if endpoint <ie> is connected to the base
        if (SE[ie] == 0):
            FF0 += - Fe[:, ie]
            TT0 += - cross(AA[:, 0:3] @ ce[:, ie], Fe[:, ie]) - Te[:, ie]

    # Calculation of joint torques/forces (eq. 3.36 and 3.37)
    tau = np.zeros(num_j)
    for i in range(1, num_j+1):

        idxi = i - 1        # Index link/joint <i> in Fjnt, Tjnt, j_type, tau
        Ez_I_i = AA[:, 3*i:3*(i+1)] @ Ez

        # If revolute joint (eq. 3.36)
        if (j_type[idxi] == 'R'):
            tau[idxi] = Tj[:, idxi] @ Ez_I_i

        # If prismatic joint (eq. 3.37)
        elif (j_type[idxi] == 'P'):
            tau[idxi] = Fj[:, idxi] @ Ez_I_i

    # Compose generalized forces (return [FF0, TT0] if single-body)
    Force = np.block([FF0, TT0, tau])

    return Force


def f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, SE, ce):
    """
    Forward dynamics: returns the accelerations given the state and any
    external input.

    !!!! Correct ee calculations based on new SE format - calc_je !!!!
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
#          tmp = calc_je(RR, AA, q, joints);        !!!! new format for SE
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


