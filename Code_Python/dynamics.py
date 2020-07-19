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


def f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, SE, ce, BB,
          j_type, Qi, cc, Ez, mass, inertia, Qe, SS):
    """
    Forward dynamics computation: returns the accelerations given the state
    (positions/orientations and velocities) and any external input.
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = len(SE)

    # Rotation matrices
    AA = calc_aa(A0, q, BB, j_type, Qi)

    # Position vectors
    RR = calc_pos(R0, AA, q, BB, j_type, cc, Ez)

    # Inertia matrice
    HH = calc_hh(RR, AA, mass, inertia, BB, j_type, cc, Ez)

    # Calculation of velocty dependent term accomplished by recursive Newton
    # Eulero inverse dynamics with accelerations and external forces set to 0.
    vd0 = np.zeros(3)
    wd0 = np.zeros(3)
    qdd = np.zeros(num_j)
    Fe  = np.zeros((3, num_e))
    Te  = np.zeros((3, num_e))
    Force0 = r_ne(RR, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te, SS, SE,
                  j_type, cc, ce, mass, inertia, Ez, Gravity, BB)

    # Forces on the generalized coordinate
    Force = np.zeros(6+num_j)
    Force[0:3] = F0
    Force[3:6] = T0
    if (num_j > 0):
        Force[6:] = tau

    # Forces on the endpoints
    Fx = np.zeros(3)
    Tx = np.zeros(3)
    taux = np.zeros(num_j)

    # Loop over all endpoints !!!! check for base SE = 0 !!!!
    for ie in range(num_e):

        i = SE[ie]                      # Link associated to endpoint

        # Endpoint associated with a link
        if (i > 0):

            # Link sequence to endpoint
            seq_link = j_num(i, BB)

            # Jacobian associated to this endpoint - shape is (6 x num_j)
            JJ_tmp = calc_je(RR, AA, q, seq_link, j_type, cc, ce[:, ie],
                             Qe[:, ie], Ez)
            JJ_tx_i = JJ_tmp[0:3, :]        # Translational component
            JJ_rx_i = JJ_tmp[3:6, :]        # Rotational component

            # Endpoint position wrt the base centroid
            A_I_i = AA[:, 3*i:3*(i+1)]
            R_0_ie = RR[:, i] - RR[:, 0] + A_I_i @ ce[:, ie]

            # Generalized forces associated with this endpoint
            Fx += Fe[:, ie]
            Tx += tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie]
            taux += JJ_tx_i.T @ Fe[:, ie] + JJ_rx_i.T @ Te[:, ie]

        # Endpoint associated with the base
        else:

            # Generalized forces associated with this endpoint
            R_0_ie = AA[:, 0:3] @ ce[:, ie]
            Fx += Fe[:, ie]
            Tx += tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie]

    # Copy values
    Force_ee = zeros(num_e)
    Force_ee[0:3] = Fx
    Force_ee[3:6] = Tx
    Force_ee[6:] = taux

    # Calculation of the acceleration - eq. 3.29 ( !!!! check signs !!!!!)
    acc = np.linalg.inv(HH) * (Force - Force0 + Force_ee)

    vd0 = acc[0:3]
    wd0 = acc[3:6]
    qdd = acc[6:]

    return vd0, wd0, qdd


def f_dyn_nb2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce, BB,
              j_type, Qi, cc, Ez, mass, inertia, Qe, SS):
    """
    Integration of the equations of motion using the Newmark-beta method and
    the Rodrigues formula to update the rotations.
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
                                      tau, SE, ce, BB, j_type, Qi, cc, Ez,
                                      mass, inertia, Qe, SS)

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
    A0_pred = rotW(w0_pred, dt) @ A0

    # 2nd step: correction
    for i in range(n_reps):

        # Get the acceleration terms for current state
        vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0_pred, A0_pred, v0_pred, w0_pred,
                                          q_pred, qd_pred, F0, T0, Fe, Te,
                                          tau, SE, ce, BB, j_type, Qi, cc,
                                          Ez, mass, inertia, Qe, SS)

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
        A0_corr = rotW(w0_corr, dt) @ A0
   
        # Next step
        R0_pred = R0_corr
        A0_pred = A0_corr
        v0_pred = v0_corr
        w0_pred = w0_corr
        q_pred  = q_corr
        qd_pred = qd_corr

    return R0_pred, A0_pred, v0_pred, w0_pred, q_pred, qd_pred 


def f_dyn_rk2(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce, BB,
              j_type, Qi, cc, Ez, mass, inertia, Qe, SS):
    """
    Integration of the equations of motion using the Runge-Kutta method and
    the Rodrigues formula to update the rotations.
    """
    # 1st step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te,
                                      tau, SE, ce, BB, j_type, Qi, cc, Ez,
                                      mass, inertia, Qe, SS)

    k1_R0 = v0 * dt
    k1_A0 = rotW(w0, dt) @ A0 - A0
    k1_q  = qd * dt
    k1_v0 = vd0_tmp * dt
    k1_w0 = wd0_tmp * dt
    k1_qd = qdd_tmp * dt

    # 2nd Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k1_R0 / 2.0, A0 + k1_A0 / 2.0,
                                      v0 + k1_v0 / 2.0, w0 + k1_w0 / 2.0,
                                      q + k1_q / 2.0, qd + k1_qd / 2.0,
                                      F0, T0, Fe, Te, tau, SE, ce, BB,
                                      j_type, Qi, cc, Ez, mass, inertia, Qe, SS)

    k2_R0 = (v0 + k1_v0 / 2.0) * dt
    k2_A0 = rotW(w0 + k1_w0 / 2.0, dt) @ A0 - A0
    k2_q  = (qd + k1_qd / 2.0) * dt
    k2_v0 = vd0_tmp * dt
    k2_w0 = wd0_tmp * dt
    k2_qd = qdd_tmp * dt

    # 3rd Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k2_R0 / 2.0, A0 + k2_A0 / 2.0,
                                      v0 + k2_v0 / 2.0, w0 + k2_w0 / 2.0,
                                      q + k2_q / 2.0, qd + k2_qd / 2.0,
                                      F0, T0, Fe, Te, tau, SE, ce, BB,
                                      j_type, Qi, cc, Ez, mass, inertia, Qe, SS)

    k3_R0 = (v0 + k2_v0 / 2.0) * dt
    k3_A0 = rotW(w0 + k2_w0 / 2.0, dt) @ A0 - A0
    k3_q  = (qd + k2_qd / 2.0) * dt
    k3_v0 = vd0_tmp * dt
    k3_w0 = wd0_tmp * dt
    k3_qd = qdd_tmp * dt

    # 4th Step
    vd0_tmp, wd0_tmp, qdd_tmp = f_dyn(R0 + k3_R0, A0 + k3_A0, v0 + k3_v0,
                                      w0 + k3_w0, q + k3_q, qd + k3_qd,
                                      F0, T0, Fe, Te, tau, SE, ce, BB,
                                      j_type, Qi, cc, Ez, mass, inertia, Qe, SS)

    k4_R0 = (v0 + k3_v0) * dt
    k4_A0 = rotW(w0 + k3_w0, dt) @ A0 - A0
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


