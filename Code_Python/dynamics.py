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
from utils import cross, rotW, tilde


def r_ne(RR, AA, v0, w0, q, qd, vd0, wd0, qdd, Fe, Te, SS, SE, BB, j_type, cc,
         ce, mass, inertia):
    """
    Inverse dynamics computation by the recursive Newton-Euler method
    (eqs. 3.30-3.39).

    Given:
    state                               RR, AA, v0, w0, q, qd
    accelerations                       vd0, wd0, qdd
    external forces and moments         Fe, Te

    Returns:
    internal forces and moments         F0, T0, tau

    i.e. the system reaction forces and moments (base + joints).
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = len(SE)             # Number of endpoints
    Ez = np.array([0.0, 0.0, 1.0])              # Joint axis direction
    Gravity = np.array([0.0, 0.0, -9.81])       # Gravity vector

    # Linear and angular velocities of all bodies
    vv, ww = kin.calc_vel(AA, v0, w0, q, qd, BB, j_type, cc)

    # Linear and angular accelerations of all bodies
    vd, wd = kin.calc_acc(AA, ww, vd0, wd0, q, qd, qdd, BB, j_type, cc)

    # Inertial forces and moments on the body centroids, (for convenience the
    # the gravitational force is also included here) - eqs. 3.30-3.31
    F_in = np.zeros((3, num_b))
    T_in = np.zeros((3, num_b))
    for i in range(num_b):

        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I_i = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T

        # Eq. 3.30 and 3.31
        F_in[:, i] = mass[i] * (vd[:, i] - Gravity)
        T_in[:, i] = In_I_i @ wd[:, i] + cross(ww[:, i], (In_I_i @ ww[:, i]))

    # Forces and moments on the joints (eqs. 3.32-3.35)
    F_jnt = np.zeros((3, num_j))
    T_jnt = np.zeros((3, num_j))

    # Start from the last link
    for i in range(num_j, 0, -1):

        idxi = i - 1            # Index joint <i> in Fjnt, Tjnt, j_type, q
        A_I_i = AA[:, 3*i:3*(i+1)]
        is_P = float(j_type[idxi] == 'P')       # = 1 if joint is prismatic

        # Vector centroid <i> to joint <i>
        L_ii = cc[:, i, i] - is_P * Ez * q[idxi]

        # Add inertial force and moment on the link centroid (the cross-product
        # is negative because the arm should go from the joint to the centroid)
        F_jnt[:, idxi] = F_in[:, i]
        T_jnt[:, idxi] = T_in[:, i] - cross(A_I_i @ L_ii, F_in[:, i])

        # Add force and moment due to connected upper links (note that a link
        # may have more than one upper connection)
        for j in range(i+1, num_j+1):

            idxj = j - 1        # Index joint <j> in j_type, q, F_jnt, T_jnt

            # Add contribution if link <j> is an upper connection of link <i>
            if (SS[i, j]):

                # Vector joint <i> to joint <j>
                L_ij = cc[:, i, j] - cc[:, i, i] + is_P * Ez * q[idxi]

                # Add joint force and moment (note that Fj and Tj are calculated
                # as joint force and moment on the j-th link, thus the force and
                # moment passed to the i-th link are equal and opposite)
                F_jnt[:, idxi] += F_jnt[:, idxj]
                T_jnt[:, idxi] += cross(A_I_i @ L_ij, F_jnt[:, idxj]) \
                                  + T_jnt[:, idxj]

        # Add external forces and moments
        for ie in range(num_e):

            # Add contribution if endpoint <ie> is connected to link <i>
            if (SE[ie] == i):

                # Vector centroid <i> to endpoint <ie>
                L_ie = ce[:, ie] - cc[:, i, i] + is_P * Ez * q[idxi]

                # Add external force and moment
                F_jnt[:, idxi] -= Fe[:, ie]
                T_jnt[:, idxi] -= (cross(A_I_i @ L_ie, Fe[:, ie]) + Te[:, ie])

    # Reaction force and moment on the base centroid (Eqs. 3.38 and 3.39)
    F0 = np.zeros(3)
    T0 = np.zeros(3)

    # Add inertial force and moment of the base
    F0 += F_in[:, 0]
    T0 += T_in[:, 0]

    # Add forces and moments from the links connected to the base
    for i in range(1, num_j+1):

        # Add if link <i> is connected
        if (SS[0, i] > 0):
            idxi = i - 1             # Index link/joint <i> in Fjnt, Tjnt
            F0 += F_jnt[:, idxi]
            T0 += cross(AA[:, 0:3] @ cc[:, 0, i], F_jnt[:, idxi]) \
                  + T_jnt[:, idxi]

    # Add external forces and moments
    for ie in range(num_e):

        # Add contribution if endpoint <ie> is connected to the base
        if (SE[ie] == 0):
            F0 -= Fe[:, ie]
            T0 -= (cross(AA[:, 0:3] @ ce[:, ie], Fe[:, ie]) + Te[:, ie])

    # Calculation of joint torques/forces (eq. 3.36 and 3.37)
    tau = np.zeros(num_j)
    for i in range(1, num_j+1):

        idxi = i - 1        # Index link/joint <i> in Fjnt, Tjnt, j_type, tau
        Ez_I_i = AA[:, 3*i:3*(i+1)] @ Ez

        # If revolute joint (eq. 3.36)
        if (j_type[idxi] == 'R'):
            tau[idxi] = T_jnt[:, idxi] @ Ez_I_i

        # If prismatic joint (eq. 3.37)
        elif (j_type[idxi] == 'P'):
            tau[idxi] = F_jnt[:, idxi] @ Ez_I_i

    # Compose generalized forces
    Force = np.block([F0, T0, tau])

    return Force


def f_dyn(R0, A0, v0, w0, q, qd, F0, T0, tau, Fe, Te, SS, SE, BB, j_type, cc,
          ce, mass, inertia, Qi, Qe):
    """
    Forward dynamics computation (chapter 3.6).

    Given:
    state                               R0, A0, v0, w0, q, qd
    internal forces and moments         F0, T0, tau
    external forces and moments         Fe, Te

    Returns:
    accelerations                       vd0, wd0, qdd
    
    i.e. the system linear and angular accelerations (base + joints).
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = len(SE)             # Number of endpoints
    Ez = np.array([0.0, 0.0, 1.0])              # Joint axis direction
    Gravity = np.array([0.0, 0.0, -9.81])       # Gravity vector
    Force = np.zeros(6+num_j)

    # Rotation matrices
    AA = kin.calc_aa(A0, q, BB, j_type, Qi)

    # Position vectors
    RR = kin.calc_pos(R0, AA, q, BB, j_type, cc)

    # Inertia matrice
    HH = kin.calc_hh(RR, AA, mass, inertia, BB, j_type, cc)

    # Calculation of velocty dependent terms using the recursive Newton-Eulero
    # inverse dynamics setting to zero all the accelerations and all the
    # external forces.
    acc0 = np.zeros(3)
    qdd0 = np.zeros(num_j)
    Fe0 = np.zeros((3, num_e))
    Force0 = r_ne(RR, AA, v0, w0, q, qd, acc0, acc0, qdd0, Fe0, Fe0, SS, SE,
                  BB, j_type, cc, ce, mass, inertia)

    # Generalized control terms on the base centroid and on the joints
    Force = np.block([F0, T0, tau])

    # Generalized terms due to the force and moments applied to the endpoints
    Fx = np.zeros(3)
    Tx = np.zeros(3)
    taux = np.zeros(num_j)

    # Loop over all endpoints
    for ie in range(num_e):

        i = SE[ie]          # Link associated to the endpoint <ie>

        # If the endpoint is associated with a link
        if (i > 0):

            # Link sequence
            seq_link = kin.j_num(i, BB)

            # Endpoint Jacobian - shape is (6 x num_j)
            JJ_tmp = kin.calc_je(RR, AA, q, seq_link, j_type, cc, ce[:, ie],
                                 Qe[:, ie])
            JJ_tx_i = JJ_tmp[0:3, :]        # Translational component
            JJ_rx_i = JJ_tmp[3:6, :]        # Rotational component

            # Endpoint position wrt the base centroid
            A_I_i = AA[:, 3*i:3*(i+1)]
            R_0_ie = RR[:, i] - RR[:, 0] + A_I_i @ ce[:, ie]

            # Generalized terms due to this endpoint
            Fx += Fe[:, ie]
            Tx += tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie]
            taux += JJ_tx_i.T @ Fe[:, ie] + JJ_rx_i.T @ Te[:, ie]

        # If the endpoint is associated with the base
        else:

            # Generalized terms due to this endpoint (no contribution to the
            # joint terms)
            R_0_ie = AA[:, 0:3] @ ce[:, ie]
            Fx += Fe[:, ie]
            Tx += tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie]

    # Assemble the endpoint contributions
    Force_ee = np.block([Fx, Tx, taux])

    # Calculation of the acceleration - eq. 3.29
    acc = np.linalg.inv(HH) @ (Force + Force_ee - Force0)

    vd0 = acc[0:3]
    wd0 = acc[3:6]
    qdd = acc[6:]

    return vd0, wd0, qdd


def f_dyn_nb(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce, BB,
             j_type, Qi, cc, mass, inertia, Qe, SS):
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
    w0_pred = w0 + k4 * (wd0 + wd0_pred)
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
        w0_corr = w0 + k4 * (wd0 + wd0_corr)
        A0_corr = rotW(w0_corr, dt) @ A0

        # Next step
        R0_pred = R0_corr
        A0_pred = A0_corr
        v0_pred = v0_corr
        w0_pred = w0_corr
        q_pred = q_corr
        qd_pred = qd_corr

    return R0_pred, A0_pred, v0_pred, w0_pred, q_pred, qd_pred


def f_dyn_rk(R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau, dt, SE, ce, BB,
             j_type, Qi, cc, mass, inertia, Qe, SS):
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
    k1_q = qd * dt
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
    k2_q = (qd + k1_qd / 2.0) * dt
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
    k3_q = (qd + k2_qd / 2.0) * dt
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
    k4_q = (qd + k3_qd) * dt
    k4_v0 = vd0_tmp * dt
    k4_w0 = wd0_tmp * dt
    k4_qd = qdd_tmp * dt

    # Compute Values at the Next Time Step
    R0_pred = R0 + (k1_R0 + 2.0 * k2_R0 + 2.0 * k3_R0 + k4_R0) / 6.0
    A0_pred = A0 + (k1_A0 + 2.0 * k2_A0 + 2.0 * k3_A0 + k4_A0) / 6.0
    q_pred = q + (k1_q + 2.0 * k2_q + 2.0 * k3_q + k4_q) / 6.0
    v0_pred = v0 + (k1_v0 + 2.0 * k2_v0 + 2.0 * k3_v0 + k4_v0) / 6.0
    w0_pred = w0 + (k1_w0 + 2.0 * k2_w0 + 2.0 * k3_w0 + k4_w0) / 6.0
    qd_pred = qd + (k1_qd + 2.0 * k2_qd + 2.0 * k3_qd + k4_qd) / 6.0

    return R0_pred, A0_pred, v0_pred, w0_pred, q_pred, qd_pred
