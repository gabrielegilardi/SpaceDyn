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
from utils import cross, rotW, tilde, rpy2dc


def r_ne(RR, AA, v0, w0, q, qd, vd0, wd0, qdd, Fe, Te, SS, SE, BB, j_type, cc,
         ce, mass, inertia, Qe):
    """
    Inverse dynamics computation by the recursive Newton-Euler method
    (eqs. 3.30-3.39).

    Given:
    state                               RR, AA, v0, w0, q, qd
    accelerations                       vd0, wd0, qdd
    external forces and moments         Fe, Te

    Returns:
    reaction forces and moments         F0, T0, tau

    i.e. the system reaction force <F0> and moment <T0> on the base, and the 
    torque/forces <tau> on the joints.
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

    # Inertial forces and moments on the body centroids (for convenience the
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
                A_i_ie = rpy2dc(Qe[:, ie]).T    # From endpoint to link/joint
                A_I_ie = A_I_i @ A_i_ie         # From endpoint to inertial
                L_i_ie = A_i_ie.T @ (ce[:, ie] - cc[:, i, i]
                                     + is_P * Ez * q[idxi])

                # Add external force and moment
                F_jnt[:, idxi] -= A_I_ie @ Fe[:, ie]
                T_jnt[:, idxi] -= A_I_ie @ (tilde(L_i_ie) @ Fe[:, ie]
                                            + Te[:, ie])

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

            A_0_ie = rpy2dc(Qe[:, ie]).T        # From endpoint to base
            A_I_ie = AA[:, 0:3] @ A_0_ie        # From endpoint to inertial
            R_0_ie = A_0_ie.T @ ce[:, ie]       # ce in endpoint frame

            F0 -= A_I_ie @ Fe[:, ie]
            T0 -= A_I_ie @ (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

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


def f_dyn(R0, A0, v0, w0, q, qd, Fe, Te, tau, SS, SE, BB, j_type, cc, ce, mass,
          inertia, Qi, Qe):
    """
    Forward dynamics computation (chapter 3.6).

    Given:
    state                               R0, A0, v0, w0, q, qd
    external forces and moments         Fe, Te, tau

    Returns:
    accelerations                       vd0, wd0, qdd

    i.e. the system linear and angular accelerations (base + joints).
    """
    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = len(SE)             # Number of endpoints
    Ez = np.array([0.0, 0.0, 1.0])              # Joint axis direction
    Force = np.zeros(6+num_j)

    # Rotation matrices
    AA = kin.calc_aa(A0, q, BB, j_type, Qi)

    # Position vectors
    RR = kin.calc_pos(R0, AA, q, BB, j_type, cc)

    # Inertia matrice
    HH = kin.calc_hh(RR, AA, mass, inertia, BB, j_type, cc)

    # Calculation of velocty dependent terms using the recursive Newton-Eulero
    # inverse dynamics setting to zero all accelerations and forces
    zero_acc = np.zeros(3)
    zero_qdd = np.zeros(num_j)
    zero_Fe = np.zeros((3, num_e))
    Force0 = r_ne(RR, AA, v0, w0, q, qd, zero_acc, zero_acc, zero_qdd, zero_Fe,
                  zero_Fe, SS, SE, BB, j_type, cc, ce, mass, inertia, Qe)

    # Generalized external forces applied on base centroid and joints
    F0 = np.zeros(3)
    T0 = np.zeros(3)

    # Loop over all endpoints
    for ie in range(num_e):

        # If the endpoint is associated with the base
        if (SE[ie] == 0):

            A_0_ie = rpy2dc(Qe[:, ie]).T        # From endpoint to base
            A_I_ie = AA[:, 0:3] @ A_0_ie        # From endpoint to inertial
            R_0_ie = A_0_ie.T @ ce[:, ie]       # ce in endpoint frame

            F0 += A_I_ie @ Fe[:, ie]
            T0 += A_I_ie @ (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

    Force = np.block([F0, T0, tau])

    # Generalized external forces applied to the link endpoints
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
            A_I_ie = A_I_i @ rpy2dc(Qe[:, ie]).T
            R_0_ie = RR[:, i] - RR[:, 0] + A_I_i @ ce[:, ie]

            # Generalized terms due to this endpoint
            Fx += A_I_ie @ Fe[:, ie]
            Tx += tilde(R_0_ie) @ A_I_ie @ Fe[:, ie] + A_I_ie @ Te[:, ie]
            taux += JJ_tx_i.T @ A_I_ie @ Fe[:, ie] \
                    + JJ_rx_i.T @ A_I_ie @ Te[:, ie]

    # Assemble the link endpoint contributions
    Force_ee = np.block([Fx, Tx, taux])

    # Calculation of the acceleration - eq. 3.29
    acc = np.linalg.inv(HH) @ (Force + Force_ee - Force0)

    vd0 = acc[0:3]
    wd0 = acc[3:6]
    qdd = acc[6:]

    return vd0, wd0, qdd


def f_dyn_nb(dt, R0, A0, v0, w0, q, qd, Fe, Te, tau, SS, SE, BB, j_type, cc,
             ce, mass, inertia, Qi, Qe):
    """
    Integration of the equations of motion using the Newmark-beta method and
    the Rodrigues formula to update the rotations.

    Given:
    time-step                           dt
    state at time <t>                   R0, A0, v0, w0, q, qd, at <t>
    external forces and moments         Fe, Te, tau

    Returns:
    state at time <t+dt>                R0, A0, v0, w0, q, qd, at <t+dt>

    Notes:
    - equations are written assuming gamma = 1/2 to guarantee the 2nd order
      accuracy in the solution.
    - beta = 0 corresponds to the central difference method.
    - beta = 1/4 corresponds to the constant acceleration method.
    - beta = 1/6 corresponds to the linear acceleration method.
    - stability for any <dt> requires beta >= 1/4.
    """
    # Newmark parameters
    n_reps = 1
    beta = 1.0/6.0
    k1 = dt
    k2 = dt * dt / 3.0
    k3 = dt * dt * beta
    k4 = dt / 2.0

    # 1st step: prediction

    # Accelerations using the current state
    vd0, wd0, qdd = f_dyn(R0, A0, v0, w0, q, qd, Fe, Te, tau, SS, SE, BB,
                          j_type, cc, ce, mass, inertia, Qi, Qe)

    # Predicted values for vd0, v0, and R0
    vd0_p = vd0
    v0_p = v0 + k4 * (vd0 + vd0_p)
    R0_p = R0 + k1 * v0 + k2 * vd0 + k3 * vd0_p

    # Predict values for wd0, w0, and A0
    wd0_p = wd0
    w0_p = w0 + k4 * (wd0 + wd0_p)
    A0_p = rotW(w0_p, dt) @ A0

    # Predicted values for qdd, qd, and q
    qdd_p = qdd
    qd_p = qd + k4 * (qdd + qdd_p)
    q_p = q + k1 * qd + k2 * qdd + k3 * qdd_p

    # 2nd step: correction
    for i in range(n_reps):

        # Accelerations using the predicted state
        vd0, wd0, qdd = f_dyn(R0_p, A0_p, v0_p, w0_p, q_p, qd_p, Fe, Te, tau,
                              SS, SE, BB, j_type, cc, ce, mass, inertia, Qi, Qe)

        # Corrected values for vd0, v0, and R0
        vd0_c = vd0
        v0_c = v0 + k4 * (vd0 + vd0_c)
        R0_c = R0 + k1 * v0 + k2 * vd0 + k3 * vd0_c

        # Corrected values for wd0, w0, and A0
        wd0_c = wd0
        w0_c = w0 + k4 * (wd0 + wd0_c)
        A0_c = rotW(w0_c, dt) @ A0

        # Corrected values for qdd, qd, and q
        qdd_c = qdd
        qd_c = qd + k4 * (qdd + qdd_c)
        q_c = q + k1 * qd + k2 * qdd + k3 * qdd_c

        # Prepare next iteration
        R0_p = R0_c
        A0_p = A0_c
        v0_p = v0_c
        w0_p = w0_c
        q_p = q_c
        qd_p = qd_c

    return R0_c, A0_c, v0_c, w0_c, q_c, qd_c


def f_dyn_rk(dt, R0, A0, v0, w0, q, qd, Fe, Te, tau, SS, SE, BB, j_type, cc,
             ce, mass, inertia, Qi, Qe):
    """
    Integration of the equations of motion using the Runge-Kutta method and
    the Rodrigues formula to update the rotations.

    Given:
    time-step                           dt
    state at time <t>                   R0, A0, v0, w0, q, qd, at <t>
    external forces and moments         Fe, Te, tau

    Returns:
    state at time <t+dt>                R0, A0, v0, w0, q, qd, at <t+dt>
    """
    # 1st step

    # Accelerations using the current state
    vd0, wd0, qdd = f_dyn(R0, A0, v0, w0, q, qd, Fe, Te, tau, SS, SE, BB,
                          j_type, cc, ce, mass, inertia, Qi, Qe)

    # Predicted values for v0 and R0
    v0_s1 = vd0 * dt
    R0_s1 = v0 * dt

    # Predict values for w0 and A0
    w0_s1 = wd0 * dt
    A0_s1 = rotW(w0, dt) @ A0 - A0

    # Predicted values for qd and q
    qd_s1 = qdd * dt
    q_s1 = qd * dt

    # 2nd Step

    # Accelerations using the state predicted by the 1st step
    vd0, wd0, qdd = f_dyn(R0 + R0_s1 / 2.0, A0 + A0_s1 / 2.0, v0 + v0_s1 / 2.0,
                          w0 + w0_s1 / 2.0, q + q_s1 / 2.0, qd + qd_s1 / 2.0,
                          Fe, Te, tau, SS, SE, BB, j_type, cc, ce, mass,
                          inertia, Qi, Qe)

    # Predicted values for v0 and R0
    v0_s2 = vd0 * dt
    R0_s2 = (v0 + v0_s1 / 2.0) * dt

    # Predict values for w0 and A0
    w0_s2 = wd0 * dt
    A0_s2 = rotW(w0 + w0_s1 / 2.0, dt) @ A0 - A0

    # Predicted values for qd and q
    qd_s2 = qdd * dt
    q_s2 = (qd + qd_s1 / 2.0) * dt

    # 3rd Step

    # Accelerations using the state predicted by the 2nd step
    vd0, wd0, qdd = f_dyn(R0 + R0_s2 / 2.0, A0 + A0_s2 / 2.0, v0 + v0_s2 / 2.0,
                          w0 + w0_s2 / 2.0, q + q_s2 / 2.0, qd + qd_s2 / 2.0,
                          Fe, Te, tau, SS, SE, BB, j_type, cc, ce, mass,
                          inertia, Qi, Qe)

    # Predicted values for v0 and R0
    v0_s3 = vd0 * dt
    R0_s3 = (v0 + v0_s2 / 2.0) * dt

    # Predict values for w0 and A0
    w0_s3 = wd0 * dt
    A0_s3 = rotW(w0 + w0_s2 / 2.0, dt) @ A0 - A0

    # Predicted values for qd and q
    qd_s3 = qdd * dt
    q_s3 = (qd + qd_s2 / 2.0) * dt

    # 4th Step

    # Accelerations using the state predicted by the 3rd step
    vd0, wd0, qdd = f_dyn(R0 + R0_s3, A0 + A0_s3, v0 + v0_s3, w0 + w0_s3,
                          q + q_s3, qd + qd_s3, Fe, Te, tau, SS, SE, BB,
                          j_type, cc, ce, mass, inertia, Qi, Qe)

    # Predicted values for v0 and R0
    v0_s4 = vd0 * dt
    R0_s4 = (v0 + v0_s3) * dt

    # Predict values for w0 and A0
    w0_s4 = wd0 * dt
    A0_s4 = rotW(w0 + w0_s3, dt) @ A0 - A0

    # Predicted values for qd and q
    qd_s4 = qdd * dt
    q_s4 = (qd + qd_s3) * dt

    # Predicted values for the next step
    R0_p = R0 + (R0_s1 + 2.0 * R0_s2 + 2.0 * R0_s3 + R0_s4) / 6.0
    A0_p = A0 + (A0_s1 + 2.0 * A0_s2 + 2.0 * A0_s3 + A0_s4) / 6.0
    q_p = q + (q_s1 + 2.0 * q_s2 + 2.0 * q_s3 + q_s4) / 6.0
    v0_p = v0 + (v0_s1 + 2.0 * v0_s2 + 2.0 * v0_s3 + v0_s4) / 6.0
    w0_p = w0 + (w0_s1 + 2.0 * w0_s2 + 2.0 * w0_s3 + w0_s4) / 6.0
    qd_p = qd + (qd_s1 + 2.0 * qd_s2 + 2.0 * qd_s3 + qd_s4) / 6.0

    return R0_p, A0_p, v0_p, w0_p, q_p, qd_p
