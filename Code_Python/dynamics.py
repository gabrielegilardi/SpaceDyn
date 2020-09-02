"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np
from kinematics import calc_aa, calc_pos, calc_vel, calc_acc, calc_hh, calc_je
from utils import cross, tilde, rpy2dc, dc2rpy

# Joint axis direction
Ez = np.array([0.0, 0.0, 1.0])


def r_ne(RR, AA, q, Yd, Ydd, Fe, Te, Conn, Prop):
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
    SS, SE, j_type = Conn[0], Conn[1], Conn[3]
    mass, inertia, cc, ce, Qe, gravity = \
        Prop[0], Prop[1], Prop[2], Prop[3], Prop[5], Prop[6]
    qd = Yd[6:]

    num_j = len(q)              # Number of joints/links
    num_b = num_j + 1           # Number of bodies
    num_e = SE.shape[1]         # Number of endpoints

    # Linear and angular velocities of all bodies
    vv, ww = calc_vel(AA, q, Yd, Conn, Prop)

    # Linear and angular accelerations of all bodies
    vd, wd = calc_acc(AA, ww, q, qd, Ydd, Conn, Prop)

    # Inertial forces and moments on the body centroids (for convenience the
    # the gravitational force is also included here) - eqs. 3.30-3.31
    F_in = np.zeros((3, num_b))
    T_in = np.zeros((3, num_b))
    for i in range(num_b):

        A_I_i = AA[:, 3*i:3*(i+1)]
        In_I_i = A_I_i @ inertia[:, 3*i:3*(i+1)] @ A_I_i.T

        # Eq. 3.30 and 3.31
        F_in[:, i] = mass[i] * (vd[:, i] - gravity)
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
            if (SE[0, ie] == i):

                # Vector centroid <i> to endpoint <ie>
                A_i_ie = rpy2dc(Qe[:, ie]).T        # Endpoint to link/joint
                A_I_ie = A_I_i @ A_i_ie             # Endpoint to inertial

                # If the external load is given wrt the local frame
                if (SE[1, ie] == 0):
                    L_i_ie = A_i_ie.T @ (ce[:, ie] - cc[:, i, i]
                                        + is_P * Ez * q[idxi])
                    F_jnt[:, idxi] -= A_I_ie @ Fe[:, ie]
                    T_jnt[:, idxi] -= A_I_ie @ (tilde(L_i_ie) @ Fe[:, ie]
                                                + Te[:, ie])

                # If the external load is given wrt the inertial frame
                else:
                    L_i_ie = A_I_ie @ (ce[:, ie] - cc[:, i, i]
                                        + is_P * Ez * q[idxi])
                    F_jnt[:, idxi] -= Fe[:, ie]
                    T_jnt[:, idxi] -= (tilde(L_i_ie) @ Fe[:, ie] + Te[:, ie])

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
        if (SE[0, ie] == 0):

            A_0_ie = rpy2dc(Qe[:, ie]).T        # Endpoint to base
            A_I_ie = AA[:, 0:3] @ A_0_ie        # Endpoint to inertial

            # If the external load is given wrt the local frame
            if (SE[1, ie] == 0):
                R_0_ie = A_0_ie.T @ ce[:, ie]
                F0 -= A_I_ie @ Fe[:, ie]
                T0 -= A_I_ie @ (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

            # If the external load is given wrt the inertial frame
            else:
                R_0_ie = A_0_ie @ ce[:, ie]
                F0 -= Fe[:, ie]
                T0 -= (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

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


def f_dyn(Y, Yd, Fe, Te, tau, Conn, Prop):
    """
    Forward dynamics computation (chapter 3.6).

    state at time <t>           Y = [R0, A0, q] & Yd = [v0, w0, qd] @ <t>
    external load               Fe, Te, tau

    Returns:
    accelerations               Ydd = [vd0, wd0, qdd] @ <t>

    i.e. the system linear and angular accelerations (base + joints).
    """
    SE = Conn[1]
    ce, Qe = Prop[3], Prop[5]
    R0, Q0, q = Y[0:3], Y[3:6], Y[6:]
    v0, w0, qd = Yd[0:3], Yd[3:6], Yd[6:]

    num_j = len(q)                              # Number of joints/links
    num_e = SE.shape[1]                         # Number of endpoints

    # Position and rotation matrices
    AA = calc_aa(Q0, q, Conn, Prop)
    RR = calc_pos(R0, AA, q, Conn, Prop)

    # Inertia matrice
    HH = calc_hh(RR, AA, Conn, Prop)

    # Calculation of velocity dependent terms using the recursive Newton-Eulero
    # inverse dynamics setting to zero all accelerations and forces
    zero_Ydd = np.zeros(6+num_j)
    zero_Fe = np.zeros((3, num_e))
    Force0 = r_ne(RR, AA, q, Yd, zero_Ydd, zero_Fe, zero_Fe, Conn, Prop)

    # Generalized external forces applied on base centroid and joints
    F0 = np.zeros(3)
    T0 = np.zeros(3)

    # Loop over all endpoints
    for ie in range(num_e):

        # If the endpoint is associated with the base
        if (SE[0, ie] == 0):

            A_0_ie = rpy2dc(Qe[:, ie]).T        # Endpoint to base
            A_I_ie = AA[:, 0:3] @ A_0_ie        # Endpoint to inertial

            # If the external load is given wrt the local frame
            if (SE[1, ie] == 0):
                R_0_ie = A_0_ie.T @ ce[:, ie]
                F0 += A_I_ie @ Fe[:, ie]
                T0 += A_I_ie @ (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

            # If the external load is given wrt the inertial frame
            else:
                R_0_ie = A_0_ie @ ce[:, ie]
                F0 += Fe[:, ie]
                T0 += (tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie])

    # Assemble all terms
    Force = np.block([F0, T0, tau])

    # Generalized external forces applied to the link endpoints
    Fx = np.zeros(3)
    Tx = np.zeros(3)
    taux = np.zeros(num_j)

    # Loop over all endpoints
    for ie in range(num_e):

        i = SE[0, ie]          # Link associated to the endpoint <ie>

        # If the endpoint is associated with a link
        if (i > 0):

            # Endpoint Jacobian - shape is (6 x num_j)
            JJ_tmp = calc_je(ie, RR, AA, q, Conn, Prop)
            JJ_tx_i = JJ_tmp[0:3, :]        # Translational component
            JJ_rx_i = JJ_tmp[3:6, :]        # Rotational component

            # Endpoint position wrt the base centroid
            A_I_i = AA[:, 3*i:3*(i+1)]
            A_I_ie = A_I_i @ rpy2dc(Qe[:, ie]).T
            R_0_ie = RR[:, i] - RR[:, 0] + A_I_i @ ce[:, ie]

            # If the external load is given wrt the local frame
            if (SE[1, ie] == 0):
                Fx += A_I_ie @ Fe[:, ie]
                Tx += tilde(R_0_ie) @ A_I_ie @ Fe[:, ie] + A_I_ie @ Te[:, ie]
                taux += + JJ_tx_i.T @ A_I_ie @ Fe[:, ie] \
                        + JJ_rx_i.T @ A_I_ie @ Te[:, ie]

            # If the external load is given wrt the inertial frame
            else:
                Fx += Fe[:, ie]
                Tx += tilde(R_0_ie) @ Fe[:, ie] + Te[:, ie]
                taux += JJ_tx_i.T @ Fe[:, ie] + JJ_rx_i.T @ Te[:, ie]

    # Assemble the link endpoint contributions
    Force_ee = np.block([Fx, Tx, taux])

    # Calculate the accelerations - eq. 3.29
    Ydd = np.linalg.inv(HH) @ (Force + Force_ee - Force0)

    return Ydd


def f_dyn_rk(dt, Y, Yd, Fe, Te, tau, Conn, Prop):
    """
    Integration of the equations of motion using the Runge-Kutta method.

    Given:
    time-step                   dt
    state at time <t>           Y = [R0, A0, q] & Yd = [v0, w0, qd] @ <t>
    external load               Fe, Te, tau

    Returns:
    state at time <t+dt>        Y = [R0, A0, q] & Yd = [v0, w0, qd] @ <t+dt>
    """
    n = len(Y)
    A0 = rpy2dc(Y[3:6]).T
 
    # Step 1
    K1 = np.block([Yd, f_dyn(Y, Yd, Fe, Te, tau, Conn, Prop)])

    # Step 2
    Y_p = Y + dt * K1[:n] / 2.0
    Yd_p = Yd + dt * K1[n:] / 2.0
    A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0     # Correct rotations
    Y_p[3:6] = dc2rpy(A0_p.T)
    K2 = np.block([Yd_p, f_dyn(Y_p, Yd_p, Fe, Te, tau, Conn, Prop)])

    # Step 3
    Y_p = Y + dt * K2[:n] / 2.0
    Yd_p = Yd + dt * K2[n:] / 2.0
    A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0     # Correct rotations
    Y_p[3:6] = dc2rpy(A0_p.T)
    K3 = np.block([Yd_p, f_dyn(Y_p, Yd_p, Fe, Te, tau, Conn, Prop)])

    # Step 4
    Y_p = Y + dt * K3[:n]
    Yd_p = Yd + dt * K3[n:]
    A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0     # Correct rotations
    Y_p[3:6] = dc2rpy(A0_p.T)
    K4 = np.block([Yd_p, f_dyn(Y_p, Yd_p, Fe, Te, tau, Conn, Prop)])

    # Assemble
    Y_p = Y + (K1[:n] + 2.0 * K2[:n] + 2.0 * K3[:n] + K4[:n]) * dt / 6.0
    Yd_p = Yd + (K1[n:] + 2.0 * K2[n:] + 2.0 * K3[n:] + K4[n:]) * dt / 6.0
    A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0     # Correct rotations
    Y_p[3:6] = dc2rpy(A0_p.T)

    return Y_p, Yd_p


def f_dyn_nb(dt, Y, Yd, Fe, Te, tau, Conn, Prop, Params):
    """
    Integration of the equations of motion using the Newmark-beta method and
    the Rodrigues formula to update the rotations.

    Given:
    time-step                   dt
    state at time <t>           Y = [R0, A0, q] & Yd = [v0, w0, qd] @ <t>
    external load               Fe, Te, tau

    Returns:
    state at time <t+dt>        Y = [R0, A0, q] & Yd = [v0, w0, qd] @ <t+dt>

    Notes:
    - gamma = 1/2 guarantees the 2nd order accuracy in the solution.
    - beta = 0 corresponds to the central difference method.
    - beta = 1/4 corresponds to the constant acceleration method.
    - beta = 1/6 corresponds to the linear acceleration method.
    - stability for any <dt> requires beta >= 1/4. ??????? (or is it <= ?)
    """
    # Newmark parameters
    reps, beta, gamma = Params['reps'], Params['beta'], Params['gamma']

    # Accelerations
    Ydd = f_dyn(Y, Yd, Fe, Te, tau, Conn, Prop)

    # Prediction (assume acc_p = acc)
    Y_p = Y + Yd * dt + 0.5 * Ydd * dt ** 2
    Yd_p = Yd + Ydd * dt

    # Correct rotations
    A0 = rpy2dc(Y[3:6]).T
    A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0
    Y_p[3:6] = dc2rpy(A0_p.T)

    # Correction
    for i in range(reps):

        # Accelerations using the predicted state
        Ydd_p = f_dyn(Y_p, Yd_p, Fe, Te, tau, Conn, Prop)

        # Prediction
        Y_p = Y + Yd * dt + ((0.5 - beta) * Ydd + beta * Ydd_p) * dt ** 2
        Yd_p = Yd + ((1.0 - gamma) * Ydd + gamma * Ydd_p) * dt

        # Correct rotations
        A0_p = rpy2dc(Y_p[3:6] - Y[3:6]).T @ A0
        Y_p[3:6] = dc2rpy(A0_p.T)

    return Y_p, Yd_p
