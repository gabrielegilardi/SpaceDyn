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
import dynamics as dyn
from utils import cross, rpy2dc, dc2rpy
import user as user


def connectivity(bodies, ee):
    """
    Builds the connectivity matrices.

    SS[i, k]        +1 = the i-th body is connected to the k-th body
                    -1 = i-th body self-connection (always present)
                     0 = no connection between i-th and k-th bodies
    SE[i]           list of bodies with endpoints
    BB[i]           lower connection of each link (base is excluded)

    Note: no own joint or previous body for the base.
    """
    num_b = len(bodies)
    num_e = len(ee)

    # Base-to-link and link-to-link connections
    SS = np.zeros((num_b, num_b), dtype=int)
    for i in range(num_b):
        for k, v in bodies[i].cc.items():
            SS[i, k] = 1
        SS[i, i] = -1               # Value for the base never used

    # Base-to-endpoint and link-to-endpoint connections
    SE = np.zeros(num_e, dtype=int)
    for i, v in ee.items():
        SE[i] = v[0]

    # Previous body connection (BB[0] for body 1, BB[1] for body 2, etc.)
    BB = np.zeros(num_b-1, dtype=int)
    for col in range(1, num_b):
        for row in range(1, col):
            if (np.abs(SS[row, col]) == 1):
                BB[col-1] = row
                break

    return SS, SE, BB


def build_cc_Qi(bodies):
    """
    Builds the matrix defining the links/joints positions and orientations.

    cc[:, i, k]     Position from the centroid of the i-th body to the joint
                    on the k-th body (own joint when k = i)

    Notes:
    - positions and orientations are given wrt the centroid frame.
    - the base does not have a relative orientation.
    """
    num_b = len(bodies)

    # Relative positions (base + links)
    cc = np.zeros((3, num_b, num_b))
    for i in range(num_b):
        for k, v in bodies[i].cc.items():
            cc[:, i, k] = np.asarray(v)

    # Relative orientations (links only)
    Qi = np.zeros((3, num_b-1))
    for i in range(1, num_b):
        Qi[:, i-1] = bodies[i].Qi

    return cc, Qi


def build_ce_Qe(ee):
    """
    Builds the matrices defining the endpoint positions and orientations.

    ce[:, i]        Position from the centroid of the SE[i]-th body to the
                    i-th endpoint
    Qe[:, i]        Orientation of the end-point on the i-th body

    Note: position and orientation are wrt the centroid frame.
    """
    num_e = len(ee)
    ce = np.zeros((3, num_e))
    Qe = np.zeros((3, num_e))
    for i, v in ee.items():
        ce[:, i] = np.asarray(v[1])           # Position
        Qe[:, i] = np.asarray(v[2])           # Orientation

    return ce, Qe


def build_mass_inertia(bodies):
    """
    Builds the model mass and inertia.

    mass = [mass_1, mass_2, ... ]                       (num_b, )
    inertia = [inertia_1, inertia_2, ... ]              (3, 3*num_b)

    Note: the inertia matrices are with respect to the body centroids.
    """
    num_b = len(bodies)
    mass = np.zeros(num_b)
    inertia = np.zeros((3, 3*num_b))
    for i in range(num_b):
        mass[i] = bodies[i].mass
        inertia[:, 3*i:3*(i+1)] = bodies[i].inertia
    return mass, inertia


def build_j_type(joints):
    """
    Builds the joint type (R or P) list.
    """
    num_j = len(joints)
    j_type = []
    for i in range(num_j):
        j_type.append(joints[i].j_type)
    return j_type


class base:

    def __init__(self, name=None, mass=0.0, inertia=np.eye(3), cc={}):
        """
        name        str         Name
        mass        scalar      Mass
        inertia     (3, 3)      Inertia matrix
        cc          dict        Position of connected links

        Notes:
        - no joint for the base, so no self-connection in <cc>.
        - orientation of its centroid is an unknown of the problem.
        """
        self.name = name
        self.mass = mass
        self.inertia = np.asarray(inertia)
        self.cc = cc

    def info(self):
        """
        """
        print()
        print('Type: base')
        print('Name = ', self.name)
        print('Mass = ', self.mass)
        print('Inertia =')
        print(self.inertia)
        print('Connections =')
        print(self.cc)
        print()

        return


class link:

    def __init__(self, name=None, mass=0.0, inertia=np.eye(3), j_type='R',
                 Qi=np.zeros(3), cc={}):
        """
        name        str         Name
        mass        scalar      Mass [kg]
        inertia     (3, 3)      Inertia matrix
        j_type      str         Joint type [R or P]
        Qi          (3, )       Link/joint relative orientation
        cc          dict        Position of connected bodies

        Notes:
        - the link and joint frames are always parallel.
        - the frame is given with respect to the previous link/joint frame.
        - joint quantities can be linear (for a prismatic joint) or angular
          (for a rotational joint)
        """
        self.name = name
        self.mass = mass
        self.inertia = np.asarray(inertia)
        self.j_type = j_type
        self.cc = cc
        self.Qi = np.asarray(Qi)

    def info(self):
        """
        """
        print()
        print('Type: link')
        print('Name = ', self.name)
        print('Mass = ', self.mass)
        print('Inertia =')
        print(self.inertia)
        print('Connections =')
        print(self.cc)
        print('Frame orientation = ', self.Qi)
        print('Joint type = ', self.j_type)
        print()

        return


class model:

    def __init__(self, name=None, bodies=[], ee={}, load=None):
        """
        name        str         Name
        bodies      num_b       List of bodies
        ee          dict        Position/orientation of end-points
        load        str         Load type
        """
        self.name = name
        self.bodies = bodies
        self.ee = ee
        self.load = load

        self.joints = self.bodies[1:]           # List of joints
        self.num_j = len(self.joints)           # Number of joints
        self.num_b = self.num_j + 1             # Number of bodies
        self.num_e = len(self.ee)

        # Connectivity
        self.SS, self.SE, self.BB = connectivity(self.bodies, self.ee)

        # Links and endpoints relative positions/orientations
        self.cc, self.Qi = build_cc_Qi(self.bodies)
        self.ce, self.Qe = build_ce_Qe(self.ee)

        # Properties
        self.mass, self.inertia = build_mass_inertia(self.bodies)
        self.j_type = build_j_type(self.joints)

    def info(self):
        """
        """
        print()
        print('MODEL')
        print('Name = ', self.name)
        print('Number of bodies = ', self.num_b)
        print('Number of joints = ', self.num_j)
        print('Number of endpoints = ', self.num_e)
        print()
        print('CONNECTIVITY')
        print('SS =')
        print(self.SS)
        print('SE = ', self.SE)
        print('BB = ', self.BB)
        print('joint type = ', self.j_type)
        print()
        print('PROPERTIES')
        print('Mass = ', self.mass)
        print('Inertia =')
        print(self.inertia)
        print('Body-to-body connections =')
        print(self.cc)
        print('Body orientations =')
        print(self.Qi)
        print('Body-to-endpoint connections =')
        print(self.ce)
        print('Endpoint orientations =')
        print(self.Qe)

        return

    def simulate(self, ts=0.0, tf=1.0, dt=0.001, rec=None, load=None):
        """
        """
        R0 = self.R0
        A0 = self.A0
        v0 = self.v0
        w0 = self.w0
        q = self.q
        qd = self.qd
        n_steps = int(np.round((tf - ts) / dt)) + 1
        n_rec = int(np.round((tf - ts) / rec)) + 1
        for i in range(1, n_steps):
            t = ts + float(i) * dt
            # print('time = ', t)
            Fe, Te, tau = user.calc_forces(t, self.num_j, self.num_e, self.load)
            R0, A0, v0, w0, q, qd = \
                dyn.f_dyn_rk(dt, R0, A0, v0, w0, q, qd, Fe, Te, tau, self.SS,
                             self.SE, self.BB, self.j_type, self.cc, self.ce,
                             self.mass, self.inertia, self.Qi, self.Qe)

        # Results
        vd0, wd0, qdd = dyn.f_dyn(R0, A0, v0, w0, q, qd, Fe, Te, tau, self.SS,
                                self.SE, self.BB, self.j_type, self.cc, self.ce,
                                self.mass, self.inertia, self.Qi, self.Qe)
        print('tf = ', t)
        print('R0 = ', R0)
        print('Q0 = ', dc2rpy(A0.T))
        print('v0 = ', v0)
        print('w0 = ', w0)
        print('vd0 = ', vd0)
        print('wd0 = ', wd0)
        print('q =', q)
        print('qd =', qd)
        print('qdd =', qdd)

    def set_init(self, R0=np.zeros(3), A0=np.eye(3), v0=np.zeros(3),
                 w0=np.zeros(3), q=np.array([]), qd=np.array([])):
        """
        Initializes quantities at the starting time.

        Note: the zero refer to the base not to the time.
        """
        # Base (position, orientation, and velocities)
        self.R0 = np.asarray(R0)
        self.A0 = np.asarray(A0)
        self.v0 = np.asarray(v0)
        self.w0 = np.asarray(w0)
        if (len(q) > 0):
            self.q = q
        else:
            self.q = np.zeros(self.num_j)
        if (len(qd)):
            self.qd = qd
        else:
            self.qd = np.zeros(self.num_j)

        # # Rotation matrices (link and joint frame are parallel)
        # self.AA = kin.calc_aa(A0, self.q, self.BB, self.j_type, self.Qi)

        # # Centroid positions
        # self.RR = kin.calc_pos(R0, self.AA, self.q, self.BB, self.j_type,
        #                        self.cc)

        # # Centroid linear and angular velocities
        # self.vv, self.ww = kin.calc_vel(self.AA, v0, w0, self.q, self.qd,
        #                                 self.BB, self.j_type, self.cc)

        # # External forces
        # Fe, Te, tau = user.calc_forces(t0, self.num_j, self.num_e, self.load)

        # # Forward dynamics
        # vd0, wd0, qdd = dyn.f_dyn(R0, A0, v0, w0, self.q, self.qd, Fe, Te,
        #                           tau, self.SS, self.SE, self.BB, self.j_type,
        #                           self.cc, self.ce, self.mass, self.inertia,
        #                           self.Qi, self.Qe)

        # # Centroid linear and angular accelerations
        # self.vd, self.wd = kin.calc_acc(self.AA, self.ww, vd0, wd0, self.q,
        #                                 self.qd, qdd, self.BB, self.j_type,
        #                                 self.cc)

    def calc_CoM(self):
        """
        Returns position, velocity, and acceleration of the system center of
        mass (CoM).
        """
        mt = self.mass.sum()

        RR_com = (self.RR * self.mass).sum(axis=1) / mt
        vv_com = (self.vv * self.mass).sum(axis=1) / mt
        vd_com = (self.vd * self.mass).sum(axis=1) / mt

        return RR_com, vv_com, vd_com

    def calc_kin_ener(self):
        """
        Returns the kinetic energy for the entire system and for each body.
        """
        TK = np.zeros(self.num_b)
        for i in range(self.num_b):
            AA = self.AA[:, 3*i:3*(i+1)]
            In = AA @ self.inertia[:, 3*i:3*(i+1)] @ AA.T
            TK[i] = self.mass[i] * (self.vv[:, i].T @ self.vv[:, i]) / 2.0 \
                    + self.ww[:, i].T @ In @ self.ww[:, i] / 2.0

        return TK.sum(), TK

    def calc_pot_ener(self):
        """
        Returns the potential energy for the entire system and for each body.
        """
        Gravity = np.array([0.0, 0.0, -9.81])       # Gravity vector

        VG = - self.mass * (Gravity.T @ self.RR)

        return VG.sum(), VG

    def calc_work(self, time=0.0):
        """
        Returns the work done by all external forces and moments (total and
        components - on the base, on the joints, on the endpoints).
        """
        # Evaluate external forces and moments
        Fe, Te, tau = user.calc_forces(time, self.num_j, self.num_e, self.load)

        # Work done by the forces/moments on the endpoints (links and base)
        WK0 = 0      # will be added later
        WKe = 0      # will be added later
        for ie in range(self.num_e):
            if (self.SE[ie] == 0):
                pass        # base
            else:
                pass        # link

        # Work done by the control forces/torques on the joints
        WKq = (tau * self.qd).sum()

        WK = np.array([WK0, WKq, WKe])

        return WK.sum(), WK

    def calc_lin_mom(self):
        """
        Returns the linear momentum for the entire system and for each body.
        """
        LM = self.mass * self.vv

        return LM.sum(axis=1), LM

    def calc_lin_mom1(self):
        """
        Returns the derivative of the linear momentum for the entire system and
        for each body.
        """
        LM1 = self.mass * self.vd

        return LM1.sum(axis=1), LM1

    def calc_ang_mom(self, P_ref=np.zeros(3)):
        """
        Returns the angular momentum for the entire system and for each body
        with respect to point P_ref.
        """
        HM = np.zeros((3, self.num_b))
        for i in range(self.num_b):
            AA = self.AA[:, 3*i:3*(i+1)]
            In = AA @ self.inertia[:, 3*i:3*(i+1)] @ AA.T
            HM[:, i] = cross(self.RR[:, i] - P_ref,
                       self.mass[i] * self.vv[:, i]) + In @ self.ww[:, i]

        return HM.sum(axis=1), HM

    def calc_ang_mom1(self, P_ref=np.zeros(3), V_ref=np.zeros(3)):
        """
        Returns the derivative of the angular momentum for the entire system
        and for each body with respect to point P_ref and velocity V_ref.
        """
        HM = np.zeros((3, self.num_b))
        for i in range(self.num_b):
            AA = self.AA[:, 3*i:3*(i+1)]
            In = AA @ self.inertia[:, 3*i:3*(i+1)] @ AA.T
            HM[:, i] = In @ self.wd[:, i] \
                + cross(self.ww[:, i], In @ self.ww[:, i]) \
                + cross(self.RR[:, i] - P_ref, self.mass[i] * self.vd[:, i]) \
                + cross(self.vv[:, i] - V_ref, self.mass[i] * self.vv[:, i])

        return HM.sum(axis=1), HM
