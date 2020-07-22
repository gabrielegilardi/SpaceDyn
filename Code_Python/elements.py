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
from utils import cross, rpy2dc
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

    def __init__(self, name=None, bodies=[], ee={}):
        """
        name        str         Name
        bodies      num_b       List of bodies
        ee          dict        Position/orientation of end-points
        Gravity     (3, )       Gravity
        Ez          (3, )       Joint axis
        """
        self.name = name
        self.bodies = bodies
        self.ee = ee

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
        print()
        print('CONNECTIVITY')
        print('SS =')
        print(self.SS)
        print('SE =')
        print(self.SE)
        print('BB =')
        print(self.BB)
        print()

        return

    def set_init(self, R0=np.zeros(3), A0=np.eye(3), v0=np.zeros(3),
                 w0=np.zeros(3), q=np.array([]), qd=np.array([])):
        """
        Initializes quantities at starting time.

        Note: the zero refer to the base not to the time.
        """
        # Base (position, orientation, and velocities)
        R0 = np.asarray(R0)
        A0 = np.asarray(A0)
        v0 = np.asarray(v0)
        w0 = np.asarray(w0)

        # Joints (rotations/displacements and velocities)
        if (len(q) > 0):
            self.q = q
        else:
            self.q = np.zeros(self.num_j)
        if (len(qd)):
            self.qd = qd
        else:
            self.qd = np.zeros(self.num_j)

        # Rotation matrices (link and joint frame are parallel)
        self.AA = kin.calc_aa(A0, self.q, self.BB, self.j_type, self.Qi)

        # Centroid positions
        self.RR = kin.calc_pos(R0, self.AA, self.q, self.BB, self.j_type,
                               self.cc)

        # Centroid velocities (linear and angular)
        self.vv, self.ww = kin.calc_vel(self.AA, v0, w0, self.q, self.qd,
                                        self.BB, self.j_type, self.cc)

        # # External forces
        # self.Fe, self.Te = user.calc_forces(self.num_j)

        # # Forward dynamics
        # vd0 = np.array([-1.7, 2.4, -4.5])
        # wd0 = np.array([0.3, -0.2, 0.13])
        # qdd = np.array([0.1, -0.3, 0.6, -1.1])
        # Force = dyn.r_ne(self.RR, self.AA, v0, w0, vd0, wd0, q, qd, qdd, self.Fe, self.Te,
        #                 self.SS, self.SE, self.j_type, self.cc, self.ce,
        #                 self.mass, self.inertia, self.BB)

        # # vd0, wd0, self.qdd = kin.f_dyn()

        # Centroid accelerations (linear and angular)
        # self.vd, self.wd = kin.calc_acc(self.AA, self.ww, vd0, wd0, self.q,
        #                                 self.qd, self.qdd, self.BB, self.j_type,
        #                                 self.cc)

    def calc_CoM(self):
        """
        Returns position, velocity, and acceleration of the system center
        of mass (CoM).
        """
        m = diag(self.mass)
        m_tot = m.sum()

        # Sum along axis 1?????
        RR_com = (self.RR @ m).sum(axis=1) / m_tot
        vv_com = (self.vv @ m).sum(axis=1) / m_tot
        vd_com = (self.vd @ m).sum(axis=1) / m_tot

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
        VG = - self.mass * (self.Gravity.T @ self.RR)

        return VG.sum(), VG

    def calc_work(self):
        """
        Returns the work done by external forces/torques.
        """
        pass

# function [ PF0, PFe, Ptau ] = calc_Work( v0, w0, vve, ww, qd, F0, T0, Fe, ...
#                                     Te, tau )

#   num_q = length(qd);
  
#   %Work forces/moments on the base
#   PF0 = F0'*v0 + T0'*w0;
  
#   %Work forces/moments on the end effector
#   PFe = zeros(1,num_q);
#   for i = 1:num_q
#     PFe(1,i) = Fe(:,i)'*vve(:,i) + Te(:,i)'*ww(:,i);
#   end
  
#   %Work forces/moments on the joints
#   Ptau = zeros(1,num_q);
#   for i = 1:num_q
#     Ptau(1,i) = tau(i)*qd(i);
#   end

# end     % End of function


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
                       self.mass[i] * self.vv[:, i] + In @ self.ww[:, i])

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
