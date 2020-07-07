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
from utils import cross


def build_cc_SS(bodies):
    """
    Builds the matrices cc and SS (connectivity between bodies/joints).

    cc[:, i, k]     Position vector from the centroid of the i-th link to the
                    k-th joint (own joint when k = i)
    SS[i, k]        +1 = the i-th body is connected to the k-th body
                    -1 = i-th body self-connection (always present)
                     0 = no connection between i-th and k-th bodies

    Note: position vectors are given wrt the centroid frame.
    """
    num_b = len(bodies)
    cc = np.zeros((3, num_b, num_b))
    SS = np.zeros((num_b, num_b), dtype=int)
    for i in range(num_b):
        for k, v in bodies[i].cc.items():
            cc[:, i, k] = np.asarray(v)
            SS[i, k] = 1
        SS[i, i] = -1
    return cc, SS


def build_ce_Qe_SE(num_b, ee):
    """
    Builds the matrices ce, Qe, and SE (connectivity between bodies and
    end-points).

    ce[:, i]        Position from the centroid of the i-th link to its end-point
    Qe[:, i]        Orientation of the end-point on the i-th link
    SE[i]           +1 = the i-th body has and end-point
                     0 = the i-th body does not have and end-point

    Note: position and orientation are wrt the centroid frame.
    """
    ce = np.zeros((3, num_b))
    Qe = np.zeros((3, num_b))
    SE = np.zeros(num_b, dtype=int)
    for i, v in ee.items():
        ce[:, i] = np.asarray(v[0:3])           # Position
        Qe[:, i] = np.asarray(v[3:6])           # Orientation
        SE[i] = 1

    return ce, Qe, SE


def build_mass_inertia(bodies):
    """
    Builds the model mass and inertia.

    mass = [mass_1, mass_2, ... ]                       (num_b, )
    inertia = [inertia_1, inertia_2, ... ]              (3, 3*num_b)

    Note: inertia_1 has indeces [:, 0:3], inertia_2 has indeces [:, 3:6], etc.
    """
    num_b = len(bodies)
    mass = np.zeros(num_b)
    inertia = np.zeros((3, 3*num_b))
    for i in range(num_b):
        mass[i] = bodies[i].mass
        inertia[:, 3*i:3*(i+1)] = bodies[i].inertia
    return mass, inertia


def build_BB(SS):
    """
    Builds vector BB.

    BB[i]       Specify the previous link

    Notes:
    - SS is not the full matrix but the sub-matrix [1:num_b, 1:num_b].
    - BB is associated with the joints, thus value BB[0] is for joint 1, value
      BB[1] for joint 2, etc.
    """
    num_j = SS.shape[1]
    BB = np.zeros(num_j, dtype=int)
    for col in range(num_j):
        for row in range(col):
            if (np.abs(SS[row, col]) == 1):
                BB[col] = row + 1
                break
    return BB


def build_j_type(joints):
    """
    Builds the joint type (R or P) list.
    """
    num_j = len(joints)
    j_type = []
    for i in range(num_j):
        j_type.append(joints[i].j_type)
    return j_type


def build_Qi(joints):
    """
    Builds the joint frames.
    """
    num_j = len(joints)
    Qi = np.zeros((3, num_j))
    for i in range(num_j):
        Qi[:, i] = joints[i].Qi
    return Qi


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

    Gravity = np.array([-9.81, 0.0, 0.0])
    Ez = np.array([0.0, 0.0, 1.0])

    def __init__(self, name=None, bodies=[], ee={}, Gravity=Gravity, Ez=Ez):
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
        self.Gravity = np.asarray(Gravity)
        self.Ez = np.asarray(Ez)

        self.joints = bodies[1:]                # List of joints
        self.num_j = len(self.joints)           # Number of joints
        self.num_b = self.num_j + 1             # Number of bodies

        # Model connectivity (Body to body and body to end-point)
        self.cc, self.SS = build_cc_SS(self.bodies)
        self.ce, self.Qe, self.SE = build_ce_Qe_SE(self.num_b, self.ee)

        # Body properties
        self.mass, self.inertia = build_mass_inertia(self.bodies)

        # Joint properties
        self.BB = build_BB(self.SS[1:, 1:])             # Body-pairs connected
        self.j_type = build_j_type(self.joints)
        self.Qi = build_Qi(self.joints)

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

    def set_init(self, R0=np.zeros(3), Q0=np.zeros(3), v0=np.zeros(3),
                 w0=np.zeros(3), q=None, qd=None):
        """
        Initializes quantities at starting time. Unknown are R0, Q0 (A0), v0,
        w0, q, qd.

        Note: the zero refer to the base not to the time.
        """
        # Base position, orientation, and velocities
        R0 = np.asarray(R0)
        Q0 = np.asarray(Q0)
        v0 = np.asarray(v0)
        w0 = np.asarray(w0)

        # Joints rotations/displacements and velocities
        if (q is not None):
            self.q = q
        else:
            self.q = np.zeros(self.num_j)
        if (qd is not None):
            self.qd = qd
        else:
            self.qd = np.zeros(self.num_j)

        # Link rotation matrices (link and joint frame are parallel)
        self.AA = kin.calc_aa(Q0, self.q, self.BB, self.j_type, self.Qi)

        # Centroid positions
        self.RR = kin.calc_pos(R0, self.AA, self.q, self.BB, self.j_type,
                               self.cc, self.Ez)

        # Centroid velocities (linear and angular)
        self.vv, self.ww = kin.calc_vel(self.AA, v0, w0, self.q, self.qd,
                                        self.BB, self.j_type, self.cc, self.Ez)

        # External forces
        self.Fe, self.Te, self.tau = dyn.calc_Forces(self.num_j)

        # Forward dynamics
        # vd0, wd0, self.qdd = kin.f_dyn()
        vd0 = np.array([1.0, 2.0, 3.0])
        wd0 = np.array([0.1, 0.2, 0.3])
        self.qdd = np.array([0.1, 0.5, 1.0, 1.5])

        # Centroid accelerations (linear and angular)
        self.vd, self.wd = kin.calc_acc(self.AA, self.ww, vd0, wd0, self.q,
                                        self.qd, self.qdd, self.BB, self.j_type,
                                        self.cc, self.Ez)
        return

    def calc_CoM(self):
        """
        Returns position, velocity, and acceleration of the system center
        of mass (CoM).
        """
        m = diag(self.mass)
        m_tot = m.sum()

        RR_com = (self.RR @ m).sum() / m_tot
        vv_com = (self.vv @ m).sum() / m_tot
        vd_com = (self.vd @ m).sum() / m_tot

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
