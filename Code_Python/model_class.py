
# Libraries and external files

import numpy as np
import kinematics as kin
import dynamics as dyn
from utils import calc_Momentum


# Helper functions

def build_cc_SS(bodies):
    """Build matrices cc and SS (used for base and links)"""
    num_b = len(bodies)
    cc = np.zeros((3,num_b,num_b))
    SS = np.zeros((num_b,num_b), dtype=int)
    for i in range(num_b):
        for k, v in bodies[i].cc.items():
            cc[:,i,k] = np.asarray(v)
            SS[i,k] = 1
        SS[i,i] = -1
    return cc, SS

def build_ce_SE(bodies):
    """Build matrices ce and SE (used for base and links). Endpoint for the base is always
       in its centroid"""
    num_b = len(bodies)
    ce = np.zeros((3,num_b))            
    SE = np.zeros(num_b, dtype=int)     
    SE[0] = 1           # Value for base is always one
    for i in range(1,num_b):
        for k, v in bodies[i].ce.items():
            ce[:,k] = np.asarray(v)
            SE[k] = 1
    return ce, SE

def build_mass_inertia(bodies):
    """Build the system mass and inertia (used for base and links)"""
    num_b = len(bodies)
    mass = np.zeros(num_b)
    inertia = np.zeros((3,3*num_b))
    for i in range(num_b):
        mass[i] = bodies[i].mass
        inertia[:,3*i:3*(i+1)] = bodies[i].inertia
    return mass, inertia

def build_Qe(bodies):
    """Build endpoint frame (used for base and links). Endpoint frame for the base is
    always coincident with the centroid frame"""
    num_b = len(bodies)
    Qe = np.zeros((3,num_b))        # Value for base is always zero
    for i in range(1,num_b):
        Qe[:,i] = bodies[i].Qe
    return Qe

def build_BB(SS):
    """Build vector BB (used for joints)"""
    num_j = SS.shape[1]
    BB = np.zeros(num_j, dtype=int)
    for col in range(0,num_j):
        for row in range(col):
            if (np.abs(SS[row,col]) == 1):
                BB[col] = row + 1
                break
    return BB

def build_j_type(joints):
    """Build joint type (used for joints)"""
    num_j = len(joints)
    j_type = []             
    for i in range(num_j):
        j_type.append(joints[i].j_type)
    return j_type

def build_Qi(joints):
    """Build joint frame (used for joints)"""
    num_j = len(joints)
    Qi = np.zeros((3,num_j))
    for i in range(num_j):
        Qi[:,i] = joints[i].Qi
    return Qi


# Class model
class model:

    # Default values for some of the parameters
    Gravity = np.array([-9.81,0.0,0.0])             # Gravity
    Ez = np.array([0.0, 0.0, 1.0])                  # Joint z-axis

    def __init__(self, bodies, joints=[]):
        """Build model and initialize properties"""

        # There is one joint for each link
        # Bodies are the links plus the base
        self.num_j = len(joints)                # Number of joints/links
        self.num_b = self.num_j + 1             # Number of bodies
        self.bodies = bodies
        self.joints = joints

        # Connectivity quantities for bodies (0 to num_b)
        self.cc, self.SS = build_cc_SS(bodies)      # Between bodies
        self.ce, self.SE = build_ce_SE(bodies)      # Between bodies and endpoints

        # Body properties (0 to num_b)
        self.mass, self.inertia = build_mass_inertia(bodies) 
        self.Qe = build_Qe(bodies) 

        # Joint properties (0 to num_j)
        self.BB = build_BB(self.SS[1:,1:])          # Body pairs connected by each joint
        self.j_type = build_j_type(joints)
        self.Qi = build_Qi(joints)


    def set_param(self, Gravity=Gravity, Ez=Ez):
        """Set some the system parameters if different from default value"""
        self.Gravity = np.asarray(Gravity)
        self.Ez = np.asarray(Ez)
        

    def set_init(self, R0=np.zeros(3), Q0=np.zeros(3), v0=np.zeros(3), w0=np.zeros(3), \
        q=None, qd=None):
        """Initialize quantities at starting time. The zeros refer to the base, not to the
           starting time. Unknown: R0, Q0 (A0), v0, w0, q, qd"""
        
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

        # Rotation matrices of all bodies (link centroid frame is equal to joint frame)
        self.AA = kin.calc_aa(Q0, self.q, self.BB, self.j_type, self.Qi)

        # Position vector of all body centroids
        self.RR = kin.calc_pos(R0, self.AA, self.q, self.BB, self.j_type, self.cc, self.Ez)

        # Linear velocity vector of all body centroids and angular velocity vector of
        # all bodies
        self.vv, self.ww = kin.calc_vel(self.AA, v0, w0, self.q, self.qd, self.BB, \
            self.j_type, self.cc, self.Ez)

        # External forces
        self.Fe, self.Te, self.tau = dyn.calc_Forces(self.num_j)

        P_ref = np.zeros(3)
        self.LM = calc_Momentum(self.RR, self.AA, self.vv, self.ww, self.mass, self.inertia, \
            P_ref)

        # Forward dynamics
        # vd0, wd0, self.qdd = kin.f_dyn()

        # Linear acceleration vector of all body centroids and angular acceleration vector
        # of all bodies
        # self.vd, self.wd = kin.calc_acc(self.AA, self.ww, vd0, wd0, self.q, self.qd, \
        #     self.qdd, self.BB, self.j_type, self.cc, self.Ez)



        #===============================

        # # Check calc_acc
        # vd0 = np.array([0.2, 0.2, 0.2])
        # wd0 = np.array([0.3, 0.3, 0.3])
        # qdd = np.array([-0.1, -0.2, -0.3, -0.4])
        # qdd = np.array([])
        # vd, wd = kin.calc_acc(self.AA, self.ww, vd0, wd0, self.q, self.qd, qdd, self.BB, \
        #     self.j_type, self.cc, self.Ez)
        
        # Check calc_je, calc_jte, calc_jre, f_kin_j, f_kin_e, j_num
        # seq_link = kin.j_num(3, self.BB)
        # Jacobian = kin.calc_je(self.RR, self.AA, self.q, seq_link, self.j_type, self.cc, \
        #     self.ce, self.Qe, self.Ez)

        # Check calc_jt, calc_jr, calc_hh
        # HH = kin.calc_hh(self.RR, self.AA, self.mass, self.inertia, self.BB, \
        #     self.j_type, self.cc, self.Ez)

        # Check r_ne (with accelerations and external forces equal to zero)
        # vd0 = np.zeros(3)
        # wd0 = np.zeros(3)
        # qdd = np.zeros(self.num_j)
        # Fe, Te, tau = dyn.calc_Forces(self.num_j)
        # Force = dyn.r_ne(self.RR, self.AA, v0, w0, vd0, wd0, self.q, self.qd, qdd, Fe, Te, \
        #     self.SS, self.SE, self.j_type, self.cc, self.ce, self.mass, self.inertia, \
        #         self.Ez, self.Gravity, self.BB)


    def update_elements(self):

        # Update body quantities
        for i in range(self.num_b):
            self.bodies[i].AA = self.AA[:,3*i:3*(i+1)]
            self.bodies[i].RR = self.RR[:,i]
            self.bodies[i].vv = self.vv[:,i]
            self.bodies[i].ww = self.ww[:,i]

        # Update joint quantities
        for i in range(self.num_j):
            self.joints[i].q = self.q[i]
            self.joints[i].qd = self.qd[i]

