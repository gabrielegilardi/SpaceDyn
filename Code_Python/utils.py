"""
"""

import numpy as np
from math import sin, cos, atan2

def rotX(theta=0.0):
    return np.array([[1.0,         0.0,         0.0],
                     [0.0, +cos(theta), +sin(theta)],
                     [0.0, -sin(theta), +cos(theta)]])


def rotY(theta=0.0):
    return np.array([[+cos(theta), 0.0, -sin(theta)],
                     [        0.0, 1.0,         0.0],
                     [+sin(theta), 0.0, +cos(theta)]])


def rotZ(theta=0.0):
    return np.array([[+cos(theta), +sin(theta), 0.0],
                     [-sin(theta), +cos(theta), 0.0],
                     [        0.0,         0.0, 1.0]])


def tilde(a):
    return np.array([[  0.0, -a[2], +a[1]],
                     [+a[2],   0.0, -a[0]],
                     [-a[1], +a[0],   0.0]])


def rpy2dc(*args):
    if (len(args) == 1):
        rpy = args[0]
        C = rotZ(rpy[2]) @ rotY(rpy[1]) @ rotX(rpy[0])
    else:
        roll = args[0]
        pitch = args[1]
        yaw = args[2]
        C = rotZ(yaw) @ rotY(pitch) @ rotX(roll)
    return C


def eul2dc(*args):
    if (len(args) == 1):
        eul = args[0]
        C = rotZ(eul[2]) @ rotX(eul[1]) @ rotZ(eul[0])
    else:
        phi = args[0]
        theta = args[1]
        psi = args[2]
        C = rotZ(psi) @ rotX(theta) @ rotZ(phi)
    return C


def cross(u, v):
    if ( u.ndim == 1):
        n = np.zeros(3)
        n[0] = u[1] * v[2] - u[2] * v[1]
        n[1] = u[2] * v[0] - u[0] * v[2]
        n[2] = u[0] * v[1] - u[1] * v[0]
    else:
        n = np.zeros([3,u.shape[1]])
        n[0,:] = u[1,:] * v[2,:] - u[2,:] * v[1,:]
        n[1,:] = u[2,:] * v[0,:] - u[0,:] * v[2,:]
        n[2,:] = u[0,:] * v[1,:] - u[1,:] * v[0,:]
    return n


def dc2rpy(C):
    rpy = np.zeros(3)
    eps = np.finfo(float).eps
    if (np.abs(C[1,0]) < eps and np.abs(C[0,0] < eps)):
        rpy[2] = 0.0
        rpy[1] = atan2(C[2,0], C[0,0])
        rpy[0] = atan2(C[1,2], C[1,1])
    else:
        rpy[2] = atan2(-C[1,0], C[0,0])
        c3 = cos(rpy[2])
        s3 = sin(rpy[2])
        rpy[1] = atan2(C[2,0], c3*C[0,0]-s3*C[1,0])
        rpy[0] = atan2(-C[2,1], C[2,2])
    return rpy


def dc2eul(C):
    eul = np.zeros(3)
    eps = np.finfo(float).eps
    if (np.abs(C[0,2]) < eps and np.abs(C[1,2] < eps)):
        eul[2] = 0.0
        eul[1] = atan2(C[1,2], C[2,2])
        eul[0] = atan2(C[0,1], C[0,0])
    else:
        eul[2] = atan2(C[0,2], C[1,2])
        c3 = cos(eul[2])
        s3 = sin(eul[2])
        eul[1] = atan2(s3*C[0,2]+c3*C[1,2], C[2,2])
        eul[0] = atan2(C[2,0], -C[2,1])
    return eul


def tr2diff(AA1, AA2):
    a = np.array(cross(AA2[:,0],AA1[:,0]) + cross(AA2[:,1],AA1[:,1]) \
        + cross(AA2[:,2],AA1[:,2])) / 2.0
    return a


def rotW(w0, dtime):
    nw0 = np.linalg.norm(w0) 
    if ( nw0 == 0.0 ):
        E0 = np.identity(3)
    else:
        theta = nw0 * dtime
        w = w0 / nw0
        ct = cos(theta)
        st = sin(theta)
        E0 = np.array([[        ct+(1-ct)*w[0]**2, (1-ct)*w[0]*w[1]-st*w[2], (1-ct)*w[0]*w[2]+st*w[1] ],
                       [ (1-ct)*w[0]*w[1]+st*w[2],        ct+(1-ct)*w[1]**2, (1-ct)*w[1]*w[2]-st*w[0] ],
                       [ (1-ct)*w[0]*w[2]-st*w[1], (1-ct)*w[1]*w[2]+st*w[0],        ct+(1-ct)*w[2]**2 ]])
    return E0


def inertia_matrix(shape):
    """Compute the moments of inertia per unit of mass wrt the link CM [m^2]"""
    In = np.zeros((3,3))
    # Thick-walled cylindrical tube: Re, Ri, H [m]
    # Set Ri = 0 for a solid cylindrical tube
    if (shape[0] == 'Cylinder'):
        Re = shape[1]
        Ri = shape[2]
        H = shape[3]
        In[0,0] = (Re*Re + Ri*Ri) / 2.0
        In[1,1] = (3*Re*Re + 3*Ri*Ri + H*H) / 12.0
        In[2,2] = (3*Re*Re + 3*Ri*Ri + H*H) / 12.0
    # Thick-walled sphere: Re, Ri [m]
    # Set Ri = 0 for a solid sphere
    elif (shape[0] == 'Sphere'):
        Re = shape[1]
        Ri = shape[2]
        In[0,0] = (2.0/5.0) * (Re**5 - Ri**5) / (Re**3 - Ri**3)
        In[1,1] = (2.0/5.0) * (Re**5 - Ri**5) / (Re**3 - Ri**3)
        In[2,2] = (2.0/5.0) * (Re**5 - Ri**5) / (Re**3 - Ri**3)
    # Default is the unit matrix
    else:
        In[0,0] = 1.0
        In[1,1] = 1.0
        In[2,2] = 1.0
    return In


def calc_Momentum(RR, AA, vv , ww, mass, inertia, P_ref=np.zeros(3)):
    """Return the linear and angular momentum"""

    # Linear momentum
    LM = mass * vv

    # Angular momentum
    n = len(mass)
    HM = np.zeros((3,n))
    for i in range(n):
        A = AA[:,3*i:3*(i+1)]
        In =  inertia[:,3*i:3*(i+1)]
        HM[:,i] = cross( (RR[:,i] - P_ref), LM[:,i] ) + A @ In @ A.T @ ww[:,i]

    return LM, HM
