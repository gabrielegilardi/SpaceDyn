"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np


def rotX(theta=0.0):
    """
    Rotation about axis X.
    """
    Rx = np.array([[1.0,            0.0,            0.0],
                   [0.0,  np.cos(theta),  np.sin(theta)],
                   [0.0, -np.sin(theta),  np.cos(theta)]])

    return Rx


def rotY(theta=0.0):
    """
    Rotation about axis Y.
    """
    Ry = np.array([[ np.cos(theta), 0.0, -np.sin(theta)],
                   [           0.0, 1.0,            0.0],
                   [ np.sin(theta), 0.0,  np.cos(theta)]])

    return Ry


def rotZ(theta=0.0):
    """
    Rotation about axis Z.
    """
    Rz = np.array([[ np.cos(theta),  np.sin(theta), 0.0],
                   [-np.sin(theta),  np.cos(theta), 0.0],
                   [           0.0,            0.0, 1.0]])

    return Rz


def tilde(a):
    """
    Converts a vector to the skew-symmetric matrix.
    """
    B = np.array([[  0.0, -a[2],  a[1]],
                  [ a[2],   0.0, -a[0]],
                  [-a[1],  a[0],   0.0]])

    return B


def rpy2dc(*args):
    """
    Converts RPY (XYZ) angles to the direction cosine matrix.
    """
    # Passed as single argument (array of three elements)
    if (len(args) == 1):
        rpy = args[0]
        C = rotZ(rpy[2]) @ rotY(rpy[1]) @ rotX(rpy[0])

    # Passed as three separate arguments
    else:
        roll = args[0]
        pitch = args[1]
        yaw = args[2]
        C = rotZ(yaw) @ rotY(pitch) @ rotX(roll)

    return C


def dc2rpy(C):
    """
    Converts a direction cosine matrix to RPY (XYZ) angles.
    """
    rpy = np.zeros(3)
    eps = np.finfo(float).eps

    # Singular solution
    if (np.abs(C[1, 0]) < eps and np.abs(C[0, 0] < eps)):
        rpy[2] = 0.0
        rpy[1] = np.atan2(C[2, 0], C[0, 0])
        rpy[0] = np.atan2(C[1, 2], C[1, 1])

    # Non singular solution
    else:
        rpy[2] = np.atan2(-C[1, 0], C[0, 0])
        c3 = np.cos(rpy[2])
        s3 = np.sin(rpy[2])
        rpy[1] = np.atan2(C[2, 0], c3 * C[0, 0] - s3 * C[1, 0])
        rpy[0] = np.atan2(-C[2, 1], C[2, 2])

    return rpy


def eul2dc(*args):
    """
    Converts Eulero (ZXZ) angles to the direction cosine matrix.
    """
    # Passed as single argument (array of three elements)
    if (len(args) == 1):
        eul = args[0]
        C = rotZ(eul[2]) @ rotX(eul[1]) @ rotZ(eul[0])

    # Passed as three separate arguments
    else:
        phi = args[0]
        theta = args[1]
        psi = args[2]
        C = rotZ(psi) @ rotX(theta) @ rotZ(phi)

    return C


def dc2eul(C):
    """
    Converts a direction cosine matrix to Eulero (ZXZ) angles.
    """
    eul = np.zeros(3)
    eps = np.finfo(float).eps

    # Singular solution
    if (np.abs(C[0, 2]) < eps and np.abs(C[1, 2] < eps)):
        eul[2] = 0.0
        eul[1] = atan2(C[1, 2], C[2, 2])
        eul[0] = atan2(C[0, 1], C[0, 0])

    # Non singular solution
    else:
        eul[2] = atan2(C[0, 2], C[1, 2])
        c3 = cos(eul[2])
        s3 = sin(eul[2])
        eul[1] = atan2(s3 * C[0, 2] + c3 * C[1, 2], C[2, 2])
        eul[0] = atan2(C[2, 0], -C[2, 1])

    return eul


def cross(u, v):
    """
    Cross product between two vectors
    """
    # Arguments are two vectors
    if (u.ndim == 1):
        n = np.zeros(3)
        n[0] = u[1] * v[2] - u[2] * v[1]
        n[1] = u[2] * v[0] - u[0] * v[2]
        n[2] = u[0] * v[1] - u[1] * v[0]

    # Arguments are two matrices with shape (3, :) --> cross product between
    # corresponding columns
    else:
        n = np.zeros([3, u.shape[1]])
        n[0, :] = u[1, :] * v[2, :] - u[2, :] * v[1, :]
        n[1, :] = u[2, :] * v[0, :] - u[0, :] * v[2, :]
        n[2, :] = u[0, :] * v[1, :] - u[1, :] * v[0, :]

    return n


def tr2diff(A1, A2):
    """
    Calculates the difference between two rotation matrices (RPY angles).

    Notes:
    - correct only for small angles.
    - fail if the angle difference is 180 degrees.
    """
    b = (cross(A2[:, 0], A1[:, 0]) + cross(A2[:, 1], A1[:, 1])
        + cross(A2[:, 2], A1[:, 2])) / 2.0

    return b


def rotW(w0, dt):
    """
    Returns a transformation matrix representing a rotation <dt> about the
    vector w0.
    """
    nw0 = np.linalg.norm(w0)
    eps = np.finfo(float).eps

    if (np.abs(nw0) < eps):
        E0 = np.identity(3)         # No rotation

    else:
        theta = nw0 * dt
        w = w0 / nw0
        ct = np.cos(theta)
        st = np.sin(theta)
        E0 = np.array([[        ct+(1-ct)*w[0]**2, (1-ct)*w[0]*w[1]-st*w[2], (1-ct)*w[0]*w[2]+st*w[1] ],
                       [ (1-ct)*w[0]*w[1]+st*w[2],        ct+(1-ct)*w[1]**2, (1-ct)*w[1]*w[2]-st*w[0] ],
                       [ (1-ct)*w[0]*w[2]-st*w[1], (1-ct)*w[1]*w[2]+st*w[0],        ct+(1-ct)*w[2]**2 ]])

    return E0


def inertia_matrix(shape, *args):
    """
    Returns the moments of inertia wrt the centroid for basic shapes.
    """
    inertia = np.zeros((3, 3))

    # Thick-walled cylindrical tube: Re, Ri, H
    # Set Ri = 0 for a solid cylindrical tube
    if (shape == 'Cylinder'):
        Re = args[0]
        Ri = args[1]
        H = args[2]
        inertia[0, 0] = (Re * Re + Ri * Ri) / 2.0
        inertia[1, 1] = (3.0 * Re * Re + 3.0 * Ri * Ri + H * H) / 12.0
        inertia[2, 2] = inertia[1, 1]

    # Thick-walled sphere: Re, Ri
    # Set Ri = 0 for a solid sphere
    elif (shape == 'Sphere'):
        Re = args[0]
        Ri = args[1]
        inertia[0, 0] = 0.4 * (Re ** 5.0 - Ri ** 5.0) / (Re ** 3.0 - Ri ** 3.0)
        inertia[1, 1] = inertia[0, 0]
        inertia[2, 2] = inertia[0, 0]

    # Default is the unit matrix
    else:
        inertia = np.eye(3)

    return inertia
