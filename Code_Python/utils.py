"""
SpaceDyn - A Toolbox for Space, Mobile, and Humanoid Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np


def rotX(alpha=0.0):
    """
    Returns the rotation matrix for a rotation <alpha> about axis X.
    """
    Rx = np.array([[1.0,            0.0,            0.0],
                   [0.0,  np.cos(alpha),  np.sin(alpha)],
                   [0.0, -np.sin(alpha),  np.cos(alpha)]])

    return Rx


def rotY(beta=0.0):
    """
    Returns the rotation matrix for a rotation <beta> about axis Y.
    """
    Ry = np.array([[ np.cos(beta), 0.0, -np.sin(beta)],
                   [          0.0, 1.0,           0.0],
                   [ np.sin(beta), 0.0,  np.cos(beta)]])

    return Ry


def rotZ(gamma=0.0):
    """
    Returns the rotation matrix for a rotation <gamma> about axis Z.
    """
    Rz = np.array([[ np.cos(gamma),  np.sin(gamma), 0.0],
                   [-np.sin(gamma),  np.cos(gamma), 0.0],
                   [           0.0,            0.0, 1.0]])

    return Rz


def rpy2dc(rpy):
    """
    Returns the rotation matrix for a set of RPY angles.

    Note: RPY (roll-pitch-yaw) angles corresponds to the X-Y-Z fixed-angles
          representation.
    """
    rpy = np.asarray(rpy).flatten()
    R = rotZ(rpy[2]) @ rotY(rpy[1]) @ rotX(rpy[0])

    return R


def dc2rpy(R):
    """
    Returns the RPY angles for the specified rotation matrix.

    Notes:
    - in the singular case (where there are an infinite number of solution) the
      solution with yaw angle equal to zero is chosen.
    - in the general case (where there are always two solutions) the positive
      solution for the pitch angle is chosen.
    """
    rpy = np.zeros(3)
    eps = np.finfo(float).eps

    # Singular solution
    if (np.abs(R[0, 0]) < eps and np.abs(R[1, 0] < eps)):
        rpy[2] = 0.0
        rpy[1] = np.arctan2(R[2, 0], R[0, 0])
        rpy[0] = np.arctan2(R[1, 2], R[1, 1])

    # Non singular solution
    else:
        rpy[2] = np.arctan2(-R[1, 0], R[0, 0])
        c2 = np.cos(rpy[2])
        s2 = np.sin(rpy[2])
        rpy[1] = np.arctan2(R[2, 0], c2 * R[0, 0] - s2 * R[1, 0])
        rpy[0] = np.arctan2(-R[2, 1], R[2, 2])

    return rpy


def eul2dc(eul):
    """
    Returns the rotation matrix for a set of Z-X-Z Euler angles.
    """
    eul = np.asarray(eul).flatten()
    R = rotZ(eul[2]) @ rotX(eul[1]) @ rotZ(eul[0])

    return R


def dc2eul(R):
    """
    Returns the Z-X-Z Euler angles for the specified rotation matrix.

    Notes:
    - in the singular case (where there are an infinite number of solution) the
      solution with the third angle equal to zero is chosen.
    - in the general case (where there are always two solutions) the positive
      solution for the X angle is chosen.
    """
    eul = np.zeros(3)
    eps = np.finfo(float).eps

    # Singular solution
    if (np.abs(R[0, 2]) < eps and np.abs(R[1, 2] < eps)):
        eul[2] = 0.0
        eul[1] = np.arctan2(R[1, 2], R[2, 2])
        eul[0] = np.arctan2(R[0, 1], R[0, 0])

    # Non singular solution
    else:
        eul[2] = np.arctan2(R[0, 2], R[1, 2])
        c2 = np.cos(eul[2])
        s2 = np.sin(eul[2])
        eul[1] = np.arctan2(s2 * R[0, 2] + c2 * R[1, 2], R[2, 2])
        eul[0] = np.arctan2(R[2, 0], -R[2, 1])

    return eul


def tilde(a):
    """
    Returns the the skew-symmetric matrix of a vector.
    """
    B_skew = np.array([[  0.0, -a[2],  a[1]],
                       [ a[2],   0.0, -a[0]],
                       [-a[1],  a[0],   0.0]])

    return B_skew


def cross(u, v):
    """
    Returns the cross product between two vectors or a set vectors.
    """
    # Arguments are two vectors --> standard cross product
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


def rotW(w, dt):
    """
    Returns a rotation matrix representing the rotation |w| x dt about the
    angular velocity vector <w> using Rodrigues' formulation.
    """
    w = np.asarray(w).flatten()
    w_norm = np.linalg.norm(w)      # Norm
    K = w / w_norm                  # Unit vector defining the rotation axis
    theta = w_norm * dt             # Rotation
    K_tilde = tilde(K)
    st = np.sin(theta)
    ct = np.cos(theta)
    R = np.eye(3) + st * K_tilde + (1 - ct) * K_tilde @ K_tilde

    return R


def inertia_matrix(shape, *args):
    """
    Returns the moments of inertia wrt the centroid for basic shapes.

    Note: inertia for unit of mass !!!!!

    Explain meaning geometric quantities
    """
    inertia = np.zeros((3, 3))

    # Thick-walled cylindrical tube: Re, Ri, H
    # Set Ri = 0 for a solid cylindrical tube
    # case Re = Ri ?????
    if (shape == 'Cylinder'):
        Re = args[0]
        Ri = args[1]
        L = args[2]
        inertia[0, 0] = (Re ** 2 + Ri ** 2) / 2.0
        inertia[1, 1] = (3.0 * Re ** 2 + 3.0 * Ri ** 2 + L ** 2) / 12
        inertia[2, 2] = inertia[1, 1]

    # Thick-walled sphere: Re, Ri
    # Set Ri = 0 for a solid sphere
    # case Re = Ri ?????
    elif (shape == 'Sphere'):
        Re = args[0]
        Ri = args[1]
        inertia[0, 0] = 0.4 * (Re ** 5 - Ri ** 5) / (Re ** 3 - Ri ** 3)
        inertia[1, 1] = inertia[0, 0]
        inertia[2, 2] = inertia[0, 0]

    # Slender bar
    elif(shape == 'Bar'):
        L = args[0]
        inertia[0, 0] = 0.0
        inertia[1, 1] = L ** 2 / 12
        inertia[2, 2] = inertia[1, 1]

    # Rectangular prism
    # hollow case with thickness???
    elif(shape == 'Prism'):
        L = args[0]
        B = args[1]
        H = args[2]
        inertia[0, 0] = (B ** 2 + H ** 2) / 12
        inertia[1, 1] = (L ** 2 + H ** 2) / 12
        inertia[2, 2] = (L ** 2 + B ** 2) / 12

    # Default is the unit matrix
    else:
        inertia = np.eye(3)

    return inertia
