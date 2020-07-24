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
from utils import cross, rotW


def calc_forces(time=0.0, num_j=0, num_e=0):
    """
    Returns any external/control force/moment applied to the bodies.

    F0          (3, )           # Control force on the base centroid
    T0          (3, )           # Control moment on the base centroid
    tau         (num_j, )       # Control torques/forces on the joints

    Fe = [Fe_0, Fe_1, ... ]     (3, num_e)      # Endpoint external forces
    Te = [Te_0, Te_1, ... ]     (3, num_e)      # Endpoint external moments

    If the endpoints are defined by:

    ee = {
      0: (link_id, position, orientation),
      1: (link_id, position, orientation),
      .....
      num_e: (link_id, position, orientation)
         }

    then Fe_0 and Te_0 are applied to the endpoint with index 0, Fe_1 and Te_1
    are applied to the endpoint with index 1, etc., up to endpoint with index
    num_e 
    """
    # Control terms
    F0 = np.array([-1.7, 2.4, -4.5])
    T0 = np.array([0.3, -0.2, 0.13])
    tau = np.array([0.1, -0.3, 0.6, -1.1])

    # Endpoint terms
    Fe = np.zeros((3, num_e))
    Te = np.zeros((3, num_e))

    # Endpoint index 0
    Fe[:, 0] = np.array([-10.3, 11.4, 20.4])
    Te[:, 0] = np.array([2.2, -4.4, 1.6])

    # Endpoint index 01
    Fe[:, 1] = 1.2 * np.array([-10.3, 11.4, 20.4])
    Te[:, 1] = -0.7 * np.array([2.2, -4.4, 1.6])

    return F0, T0, tau, Fe, Te

