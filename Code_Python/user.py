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


def calc_forces(num_e):
    """
    Returns any external force/moment applied on the bodies.

    Fe = [Fe_0, Fe_1, ... ]         (3, num_e)
    Te = [Te_0, Te_1, ... ]         (3, num_e)

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
    Fe = np.zeros((3, num_e))           # Forces on the bodies
    Te = np.zeros((3, num_e))           # Moments on the bodies

    # Endpoint index 0
    Fe[:, 0] = np.array([-10.3, 11.4, 20.4])
    Te[:, 0] = np.array([2.2, -4.4, 1.6])

    # Endpoint index 0
    Fe[:, 1] = 1.2 * np.array([-10.3, 11.4, 20.4])
    Te[:, 1] = -0.7 * np.array([2.2, -4.4, 1.6])

    return Fe, Te
