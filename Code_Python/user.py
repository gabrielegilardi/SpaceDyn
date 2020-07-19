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


def calc_Forces(num_j=0):
    """
    Returns the external forces/moments/torques applied to base, links, and
    joints.
    """
    num_b = num_j + 1
    Fe = np.zeros((3, num_b))           # Forces on base and link end-points
    Te = np.zeros((3, num_b))           # Moments on base and link end-points
    tau = np.zeros(num_j)               # Torques/forces on joints

    return Fe, Te, tau
