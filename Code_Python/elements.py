"""
SpaceDyn - A Toolbox for Space and Mobile Robots.

Copyright (c) 2020 Gabriele Gilardi

Python version of:
    The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
    (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
    Tohoku University, Japan.
"""

import numpy as np


# Link class
class link:

    def __init__(self, name=None, mass=0.0, inertia=np.eye(3), cc={}):
        """
        name        str         Name
        mass        scalar      Mass [kg]
        inertia     (3, 3)      Inertia matrix [kgm^2]
        cc          dict        Identifiers/positions of connected bodies [m]
        """
        self.name = name
        self.mass = mass
        self.inertia = np.asarray(inertia)
        self.cc = cc

        self.AA = np.array([])          # Rotation matrix
        self.RR = np.array([])          # Position [m]
        self.vv = np.array([])          # Velocity [m/s]
        self.ww = np.array([])          # Acceleration [m/s^2]

    def info(self):
        """
        """
        print()
        print('Type: LINK')
        print('Name = ', self.name)
        print('Mass = ', self.mass)
        print('Inertia =')
        print(self.inertia)
        print('Connections =')
        print(self.cc)
        print('Rotation matrix =')
        print(self.AA)
        print('Position = ', self.RR)
        print('Velocity = ', self.vv)
        print('Acceleration = ', self.ww)
        print()

        return


# Joint class
class joint:

    def __init__(self, name=None, j_type='R', Qi=np.zeros(3)):
        """
        name        str         Name
        j_type      str         Joint type [R or P]
        Qi          (3, )       Joint (link) relative orientation [rad]

        Notes:
        - the joint and link frame are always parallel.
        - the frame is given with respect to the previous joint/link frame.
        """
        self.name = name
        self.j_type = j_type
        self.Qi = np.asarray(Qi)

        self.q = np.array([])           # Rotation [rad]
        self.qd = np.array([])          # Angular velocity [rad/s]
        self.qdd = np.array([])         # Angular acceleration [rad/s^2]

    def info(self):
        """
        """
        print()
        print('Element type: JOINT')
        print('Name = ', self.name)
        print('Type = ', self.j_type)
        print('Orientation = ', self.Qi)
        print('Rotation = ', self.q)
        print('Angular velocity = ', self.qd)
        print('Angular acceleration = ', self.qdd)
        print()

        return
