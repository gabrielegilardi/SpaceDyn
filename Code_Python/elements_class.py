
# Libraries and external files

import numpy as np


# Base class
class base:

    def __init__(self, mass=0.0, inertia=np.eye(3), cc={}, name=None):

        self.mass = mass                        # Mass [kg]
        self.inertia = np.asarray(inertia)      # Moments of inertia [kgm^2]
        self.cc = cc
        self.name = name

        self.AA = np.zeros((3,3))
        self.RR = np.zeros(3)
        self.vv = np.zeros(3)
        self.ww = np.zeros(3)



# Link class
class link(base):

    def __init__(self, mass=0.0, inertia=np.eye(3), cc={}, ce={}, Qe=np.zeros(3), name=None):

        super().__init__(mass, inertia, cc, name)
        self.ce = ce
        self.Qe = np.asarray(Qe)                # End-point (relative) orientation


# Joint class
class joint:

    def __init__(self, j_type='R', Qi=np.zeros(3), name=None):
        
        self.j_type = j_type                    # Joint type (R or P)
        self.Qi = np.asarray(Qi)                # Joint (relative) orientation
        self.name = name

        self.q = np.zeros(3)
        self.qd = np.zeros(3)


