
import numpy as np
import elements_class as element
import model_class as mc
from utils import *

# Foot (base)
In = 1.5 * np.eye(3)
foot = element.base(2.3, In, cc={1: [0.05, 0.0, 0.0]}, name='foot')

# Leg
a = {1: [-0.2, 0.0, 0.0], 2: [ 0.2, 0.0, 0.0], 4: [ 0.2, 0.0, 0.0]}
leg = element.link(0.4, 0.5*In, cc=a, Qe=[0.1, 0.2, 0.3], name='leg')
j_leg = element.joint('R', Qi=[0.1, 0.2, 0.3], name='joint foot-leg')

# Trunk
a = {2: [-0.2, 0.0, 0.0], 3: [ 0.2, 0.0, 0.0]}
trunk = element.link(1.4, 0.3*In, cc=a, name='trunk')
j_trunk = element.joint('R', name='joint leg-trunk')

# Upper arm
a = {3: [-0.2, 0.0, 0.0]}
upper_arm = element.link(0.7, 2.3*In, cc=a, ce={3: [0.2, 0.0, 0.0]}, name='upper arm')
j_upper_arm = element.joint('R', name='joint trunk-upper arm')

# Lower arm
a = {4: [-0.2, 0.0, 0.0]}
lower_arm = element.link(0.8, 1.5*In, cc=a, ce={4: [0.2, 0.0, 0.0]}, name='lower arm')
j_lower_arm = element.joint('R', name='joint leg-lower arm')

# robot.set_param(Gravity=[-9.81, 0.0, 0.7], Ez=[0.0, 1.0, 1.0])
single = False
if (single):
    foot = element.base(2.3, 1.5*np.eye(3), name='foot')
    robot = mc.model([foot])
    robot.set_init(R0=[1.0, 2.0, 3.0], Q0=[0.1, 0.2, 0.3], v0=[1.0, 2.0, 3.0], w0=[0.1, 0.2, 0.3])

else:
    bodies = [foot, leg, trunk, upper_arm, lower_arm]
    joints = [j_leg, j_trunk, j_upper_arm, j_lower_arm]
    robot = mc.model(bodies, joints)
    robot.set_init(R0=[1.0, 2.0, 3.0], Q0=[0.1, 0.2, 0.3], v0=[1.0, 2.0, 3.0], \
        w0=[0.1, 0.2, 0.3], q=[0.1, 0.5, 1.0, 1.5], qd = [ 0.01, 0.02, 0.03, 0.04 ])

robot.update_elements()

R0=np.array([[ 1.7,  1.7], [-2.0, -2.0 ], [3.4,  3.4 ]])
Q0=np.array([[-0.1, -0.1], [ 0.23, 0.23], [0.76, 0.76]])
print(R0,type(R0))
print(Q0,type(Q0))
dd = R0 - np.reshape(Q0[:,0], (3,1))
print(dd,dd.ndim)


