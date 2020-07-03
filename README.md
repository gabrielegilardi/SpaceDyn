# SpaceDyn - A toolbox for space, mobile, and humanoid robotic systems

This is a Python version of the [original Matlab code](http://www.astro.mech.tohoku.ac.jp/spacedyn/) from the Space Robotics Lab. at Tohoku University (Sendai, Japan). New functions useful for control and for humanoid robots have also been added.

!!! Work in progress !!!

## Implemented functions

- File *elements.py*:
  - init (new function - class Link)
  - info (new function - class Link)
  - init (new function - class Joint)
  - info (new function - class Joint)

- File *model.py*:
  - build_cc_SS (connectivity between bodies/joints)
  - build_ce_Qe_SE (connectivity between bodies/end-points)
  - build_mass_inertia (system mass and inertia)
  - build_BB (branch link sequence)
  - build_j_type (joint type)
  - build_Qi (joint frames)
  - init (new function - class Model)
  - info (new function - class Model)
  - set_param (new function - initialize system param)
  - set_init (new function - initialize system)
  - calc_CoM (new function - position, velocity, and acceleration of the CoM)

- File *kinematics.py*:
  - j_num (link sequence to an end-point)
  - f_kin_e (position and orientation of an end-point)
  - f_kin_j (position and orientation of all joints in a link sequence)
  - calc_jte (translational Jacobian of an end-point)
  - calc_jre (rotational Jacobian of an end-point)
  - calc_je (Jacobian of an end-point)
  - calc_aa (link rotation matrices)
  - calc_pos (link centroid positions)
  - calc_vel (link centroid velocities)
  - calc_acc (link centroid accelerations)
  - calc_jt (translational jacobian wrt link centroids)
  - calc_jr (rotational jacobian wrt link centroids)
  - calc_hh (system inertia matrix)

- File *dynamics.py*:
  - calc_Forces (user-defined external forces)
  - r_ne (inverse dynamics by recursive Newton-Euler method)
  - f_dyn_nb2 (integration using Newmark-beta method)

- File *utils.py*:
  - rotX
  - rotY
  - rotZ
  - tilde
  - rpy2dc
  - dc2rpy
  - eul2dc
  - dc2eul
  - cross
  - tr2diff
  - rotW
  - inertia_matrix

- Functions done in Matlab but not ported in Python yet:
  - f_dyn (forward dynamics)
  - f_dyn_rk2 (integration using Runge-Kutta)
  - i_dyn (inverse dynamics)
  - i_kine (inverse kinematics)
  - calc_gh (generalized inertia matrix)
  - calc_gj (generalized Jacobian)
  - calc_Lin_Mom (new function - system linear momentum)
  - calc_Ang_Mom (new function - system angular momentum)
  - calc_der_Lin_Mom (new function - system linear momentum 1st derivative)
  - calc_der_Ang_Mom (new function- system angular momentum 1st derivative)
  - calc_Kin_Ene (new function- system kinetic energy)
  - calc_Pot_Ene (new function - system potential energy)
  - calc_Work (new function - system work)
  - inv_Kin_Func (new function - inverse kinematics for system CoM)
  - Contact_Forces (new function - contact forces)
  - PSO (new function - particle swarm optimizer)
  - Animation (new function - animate system)
  - SaveResults (new function - save results)
