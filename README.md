# SpaceDyn - A toolbox for Space, Mobile, and Humanoid Robotic Systems

This is a Python version of the [original Matlab code](http://www.astro.mech.tohoku.ac.jp/spacedyn/) from the Space Robotics Lab. at Tohoku University (Sendai, Japan). New functions useful for control and for humanoid robots have also been added.

!!! Work in progress !!!

## Implemented functions

- File *elements.py*:
  - init (new - class Base, tested)
  - info (new - class Base, tested)
  - init (new - class Link, tested)
  - info (new - class Link, tested)
  - init (new - class Model, tested)
  - info (new - class Model, tested)
  - connectivity (new - connectivity between bodies and endpoints, tested)
  - build_cc_Qi (new - position/orientation links, tested)
  - build_ce_Qe (new - position/orientation end-points, tested)
  - build_mass_inertia (new - system mass and inertia, tested)
  - set_init (new - initialize system)
  - calc_CoM (new - position, velocity, and acceleration of the CoM)
  - calc_kin_ener (new - kinetic energy for system/links)
  - calc_pot_ener (new - potential energy for system/links)
  - calc_lin_mom (new - linear momentum for system/links)
  - calc_lin_mom1 (new - system linear momentum 1st derivative)
  - calc_ang_mom (new - angular momentum for system/links)
  - calc_ang_mom1 (new - system angular momentum 1st derivative)

- File *kinematics.py*:
  - j_num (link sequence to an end-point, tested)
  - f_kin_e (position and orientation of an end-point, tested)
  - f_kin_j (position and orientation of all joints in a link sequence, tested)
  - calc_je (Jacobian of an end-point, tested)
  - calc_aa (link rotation matrices, tested)
  - calc_pos (link centroid positions, tested)
  - calc_vel (link centroid velocities, tested)
  - calc_acc (link centroid accelerations)
  - calc_jt (translational jacobian wrt link centroids)
  - calc_jr (rotational jacobian wrt link centroids)
  - calc_hh (system inertia matrix)

- File *dynamics.py*:
  - calc_Forces (new - user-defined external forces)
  - r_ne (inverse dynamics by recursive Newton-Euler method)
  - f_dyn_nb2 (integration using Newmark-beta method)
  - f_dyn_rk2 (integration using Runge-Kutta)
  - f_dyn (forward dynamics)

- File *user.py*:
  - calc_Forces (new - user-defined external forces)

- File *utils.py*:
  - rotX (rotation about X axis, tested)
  - rotY (rotation about Y axis, tested)
  - rotZ (rotation about Z axis, tested)
  - tilde (skew-symmetric forn of a vector, tested)
  - rpy2dc (RPY angles to rotation matrix, tested)
  - dc2rpy (rotation matrix to RPY angles, tested)
  - eul2dc (Euler angles to rotation matrix, tested)
  - dc2eul (rotation matrix to Euler angles, tested)
  - cross (single/multi-vector cross product, tested)
  - rotW (rotation vector to rotation matrix, tested)
  - inertia (new - moments of inertia for basic shapes, tested)

- Functions done in Matlab but not ported in Python yet:
  - i_dyn (inverse dynamics)
  - i_kine (inverse kinematics)
  - calc_gh (generalized inertia matrix)
  - calc_gj (generalized Jacobian)
  - calc_Work (new - work for system/links)
  - inv_Kin_Func (new - inverse kinematics for system CoM)
  - Contact_Forces (new - contact forces)
  - PSO (new - particle swarm optimizer)
  - Animation (new - animate system)
  - SaveResults (new - save results)
  - Joint_Limits (new - allows to impose limits on the joints)
