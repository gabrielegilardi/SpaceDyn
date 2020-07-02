# The Spacedyn - A toolbox for space and mobile robots

This is a Python version of the [original Matlab code](http://www.astro.mech.tohoku.ac.jp/spacedyn/) from the Space Robotics Lab. (Prof. Kazuya Yoshida) at Tohoku University, Sendai, Japan. Several new functions have been added.

!!! Work in progress !!!

## Implemented functions

- File *elements.py*:
  - init (class Link)
  - info (class Link)
  - init (class Joint)
  - info (class Joint)

- File *model.py*:
  - build_cc_SS
  - build_ce_Qe_SE
  - build_mass_inertia
  - build_BB
  - build_j_type
  - build_Qi
  - init (class Model)
  - info (class Model)
  - set_param
  - set_init

- File *kinematics.py*:
  - j_num
  - f_kin_e
  - f_kin_j
  - calc_jte
  - calc_jre
  - calc_je
  - calc_aa
  - calc_pos
  - calc_vel
  - calc_acc
  - calc_jt
  - calc_jr
  - calc_hh

- File *dynamics.py*:
  - calc_Forces
  - r_ne
  - f_dyn_nb2

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
  - f_dyn
  - f_dyn_rk2
  - i_dyn
  - i_kine
  - calc_gh
  - calc_gj
  - cacl_pos_CoM (new function)
  - calc_vel_CoM (new function)
  - calc_acc_CoM (new function)
  - calc_Lin_Mom (new function)
  - calc_Ang_Mom (new function)
  - calc_der_Lin_Mom (new function)
  - calc_der_Ang_Mom (new function)
  - calc_Kin_Ene (new function)
  - calc_Pot_Ene (new function)
  - calc_Work (new function)
  - inv_Kin_Func (new function)
  - contact_model (new function)
  - PSO (new function)
  - Animation (new function)
