# The Spacedyn - A toolbox for space and mobile robots

This is a Python version of the [original Matlab code](http://www.astro.mech.tohoku.ac.jp/spacedyn/). Several new functions have also been added.

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
