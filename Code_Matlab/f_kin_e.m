%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_KIN_E	Foward Kinematics
%		Caluculate the Position/Orientation of the end point
%		specified by connection vector 'joints.'
%		POS : 3x1 vector, ORI : 3x3 matrix.
%
% 1998.1.12 (C)Space Robotics Lab
% A.Kurosu
% 1999.2.4 K.Yoshida
% 2001.9.11 H.Hamano
%

function [ POS_e , ORI_e ] = f_kin_e( RR, AA, joints )

global Qe
global ce

% Check number of the corresponding joints
n = length(joints);
k = joints(n);

% Calcurate coordinate trasformation matrix of Effector
A_I_i = AA(:,k*3-2:k*3);
A_i_EE = rpy2dc(Qe(:,k))';
ORI_e = A_I_i*A_i_EE;

% Calculate position vector of Effector
POS_e = RR(:,k) + A_I_i*ce(:,k);


%%%EOF

