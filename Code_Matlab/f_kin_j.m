%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 11, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_KINE_J	Forward Kinematics:
%		Caluculate the Position/Orientation of the joints
%		corresponding the end point specified by connection vector 'joints.'
%		POS : 3x1 vector, ORI : 3x3 matrix,
%		POS_j : 3xn, ORI : 3x3n,
%		where n : number of joints between the End-point to link 0.
%
% 1998 1.9 (C)Space Robotics Lab
% A.Kurosu
% 1999.2.4 K.Yoshida
% 2001.9.11 H.Hamano
%

function [ POS_j , ORI_j ] = f_kin_j( RR, AA, q, joints )

global J_type
global cc
global Ez


% Check the number of the corresponding joints
n = length( joints );


% Calculation of Orientation and Position of each joints
POS_j = [];
ORI_j = [];

for i = 1 : 1 : n
   
   PorR = ( J_type(joints(i)) == 'P' );
   ORI_tmp = AA(:,joints(i)*3-2:joints(i)*3);
   POS_tmp = RR(:,joints(i)) + ORI_tmp*( cc(:,joints(i),joints(i)) - PorR*Ez*q(joints(i)) );
   
   POS_j = [POS_j POS_tmp];
   ORI_j = [ORI_j ORI_tmp];
   
end


%%%EOF
