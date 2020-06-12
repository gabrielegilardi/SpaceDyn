%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% CALC_GJ	Calculate the Generalized Jacobian.
%
%		CALC_GJ returns the generalized jacobian, GJ (6xn).
%
% 1998 (C)Space Robotics Lab, by Koichi Fujishima
% 2001.9.13 H.Hamano
%

function GJ = calc_gj( R0, RR, A0, AA, q, num_e )

num_q = length(q);

% Calculate inertia matrices, HH
HH = calc_hh( R0, RR, A0, AA );

% Find joint connection from the end-link to the 0-th link
joints = j_num( num_e );

% calculate Jacobian and inertia matrices
Jm = calc_je( RR, AA, q, joints );

[ pe, tmp1 ] = f_kin_e( RR, AA, joints );
Js = [   eye(3,3)  -tilde(pe-R0);
       zeros(3,3)       eye(3,3) ];

Hs = HH(1:6,1:6);
Hm = HH(1:6,7:6+num_q);

% Calculate the Generalized Jacobian
GJ = Jm - Js*inv(Hs)*Hm;

%%% EOF
