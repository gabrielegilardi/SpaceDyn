%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Dec. 13, 2005, Last modification by H.Nakanishi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% INT_EU2	Integration of the system motion by a simple Euler method
%
%		INT_EU returns velocity & position vectors 
%		of each link of the system in the inertia frame.
%		It uses a simple, most primitive  Euler integration.
%
%		(The Rodorigues fomula for infinitesimal rotation
%		is used to update from A0 to A0_n.
%		This seems practically better than the fomula
%		delta_C0 = tilde(w0')*C0.)
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]=int_eu3(R0,A0,v0,w0,vd0,wd0,q,qd,qdd)
%	
%		NOTE:v0,w0,vd0,wd0 are in the inertia frame.
%
% 1999 10/4 K.Yoshida
%

function [ R0, A0, v0, w0, q, qd ] = int_eu2( R0, A0, v0, w0, vd0, wd0, q, qd, qdd )

global d_time

% Number of links
num_q = length( q );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integlation of q and qd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Single body
if num_q == 0
   q_n = [];
   qd_n = [];
   
% Multi bodies
else
   qd_n = qd + qdd * d_time;
   q_n  = q  + qd_n  * d_time;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of v0, R0, and A0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0_n = v0 + vd0 * d_time;
w0_n = w0 + wd0 * d_time;
R0_n = R0 + v0_n  * d_time;
A0_n = aw( w0_n ) * A0;


% outputs

R0 = R0_n;
A0 = A0_n;
v0 = v0_n;
w0 = w0_n;
q  = q_n;
qd = qd_n;

%%%EOF
