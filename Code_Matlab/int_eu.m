%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% INT_EU	Integration of the system motion by a simple Euler method
%
%		INT_EU returns velocity & position vectors 
%		of each link of the system in the inertia frame.
%		It uses a simple, most primitive  Euler integration.
%
%		(The fomula of delta_CO = tilde(w0')*C0 is used.
%		But its integration is subject to error accumulation,
%		yielding non-normal, non-orthogonal C0.
%		int_eu2 is recommended for large attitude motion.) 
%
%		A constant time step is used.
%
%		[R0,A0,q,v0,w0,qd]=int_eu(R0,A0,q,v0,w0,qd,vd0,wd0,qdd)
%	
%		NOTE:v0,w0,vd0,wd0 are in the inertia frame.
%
% 1999 10/4 K.Yoshida
%


function [ R0, A0, v0, w0, q, qd ] = int_eu( R0, A0, v0, w0, vd0, wd0, q, qd, qdd )

global d_time

% Number of links
num_q = length( q );

% Single body
if num_q == 0
   q_n = [];
   qd_n = [];
   
% Multi bodies
else
   qd_n = qd + qdd * d_time;
   q_n  = q  + qd  * d_time;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of v0, R0, and A0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R0_n = R0 + v0  * d_time;
v0_n = v0 + vd0 * d_time;
w0_n = w0 + wd0 * d_time;

% Note that C0 is the direction cosines of body 0
% C0 = ( A0 )^T

C0 = ( A0 )';

% dC0 is the time derivative of C0
dC0 = (tilde(w0))' * C0;
C0_n = C0 + dC0 * d_time;

% outputs
R0 = R0_n;
A0 = ( C0_n )';
v0 = v0_n;
w0 = w0_n;
q  = q_n;
qd = qd_n;

%%%EOF
