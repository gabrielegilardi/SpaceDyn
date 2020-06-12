%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    version 1.2 // Feb 5, 2006, Last modification by H.Nakanishi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_DYN_TR	Compute the Forward Dynamics of the system
%		with Integration by Trapezoidal Rule (2nd order Adams-Moulton) method
%
%		A constant time step is used.
%
%		[ R0, A0, v0, w0, vd0, wd0, q, qd, qdd ] = 
%                           f_dyn_tr( R0, A0, v0, w0, vd0, wd0, q, qd, qdd, F0, T0, Fe, Te, tau )
%
%		NOTE: v0,w0,vd0,wd0 are in the inertia frame.
%

function [ R0, A0, v0, w0, vd0, wd0, q, qd, qdd ] = f_dyn_tr( R0, A0, v0, w0, vd0, wd0, q, qd, qdd, F0, T0, Fe, Te, tau )

global d_time

% Number of links
num_q = length( q );

% Definition

vd0_beforestep = vd0;
wd0_beforestep = wd0;
qdd_beforestep = qdd;

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integlation of q and qd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Single body
if num_q == 0
   q_n = [];
   qd_n = [];
   
% Multi bodies
else
   qd_n = qd + (qdd_beforestep + tmp_qdd) * d_time / 2;
   q_n  = q  + (qd + qd_n) * d_time / 2;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration of v0, R0, and A0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v0_n = v0 + (vd0_beforestep + tmp_vd0) * d_time / 2;
w0_n = w0 + (wd0_beforestep + tmp_wd0) * d_time / 2;
R0_n = R0 + (v0 + v0_n) * d_time / 2;
A0_n = (aw(w0) + aw(w0_n)) * A0 / 2;


% outputs

R0 = R0_n;
A0 = A0_n;
v0 = v0_n;
w0 = w0_n;
vd0 = tmp_vd0;
wd0 = tmp_wd0;
q  = q_n;
qd = qd_n;
qdd = tmp_qdd;

%%%EOF
