%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_DYN_RK	Compute the Forward Dynamics of the system
%		with Integration by Runge-Kutta method
%		(The fomula of delta_CO = tilde(w0')*C0 is used.
%		But its integration is subject to error accumulation,
%		yielding non-normal, non-orthogonal C0.
%		f_dyn_rk2 is recommended for large attitude motion.) 
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]
%		   = f_dyn_rk(R0,A0,q,v0,w0,q,qd,F0,T0,Fe,Te,tau)
%
%		NOTE: v0,w0,vd0,wd0 are in the inertia frame.
%
% 1999 10/4 K.Yoshida
%

function [ R0, A0, v0, w0, q, qd ] = f_dyn_rk( R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau )

global d_time

% Definition

C0 = A0';

% 1st Step

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);

k1_R0 = d_time * v0;
k1_C0 = d_time * tilde(w0)'*C0;
k1_A0 = k1_C0';
k1_q  = d_time * qd;

k1_v0 = d_time * tmp_vd0;
k1_w0 = d_time * tmp_wd0;
k1_qd = d_time * tmp_qdd;

% 2nd Step

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0+k1_R0/2,A0+k1_A0/2,v0+k1_v0/2,w0+k1_w0/2,q+k1_q/2,qd+k1_qd/2,F0,T0,Fe,Te,tau);

k2_R0 = d_time * ( v0 + k1_v0/2 );
k2_C0 = d_time * ( tilde( w0 + k1_w0/2 )' * C0 );
k2_A0 = k2_C0';
k2_q  = d_time * ( qd + k1_qd/2 );

k2_v0 = d_time * tmp_vd0;
k2_w0 = d_time * tmp_wd0;
k2_qd = d_time * tmp_qdd;

% 3rd Step

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0+k2_R0/2,A0+k2_A0/2,v0+k2_v0/2,w0+k2_w0/2,q+k2_q/2,qd+k2_qd/2,F0,T0,Fe,Te,tau);
k3_R0 = d_time * ( v0 + k2_v0/2 );
k3_C0 = d_time * ( tilde( w0 + k2_w0/2 )'* C0 );
k3_A0 = k3_C0';
k3_q  = d_time * ( qd + k2_qd/2 );

k3_v0 = d_time * tmp_vd0;
k3_w0 = d_time * tmp_wd0;
k3_qd = d_time * tmp_qdd;

% 4th Step

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0+k3_R0,A0+k3_A0,v0+k3_v0,w0+k3_w0,q+k3_q,qd+k3_qd,F0,T0,Fe,Te,tau);
k4_R0 = d_time * ( v0 + k3_v0 );
k4_C0 = d_time * ( tilde( w0 + k3_w0 )'* C0 );
k4_A0 = k4_C0';
k4_q  = d_time * ( qd + k3_qd );

k4_v0 = d_time * tmp_vd0;
k4_w0 = d_time * tmp_wd0;
k4_qd = d_time * tmp_qdd;

% Compute Values at the Next Time Step

R0_next = R0 + ( k1_R0 + 2*k2_R0 + 2*k3_R0 + k4_R0 )/6;
C0_next = C0 + ( k1_C0 + 2*k2_C0 + 2*k3_C0 + k4_C0 )/6;
A0_next = C0_next';
q_next  = q + ( k1_q + 2*k2_q + 2*k3_q + k4_q )/6;
v0_next = v0 + ( k1_v0 + 2*k2_v0 + 2*k3_v0 + k4_v0 )/6;
w0_next = w0 + ( k1_w0 + 2*k2_w0 + 2*k3_w0 + k4_w0 )/6;
qd_next = qd + ( k1_qd + 2*k2_qd + 2*k3_qd + k4_qd )/6;

% Solution

R0 = R0_next;
A0 = A0_next;
q  = q_next;
v0 = v0_next;
w0 = w0_next;
qd = qd_next;


%%% EOF
