%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% F_DYN_NB	Compute the Forward Dynamics of the system
%		with Integration by Newmark beta method.
%		Newmark beta is an implicit integration method,
%		considered robust to the energy accumulation.
%
%		(The Rodorigues fomula for infinitesimal rotation
%		is used to update from A0 to A0_n.
%		This seems practically better than the fomula
%		delta_C0 = tilde(w0')*C0.)
%
%		A constant time step is used.
%
%		[R0,A0,v0,w0,q,qd]
%			= f_dyn_nb2(R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau)
%
%		NOTE: v0,w0,vd0,wd0 are in the inertia frame.
%
%
% 1998.12.22  K.HASHIZUME
% 1999.10.4   K.Yoshida

function [ R0, A0, v0, w0, q, qd ] = f_dyn_nb2( R0, A0, v0, w0, q, qd, F0, T0, Fe, Te, tau )

global d_time

% Number of links
num_q = length( q );

% Definition

B  = 1/6 ;    % beta coefficient
K1 = d_time ;               %
K2 = d_time * d_time / 3 ;  %   Definition the coefficient of
K3 = d_time * d_time * B ;  %   the Newmark beta formula
K4 = d_time / 2 ;           %

% 1st Step

[ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);

% 1st leading R and v

Rdd0       = tmp_vd0 ;
Rdd0_next1 = Rdd0 ;       % Rdd(t+d_t)=Rdd(t)
R0_next1   = R0 + K1*v0 + K2*Rdd0 + K3*Rdd0_next1 ; % v0=Rd0
Rd0_next1  = v0 + K4*Rdd0 + K4*Rdd0_next1 ;
v0_next1   = Rd0_next1 ;  % v0=Rd0

% 1st leading q

qdd       = tmp_qdd ;
qdd_next1 = qdd ;         % qdd(t+d_t)=qdd(t)
q_next1   = q + K1*qd + K2*qdd + K3*qdd_next1 ;
qd_next1  = qd + K4*qdd + K4*qdd_next1 ;

% 1st leading A0 and w

wd0       = tmp_wd0 ;
wd0_next1 = wd0 ;
w0_next1  = w0 + K4*wd0 + K4*wd0_next1 ;
A0_next1  = aw(w0_next1) * A0;

% Definition the number of Repetition

for i = 1 : 1 : 1 ;  % number of repetiton (n+1)    ( In this case n=1 )
   
   % 2nd Step
   
   [ tmp_vd0,tmp_wd0,tmp_qdd ] = f_dyn( R0_next1,A0_next1,v0_next1,w0_next1,q_next1,qd_next1,F0,T0,Fe,Te,tau) ;
   
   % 2nd leading R and v
   
   Rdd0      = tmp_vd0 ;
   Rdd0_next = Rdd0 ;
   R0_next   = R0 + K1*v0 + K2*Rdd0 + K3*Rdd0_next ;
   Rd0_next  = v0 + K4*Rdd0 + K4*Rdd0_next ;
   v0_next   = Rd0_next ;
   
   % 2nd leading q
   
   qdd      = tmp_qdd ;
   qdd_next = qdd ;
   q_next   = q + K1*qd + K2*qdd + K3*qdd_next ;
   qd_next  = qd + K4*qdd + K4*qdd_next ;
   
   % 2nd leading C and w
   
   wd0      = tmp_wd0 ;
   wd0_next = wd0 ;
   w0_next  = w0 + K4*wd0 + K4*wd0_next ;
   A0_next  = aw(w0_next1) * A0;
   
   % Ready for next step
   
   R0_next1 = R0_next ;
   A0_next1 = A0_next ;
   v0_next1 = v0_next ;
   w0_next1 = w0_next ;
   q_next1  = q_next ;
   qd_next1 = qd_next ;
   
end

% Solution

R0 = R0_next1 ;
A0 = A0_next1 ;
v0 = v0_next1 ;
w0 = w0_next1 ;
q  = q_next1 ;
qd = qd_next1 ;


%%% EOF
