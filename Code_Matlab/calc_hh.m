%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 12, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_HH	Calculate the Inertia Matrices H.
%
%		CALC_HH returns the inertia matrices HH (6+n)x(6+n).
%
%
% (C)Space Robotics Lab, by Koichi Fujishima
% 1999.2.4  modified by K.Yoshida
% 2001.9.12 modified by H.Hamano
%

function HH = calc_hh( R0, RR, A0, AA )

global m0 m inertia0 inertia
global num_q

mass = m0 + sum(m);

% Calculation of partial translational & rotational jacobian
JJ_t = calc_jt( RR, AA );
JJ_r = calc_jr( AA );


% Calculation of HH matrix
wE = mass * eye(3,3);

JJ_tg = zeros(3,num_q);
HH_w = zeros(3,3);
HH_wq = zeros(3,num_q);
HH_q  = zeros(num_q,num_q);


% If a Single body,
if num_q == 0
   HH_w = A0*inertia0*A0';
   HH = [         wE  zeros(3,3);
          zeros(3,3)        HH_w ];
   
   
% Multi body system
else
   
   % Calculation of the position of gravity centroid, Rg
   Rm = zeros(3,1);
   
   for i = 1 : num_q
      
      Rm = Rm + m(i) * RR(:,i);
      
   end
   
   Rm = Rm + m0 * R0;
   Rg = Rm / mass;
   
   wr0g = (Rg-R0) * mass;

   
   for i = 1 : num_q
      
      r0i = RR(:,i) - R0;
      A_I_i = AA(:,i*3-2:i*3);
      JJ_tg = JJ_tg + m(i)*JJ_t(:,(i-1)*num_q+1:i*num_q);
      
      HH_w  = HH_w + A_I_i*inertia(:,i*3-2:i*3)*A_I_i' ...
         + m(i)*(tilde(r0i))'*(tilde(r0i));
      HH_wq = HH_wq ...
         + (A_I_i*inertia(:,i*3-2:i*3)*A_I_i')*JJ_r(:,(i-1)*num_q+1:i*num_q) ...
         + m(i) * (tilde(r0i)) * JJ_t(:,(i-1)*num_q+1:i*num_q);
      HH_q  = HH_q ...
         + JJ_r(:,(i-1)*num_q+1:i*num_q)'* ...
         (A_I_i*inertia(:,i*3-2:i*3)*A_I_i')* ...
         JJ_r(:,(i-1)*num_q+1:i*num_q) ...
         + m(i) * JJ_t(:,(i-1)*num_q+1:i*num_q)'*JJ_t(:,(i-1)*num_q+1:i*num_q);
      
   end
   
   HH_w = HH_w + A0*inertia0*A0';
   HH = [          wE  tilde(wr0g)'  JJ_tg;
          tilde(wr0g)          HH_w  HH_wq;
               JJ_tg'        HH_wq'   HH_q ];
   
end


%%% EOF
