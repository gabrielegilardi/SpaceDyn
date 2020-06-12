%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	March 9, 2005, Last modification by H.Nakanishi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_GH	Calculate the Generalized Inertia Matrices GH.
%
%		CALC_GH returns the inertia matrices GH (n)x(n).
%
%
% (C)Space Robotics Lab, by Hiroki Nakanishi
% March 9, 2005 made from calc_HH by Hiroki Nakanishi.
%

function GH = calc_gh( R0, RR, A0, AA )

global m0 m inertia0 inertia
global num_q

mass = m0 + sum(m);               % Total mass

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
           
   HH_m = [JJ_tg;
           HH_wq];
   HH_s = [          wE  tilde(wr0g)';
          tilde(wr0g)          HH_w];
           
  GH = HH_q - HH_m' * inv(HH_s) * HH_m;
end


%%% EOF
