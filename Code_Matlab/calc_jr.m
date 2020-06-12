%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 12, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CALC_JR	Rotational Jacobians w.r.t. link centroid
%
%		'97 6/16 T.hiraoka
%       '99 2/4  K.Yoshida
%      2001.9.12 H.Hamano
%

function JJ_r = calc_jr( AA )

global BB J_type
global num_q Ez


% Calculation of translational jacobians
JJ_r = zeros(3,num_q*num_q);

% If a Single body,
if num_q == 0
   
   JJ_r = [];

% If a Multi body system,
else
   
   for i = 1 : num_q
      
      A_I_i = AA(:,i*3-2:i*3);
         
      % Rotational joint
      if J_type(i) == 'R'
         JJ_r(:,(i-1)*num_q+i) = A_I_i*Ez;
         
      % Prismatic joint
      else
         JJ_r(:,(i-1)*num_q+i) = [ 0 0 0 ]' ;
         
      end
      
      j = BB(i);
      
      if j ~= 0
         JJ_r(:,(i-1)*num_q+1:(i-1)*num_q+i-1) = JJ_r(:,(j-1)*num_q+1:(j-1)*num_q+i-1);
         
      end
      
   end
   
end


%%%EOF
