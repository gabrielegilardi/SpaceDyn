%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JRE	Translational Jacobians (3 by nr. of joints in joints)
%           for the end point specified by connection vector 'joints'
%
%		1998.1.12 A.Kurosu
%       1998.2.4  K.Yoshida
%		2001.9.11 H.Hamano
%

function JJ_re = calc_jre( AA, joints )

global J_type
global num_q Ez


% Check number of joints
n = length(joints);


% If a Single body,
if num_q == 0
   
   JJ_re = [];
   
% If a Multi body system,
else
   
   JJ_re = [];
   
   for i = 1 : 1 : n
      
      A_I_i = AA(:,joints(i)*3-2:joints(i)*3);
      
      % Rotational joint
      if J_type(joints(i)) == 'R'
         JJ_re = [ JJ_re A_I_i*Ez ];
         
      % Prismatic joint
      else
         JJ_re = [ JJ_re [ 0 0 0 ]' ];
         
      end
      
   end
   
end


%%%EOF
