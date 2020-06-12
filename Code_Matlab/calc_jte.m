%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 11, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JTE	Translational Jacobians (3 by nr. of joints in joints)
%           for the point specified by connection vector 'joints.'
%
%		1998.1.12 A.Kurosu
%       1998.2.4  K.Yoshida
%       2001.9.11 H.Hamano

function JJ_te = calc_jte( RR, AA, q, joints )

global J_type
global num_q Ez


% Check number of joint
n = length(joints);

% If a Single body,
if num_q == 0
   
   JJ_te = [];

% If a Multi body system,
else
   
   % Calculation of Joint Position
   [ POS_j, ORI_j ] = f_kin_j( RR, AA, q, joints );
   
   % Calculation of Effector Position
   [ POS_e, ORI_e ] = f_kin_e( RR, AA, joints );
   
   JJ_te = [];
   
   for i = 1 : 1 : n
      
      A_I_i = AA(:,joints(i)*3-2:joints(i)*3);
      
      % Rotational joint
      if J_type(joints(i)) == 'R'
         
         temp = cross( (A_I_i*Ez) , ( POS_e - POS_j(:,i) ) );
         JJ_te = [ JJ_te temp ];
         
      % Prismatic joint
      else
         
         JJ_te = [ JJ_te A_I_i*Ez ];
         
      end
      
   end
   
end


%%%EOF
