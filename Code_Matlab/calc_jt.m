%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 12, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_JT	Translational Jacobians w.r.t. link centroid
%
%		'97 6/16 T.hiraoka
%       '99 2/4  K.Yoshida
%      2001.9.12 H.Hamano
%

function JJ_t = calc_jt( RR, AA )

global BB J_type
global cc
global num_q Ez


% Calculation of translational jacobians
JJ_t = zeros(3,num_q*num_q);

% If a Single body,
if num_q == 0
   
   JJ_t = [];

% If a Multi body system,
else
   
   for i = 1 : num_q
      
      j = i;
      
      while ( j > 0 )
         
         A_I_j  = AA(:,j*3-2:j*3);
         
         % Rotational joint
         if J_type(j) == 'R'
            JJ_t(:,(i-1)*num_q+j) = cross( (A_I_j*Ez) , ( RR(:,i)-RR(:,j)-A_I_j*cc(:,j,j) ) );
            
         % Prismatic joint
         else
            JJ_t(:,(i-1)*num_q+j) = A_I_j*Ez ;
            
         end
         
         j = BB(j);
         
      end
      
   end
   
end


%%%EOF
