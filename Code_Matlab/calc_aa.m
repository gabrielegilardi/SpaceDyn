%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1 // Oct 4, 1999, Last modification by K.Yoshida
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_AA 	Calculate the coordinate transfrom matrices in the Robotics convention.
%
%		CALC_TRN( A0 , q ) returns the coordinate tranformation matrices AA.
%		AA is a collection of A_I_1, A_I_2, ... A_I_n.
%

function AA = calc_aa( A0 , q )

global BB J_type Qi

num_q = length(q);

% If a Single body,
if num_q == 0
   
   AA = [];
   
% If a Multi body system,
else
  
   % Calculation of coordinate transformation matrices
   A_I_0 = A0;
   for i = 1 : num_q
      
      % Current (i-th) link connects to the 0-th link
      if BB(i) == 0
         % Rotational joint
         if J_type(i) == 'R'
            A_0_i = (rpy2dc( Qi(1,i) , Qi(2,i) , Qi(3,i)+q(i) )');
         % Prismatic joint
         else
            A_0_i = (rpy2dc( Qi(:,i) )');
         end
         AA(:,i*3-2:i*3) = A_I_0*A_0_i;
      % Current (i-th) link doesn't connect to the 0-th link
      else
         % Rotational joint
         if J_type(i) == 'R'
            A_BB_i = (rpy2dc( Qi(1,i), Qi(2,i), Qi(3,i)+q(i) )');
         % Prismatic joint
         else
            A_BB_i = (rpy2dc( Qi(:,i) )');
         end
         AA(:,i*3-2:i*3) = AA(:,BB(i)*3-2:BB(i)*3)*A_BB_i;
      end
   end
end


%%%EOF
