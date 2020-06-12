%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 17, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_POS	Calculate the position vectors of each link
%
%		CALC_POS returns the position vectors RR in the Inertia frame.
%
%		'97 6/12 T.Hiraoka
%       '99 2/4  K.Yoshida
%      2001.9.17 H.Hamano
%

function RR = calc_pos( R0, A0, AA, q )

global BB J_type
global c0 cc
global num_q Ez

% If a Single body,
if num_q == 0
   RR = [];
   
% If a Multi body system,
else
   % Calculation of position vectors
   for i = 1 : num_q
      A_I_i = AA(:,i*3-2:i*3);      
      % Current (i-th) link connects to the 0-th link
      if BB(i) == 0         
         % Rotational joint
         if J_type(i) == 'R'
            RR(:,i) = R0(:) + A0*c0(:,i) - A_I_i*cc(:,i,i);            
         % Prismatic joint
         else
            RR(:,i) = R0(:) + A0*c0(:,i) + A_I_i*( Ez*q(i)-cc(:,i,i) );            
         end
         
      % Current (i-th) link doesn't connect to the 0-th link
      else         
         A_I_BB = AA(:,BB(i)*3-2:BB(i)*3);        
         % Rotational joint
         if J_type(i) == 'R'
            RR(:,i) = RR(:,BB(i)) + A_I_BB*cc(:,BB(i),i) - A_I_i*cc(:,i,i);
         % Prismatic joint
         else
            RR(:,i) = RR(:,BB(i)) + A_I_BB*cc(:,BB(i),i) + A_I_i*( Ez*q(i)-cc(:,i,i) );            
         end         
      end      
   end   
end


%%%EOF
