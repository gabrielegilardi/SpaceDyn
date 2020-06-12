%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 17, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_VEL	Calculate the velocity of each link
%
%		CALC_VEL returns the link velocity in the inertial frame
%		for link 1 to n.
%
%	'97 6/12 T.Hiraoka
%   '99 2/4  K.Yoshida
%  2001.9.17 H.Hamano

function [ vv,ww ] = calc_vel( A0, AA, v0, w0, q, qd )

global BB J_type
global c0 cc
global num_q Ez

% If a Single body,
if num_q == 0
   vv = [];
   ww = [];

% If a Multi body system,
else
   
   % Calculation of velocity vectors vv,ww
   for i = 1 : num_q
      
      % Check the link connection: Is the lower one of this link, 0 ?
      if BB(i) == 0
         
         % Current (i-th) link connects to the 0-th link
         A_I_i = AA(:,i*3-2:i*3);
         
         % Rotational joint
         if J_type(i) == 'R'
            ww(:,i) = w0(:) + A_I_i*Ez*qd(i);
            vv(:,i) = v0(:) ...
               + cross( w0(:),(A0*c0(:,i)) ) ...
               - cross( ww(:,i),(A_I_i*cc(:,i,i)) );
            
         % Prismatic joint
         else
            ww(:,i) = w0(:);
            vv(:,i) = v0(:) ...
               + cross( w0(:),(A0*c0(:,i)) ) ...
               - cross( ww(:,i),(A_I_i*cc(:,i,i)) ) ...
               + cross( ww(:,i),(A_I_i*Ez*q(i)) ) ...
               + A_I_i*Ez*qd(i);
            
         end
         
      % Current (i-th) link doesn't connect to the 0-th link
      else
         
         A_I_BB = AA(:,BB(i)*3-2:BB(i)*3);
         A_I_i  = AA(:,i*3-2:i*3);
         
         % Rotational joint
         if J_type(i) == 'R'
            ww(:,i) = ww(:,BB(i)) + A_I_i*Ez*qd(i);
            vv(:,i) = vv(:,BB(i)) ...
               + cross( ww(:,BB(i)),(A_I_BB*cc(:,BB(i),i)) ) ...
               - cross( ww(:,i),(A_I_i*cc(:,i,i)) );
            
         % Prismatic joint
         else
            ww(:,i) = ww(:,BB(i));
            vv(:,i) = vv(:,BB(i)) ...
               + cross( ww(:,BB(i)),(A_I_BB*cc(:,BB(i),i)) ) ...
               - cross( ww(:,i),(A_I_i*cc(:,i,i)) ) ...
               + cross( ww(:,i),(A_I_i*Ez*q(i)) ) ...
               + A_I_i*Ez*qd(i);
            
         end
         
      end
      
   end
   
end


%%%EOF
