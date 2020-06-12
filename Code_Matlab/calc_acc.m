%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 13, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% CALC_ACC 	Calculate the accelation of links
%
%		CALC_ACC returns the accelation in the Inertia frame
%		for link 1 to n.
%
%
%	'97 6/13 T.Hiraoka
%   '99 10/4 K.Yoshida

function [ vd , wd ] = calc_acc( A0, AA, w0, ww, vd0, wd0, q, qd, qdd )

global BB J_type
global c0 cc
global num_q Ez


% If Single body
if num_q == 0
   
   vd = [];
   wd = [];
   
% If Multi body system
else
   
   % Calcuration of coordinate transfromation matrices
   A_I_0 = A0;
   
   % Calculation of acceletion vectors vd,wd
   for i = 1 : num_q
      
      % Check the link connection: Is the lower one of this link, 0 ?
      if BB(i) == 0
         
         % If the i-th link connects with 0-th link
         A_I_i = AA(:,i*3-2:i*3);
         
         % Rotational joint
         if J_type(i) == 'R'
            wd(:,i) = wd0(:) ...
               + cross(ww(:,i),(A_I_i*Ez*qd(i))) ...
               + A_I_i*Ez*qdd(i);
            vd(:,i) = vd0(:) ...
               + cross( wd0(:),(A_I_0*c0(:,i)) ) ...
               + cross( w0(:),cross(w0(:),(A_I_0*c0(:,i))) )  ...
               - cross( wd(:,i),(A_I_i*cc(:,i,i)) ) ...
               - cross( ww(:,i), cross(ww(:,i),(A_I_i*cc(:,i,i))) );
            
         % Prismatic joint
         else
            wd(:,i) = wd0(:);
            vd(:,i) = vd0(:) ...
               + cross( wd0(:),(A_I_0*c0(:,i)) ) ...
               + cross( w0(:),cross(w0(:),(A_I_0*c0(:,i))) ) ...
               + cross( wd(:,i),(A_I_i*Ez*q(i)) ) ...
               + cross( ww(:,i),cross(ww(:,i),(A_I_i*Ez*q(i))) ) ...
               + 2*cross(ww(:,i),(A_I_i*Ez*qd(i))) + (A_I_i*Ez*qdd(i)) ...
               - cross( wd(:,i),(A_I_i*cc(:,i,i)) ) ...
               - cross( ww(:,i),cross(ww(:,i),(A_I_i*cc(:,i,i))) ) ;
            
         end
         
         
      % Current (i-th) link doesn't have connection with the 0-th link
      else
         A_I_BB = AA(:,BB(i)*3-2:BB(i)*3);
         A_I_i  = AA(:,i*3-2:i*3);
         
         % Rotational joint
         if J_type(i) == 'R'
            wd(:,i) = wd(:,BB(i)) ...
               + cross( ww(:,i),(A_I_i*Ez*qd(i)) ) + (A_I_i*Ez*qdd(i));
            vd(:,i) = vd(:,BB(i)) ...
               + cross( wd(:,BB(i)),(A_I_BB*cc(:,BB(i),i)) ) ...
               + cross( ww(:,BB(i)),cross(ww(:,BB(i)),(A_I_BB*cc(:,BB(i),i))) ) ...
               - cross( wd(:,i),(A_I_i*cc(:,i,i)) ) ...
               - cross( ww(:,i),cross(ww(:,i),(A_I_i*cc(:,i,i))) );
         
         % Prismatic joint
         else
            wd(:,i) = wd(:,BB(i));
            vd(:,i) = vd(:,BB(i)) ...
               + cross( wd(:,BB(i)),(A_I_BB*cc(:,BB(i),i)) ) ...
               + cross( ww(:,BB(i)),cross(ww(:,BB(i)),(A_I_BB*cc(:,BB(i),i))) ) ...
               + cross( wd(:,i),(A_I_i*Ez*q(i)) ) ...
               + cross( ww(:,i),cross(ww(:,i),(A_I_i*Ez*q(i))) ) ...
               + 2*cross( ww(:,i),(A_I_i*Ez*qd(i)) ) + (A_I_i*Ez*qdd(i)) ...
               - cross( wd(:,i),(A_I_i*cc(:,i,i)) ) ...
               - cross( ww(:,i),cross(ww(:,i),(A_I_i*cc(:,i,i))) );
            
         end
         
      end
      
   end
   
end


%%%EOF
