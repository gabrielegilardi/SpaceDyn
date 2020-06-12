%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Sep. 17, 2001, Last modification by H.Hamano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% I_KINE	Inverse Kinematics
%		Calculate joint angle solution (q_sol) for
%		given end point position/orientation.
%		It takes an iterative approach from a specified
%		initial posture (q_init). 
%               
%
% 1998.1.13 (C)Space Robotics Lab, by Akihide Kurosu
% 1998.9.16 K.Yoshida
% 2001.9.17 H.Hamano

function q_sol = i_kine( R0, A0, POS_e, ORI_e, q_init, num_e )


% Check matrix size of POS_e and ORI_e
% if (size(POS_e)~=[3,1]) || (size(ORI_e)~=[3,3])
%    error('Dimensiones of Position and Orientation does not meet the requirement');
% end

% Set some values
loop_limit = 1000;
norm_limit = 1e-8;
nm    = 1;
count = 0;
gain  = 0.1*ones(6,1);
q = q_init;

% Start convergent calculation
while ( nm > norm_limit )
   
   % Calculate the orientation/position of all link
   AA = calc_aa( A0, q );
   RR = calc_pos( R0, A0, AA, q );
   
   % Calculate the joint connection from base to endpoint
   joints = j_num( num_e );
   
   % Present endpoint position/orientation (forward kinematics)
   [now_p, now_o] = f_kin_e( RR, AA, joints );
   
   % Calculate the error between present and goal position/orientation
   err_p = POS_e - now_p;
   err_o = tr2diff(ORI_e, now_o);
   err = [ err_p; err_o ];
   err = err.*gain;
   
   % Calculate the Jacobian matrix
   Jacob = calc_je( RR, AA, q, joints );
   
   % Calculate joint angular velocity using Jacobian
   qd = pinv( Jacob ) * err;
   
   % Next Joint angle
   q = q + qd;
   
   % Calculate the norm of joint velocity
   nm = norm(qd);
   
   % Check the loop condition
   count = count + 1;
   if ( count > loop_limit )
      error('Solution does not converge');
   end
   
end

q_sol = q;

%%%EOF
