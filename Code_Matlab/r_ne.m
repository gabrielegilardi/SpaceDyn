%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The SpaceDyn, a MATLAB toolbox for Space and Mobile Robots.
% (C)1998 The Space Robotics Lab. directed by Kazuya Yoshida,
% Tohoku University, Japan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   version 1.2 // Jul. 20, 2006, Last modification by H.Nakanishi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% R_NE		Inverse Dynamic computation by the Recursive Newton-Euler method
%
%		R_NE returns a generalized force, which consists of
%		the reaction forces FF0, TT0 on the link 0, and torque tau
%		of each joint.
%
%
%		'97 11/25 T.Hiraoka
%		'98 2/13 K.Fujishima
%       '00 9/6 H.Hamano
%       '06 7/20 H.Nakanishi

function Force = r_ne( R0, RR, A0, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te )

global S0 SS SE J_type
global c0 cc ce
global m0 m inertia0 inertia
global Ez Gravity
global FFs TTs Fjs Tjs FF0s TT0s

% Number of links
num_q = length( q );

% Calculation of coordinate transfromation matrices
A_I_0 = A0;

% Calculation of velocity vectors vv,ww
% NOTE:	vv,ww are in the Inertial frame

[ vv,ww ] = calc_vel( A0, AA, v0, w0, q, qd );

% Calculation of acceletion vectors vd,wd
% NOTE:	vd,wd are in the Inertial frame

[ vd , wd ] = calc_acc( A0, AA, w0, ww, vd0, wd0, q, qd, qdd );


% Calculation of inertia force & torque of link 0
% NOTE:	FF,TT(FF0,TT0) are in the Inertial frame.

FF0 = m0 * (vd0-Gravity);
TT0 = (A_I_0*inertia0*A_I_0')*wd0 ...
   + cross( w0 , ((A_I_0*inertia0*A_I_0')*w0) );

% Calculation of inertia force & torque of link i
% from link 1 to n
% Single or multi body ?

if num_q == 0
   % If a Single body,
   FF = [];
   TT = [];
   
   % If a Multi body system,
else
   
   for i = 1 : num_q
      A_I_i  = AA(:,i*3-2:i*3);
      in_i = inertia(:,i*3-2:i*3);
      FF(:,i) = m(i) * (vd(:,i)-Gravity) ;
      TT(:,i) = (A_I_i*in_i*A_I_i')*wd(:,i) ...
         + cross( ww(:,i) , ((A_I_i*in_i*A_I_i')*ww(:,i)) );
   end
   
end

TT


% Equilibrium of forces & torques on each link
% On the i-th link

Fj = zeros(3,num_q);
Tj = zeros(3,num_q);

if (num_q ~= 0)
   
   % Multi body system
   % from link n to 1
   
   for i = num_q : -1 : 1
      
      F = zeros(3,1);
      T = zeros(3,1);
      
      for j=i+1:num_q
         
         F =  F + SS(i,j)*Fj(:,j);
         
      end
      
      Fj(:,i) = FF(:,i) + F - SE(i)*Fe(:,i);
      
      for j=i+1:num_q
         
         A_I_i = AA(:,i*3-2:i*3);
         T =  T ...
            + SS(i,j)*( cross(A_I_i*(cc(:,i,j)-cc(:,i,i)+(J_type(i)=='P')*Ez*q(i)),Fj(:,j) ) + Tj(:,j) );
         
      end
      
      if J_type(i) == 'R'
         % Rotational joint
         Tj(:,i) = TT(:,i) + T ...
            - cross( A_I_i*cc(:,i,i),FF(:,i) ) ;
         
      else
         % Prismatic joint
         Tj(:,i) = TT(:,i) + T ...
            + cross( A_I_i*(Ez*q(i))-A_I_i*cc(:,i,i) , FF(:,i) );
         
      end
      
      Tj(:,i) = Tj(:,i) ...
         - SE(i)*( cross(A_I_i*(ce(:,i)-cc(:,i,i)+(J_type(i)=='P')*Ez*q(i)),Fe(:,i) ) + Te(:,i) );
      
   end
   
   Fj
   
   % Equilibrium on the link 0
   
   F = zeros(3,1);
   T = zeros(3,1);
   
   for i=1:num_q
      
      if (S0(i) ~= 0)
         F =  F + S0(i)*Fj(:,i);
         
      end
      
   end
   
   FF0 = FF0 + F;
   
   for i=1:num_q
      
      if (S0(i) ~= 0)
         T = T + S0(i)*( cross( (A_I_0*c0(:,i)) , Fj(:,i) ) + Tj(:,i) );
         
      end
      
   end
   
   TT0 = TT0 + T;
   
end


% Calculation of torques of each joint

% Single body

if num_q == 0
   tau = zeros(0);
   
else
   % Multi body system
   
   for i = 1 : num_q
      
      A_I_i = AA(:,i*3-2:i*3);
      
      if J_type(i) == 'R'
         % Rotational joint
         tau(i,1) = Tj(:,i)'*(A_I_i*Ez);
         
      else
         % Prismatic joint
         tau(i,1) = Fj(:,i)'*(A_I_i*Ez);
         
      end
      
   end
   
end


% Compose a generalized force

% Single or multi body ?

if num_q == 0
   % Single body,
   Force = [ FF0' TT0' ]';
   
else
   % Multi body system,
   Force = [ FF0' TT0' tau' ]';
   
end

FFs = FF;
TTs = TT;
Fjs = Fj;
Tjs = Tj;
FF0s = FF0;
TT0s = TT0;


%%%EOF
