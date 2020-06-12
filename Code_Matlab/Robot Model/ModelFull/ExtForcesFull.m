
% Determine all the external forces/torques 

function [F0, T0, Fe, Te, tau] = ExtForcesFull(time,R0,A0,v0,w0,q,qd)

global m m0
global Pc x_plane Pcon Fcon delta delta1 delta1t
global Kcon Dcon Ncon mu_c mu_d Kp min_delta1t q_ini qd_ini

% Initialize
num_q = length(q);
F0 = [ 0.0	 0.0   0.0 ]';    % Forces/torques on the main body
T0 = [ 0.0	 0.0   0.0 ]';
Fe = zeros(3,num_q);          % Forces/torques on the end points
Te = zeros(3,num_q);  
tau = zeros(num_q,1);         % Torques on the joints

% Forces on the feet (links 3 and 4) due to contact
npc = size(Pc,2);
Pcon = zeros(3,npc);
delta = zeros(1,npc);
delta1 = zeros(1,npc);
delta1t = zeros(1,npc);
Fcon = zeros(3,npc);
AA = calc_aa(A0,q);
RR = calc_pos(R0,A0,AA,q);
[vv,ww] = calc_vel(A0,AA,v0,w0,q,qd);
for i = 1:4
  if ( i <= 2 )
    idx = 3;
  else
    idx = 6;
  end
  k = 3*(idx-1)+1;
  dist = AA(:,k:k+2)*Pc(:,i);
  Pcon(:,i) = RR(:,idx) + dist;  % Current position
  delta(i) = x_plane - Pcon(1,i);    % Penetration
  Vcon = vv(:,idx) + cross(ww(:,idx),dist);  % Speed
  delta1(i) = -Vcon(1);
  delta1t(i) = -Vcon(2);
  % If positive there is contact
  if ( delta(i) > 0 )
    Fn = Kcon*delta(i)^Ncon + Dcon*delta(i)^Ncon*delta1(i);
    if ( abs(delta1t(i)) > min_delta1t ) 
      Ft = mu_c*Fn*sign(delta1t(i)) + mu_d*delta1t(i);
    else
      Ft = 0;
    end
    Fcon(:,i) = [ Fn  Ft  0 ]';
    Fe(:,idx) = Fe(:,idx) + Fcon(:,i);
    Te(:,idx) = Te(:,idx) + cross(dist,Fcon(:,i));
  else
    delta(i) = 0;
    delta1(i) = 0;
    delta1t(i) = 0;
  end
end

% Torques on the joints
Kd = 2*0.5*sqrt(m.*Kp);
tau = - Kp.*( q - q_ini ) - Kd.*( qd - qd_ini );

% EOF