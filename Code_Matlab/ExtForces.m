
% Determine all the external forces/torques 

function [F0, T0, Fe, Te, tau] = ExtForces(time,R0,A0,v0,w0,q,qd)

global m m0 inertia inertia0 q_ini qd_ini q_lim d_time R0ini
global Pc x_plane Pcon Fcon delta delta1 delta1t
global tau_save CoP_old CoP_old_old p0_des_old p1_des_old 
global Kcon Dcon Ncon mu_c mu_d min_delta1t Kp tau_lim param

% Data needed
num_q = length(q);

% Forces/torques on the end points
Fe = zeros(3,num_q);                 
Te = zeros(3,num_q); 

% Forces/torques on the main body
F0 = [ 0.0	 0.0   0.0 ]';                                
T0 = [ 0.0	 0.0   0.0 ]';
Kcon = 5e6;
Dcon = 5e6;
Ncon = 1.5;
mu_c = 0.0;
mu_d = 500.0;
min_delta1t = 0;
npc = size(Pc,2);
Pcon = zeros(3,npc);
delta = zeros(1,npc);
delta1 = zeros(1,npc);
delta1t = zeros(1,npc);
Fcon = zeros(3,npc);
for i = 1:npc
  dist = A0*Pc(:,i);
  Pcon(:,i) = R0 + dist;  % Current position
  delta(i) = x_plane - Pcon(1,i);    % Penetration
  Vcon = v0 + cross(w0,dist);  % Speed
  delta1(i) = -Vcon(1);
  delta1t(i) = -Vcon(2);
  % If positive there is contact
  if ( delta(i) > 0 )
    Fn = Kcon*delta(i)^Ncon + Dcon*delta(i)^Ncon*delta1(i);
    if ( abs(delta1t(i)) > min_delta1t ) 
      Ft = mu_c*Fn*sign(delta1t(i)) + mu_d*delta1t(i);
    else
      Ft = 0;  % Use penalty function up to static friction
    end
    Fcon(:,i) = [ Fn  Ft  0 ]';
    F0 = F0 + Fcon(:,i);
    T0 = T0 + cross(dist,Fcon(:,i));
  else
    delta(i) = 0;
    delta1(i) = 0;
    delta1t(i) = 0;
  end
end

% Joint limits
% Add tau_lim to tau to impose limits
tau_lim = zeros(num_q,1);
for i = 1:size(q_lim,1)
  idj = q_lim(i,1);
  if ( idj == 0 )
    continue
  end
  q_lim_do = q_lim(i,2);
  q_lim_up = q_lim(i,3);
  K_lim = q_lim(i,4);
  D_lim = q_lim(i,5);
  if ( q(idj) < q_lim_do )
    delta_q = q_lim_do - q(idj);
    tau_lim(idj) = + K_lim*delta_q^Ncon - D_lim*delta_q^Ncon*qd(idj);
  end
  if ( q(idj) > q_lim_up )
    delta_q = q(idj) - q_lim_up;
    tau_lim(idj) = - K_lim*delta_q^Ncon - D_lim*delta_q^Ncon*qd(idj);
  end
end

global paramPSO 
persistent R0des Q0des qdes tdes

if (time == 0)
  tdes = 0;
  qdes  = q;
  R0des = R0;
  Q0des = dc2rpy(A0');
end

% if (time <= 9)
%   qdes(4) = q_ini(4) - q(1) - pi*time/9;
% else
%   qdes(4) = q_ini(4) - q(1) - pi;
% end
qdes(4) = ( q_ini(4) - q(1) - pi*sin_ramp(time,9) );

dtPSO = 0.1;
if ( rem(time,dtPSO) == 0 )

  rng(10);
  paramPSO.type = 0;
  if (time <= 3)
    paramPSO.POS = [ 0.4666-time*0.05  0.0261-time*0.05  0 ]';
  elseif (time > 3 && time <= 6)
    paramPSO.POS = [ 0.3166+(time-3)*0.05  -0.1239  0 ]';
  elseif (time > 6 && time <= 9)
    paramPSO.POS = [ 0.4666  -0.1239-(time-6)*0.05  0 ]';
  else
%     paramPSO.POS = [ 0.4666   0.0261  0 ]';
    paramPSO.POS = [ 0.4666  -0.2739  0 ]';
  end
  bound = 0.5;
%   % All joints limited
%   LB = [ R0ini 0 0    0 0 0    -pi/6 -pi/2 -pi -pi ];
%   UB = [ R0ini 0 0    0 0 0    +pi/6 +pi/2 +pi +pi ];
%   % All joints near-limited
%   LB = [ R0ini R0des(2)-0.01 0    0 0 0    (qdes-bound)' ];
%   UB = [ R0ini R0des(2)+0.01 0    0 0 0    (qdes+bound)' ];
%   % Joint 4 horizontal
%   LB = [ R0ini R0des(2)-bound/100 0    0 0 0    (qdes(1:3)-bound)'  q_ini(4)-q(1) ];
%   UB = [ R0ini R0des(2)+bound/100 0    0 0 0    (qdes(1:3)+bound)'  q_ini(4)-q(1) ];
%   % Joint 4 rotating 180 deg
  LB = [ R0ini R0des(2)-bound/100 0    0 0 0    (qdes(1:3)-bound)'  qdes(4) ];
  UB = [ R0ini R0des(2)+bound/100 0    0 0 0    (qdes(1:3)+bound)'  qdes(4) ];
  [X,Fpso] = PSO(@inv_Kin_Func,LB,UB,1,10,1e-8);
  if ( length(Fpso) == 100 )
    time
    Fpso(end)
  end

  tdes = time;
  R0des = X(1:3)';
  Q0des = X(4:6)';
  qdes  = X(7:end)';

end

% Control
Kp = 1000*[ 1  1  0  0 ]';
Kd = 2*0.5*sqrt(m.*Kp);
%tau = - Kp.*( q - qdes ) - Kd.*qd;
tau = - Kp.*( q - q_ini ) - Kd.*qd;

% Kp0 = 1000;
% Kd0 = 2*0.5*sqrt(m0*Kp0);
% tau0_2 = - Kp0*( R0(2) - R0des(2) ) - Kd0*v0(2);
% F0 = F0 + [ 0  tau0_2*sin_ramp(time-tdes,0.1)  0 ]';

tau = tau + tau_lim;

return

% KL = 0.1;
% DL = 0.1;
% c0_des = [ 0.51  0.045  0 ]';
% c1_des = [ 0  0  0 ]';
% bound = 10;
% dt_ctrl = 0.02;
% 
% KH = 0.1;
% DH = 0.1;
% p0_des = 0;
% p1_des = 0;
% W = [ 1  1 ];
% 
% tau = tau_save;
% 
% AA = calc_aa(A0,q);
% RR = calc_pos(R0,A0,AA,q);
% [vv,ww] = calc_vel(A0,AA,v0,w0,q,qd);
% [vd0,wd0,qdd] = f_dyn(R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);
% mass = m0 + sum(m);
% CoM  = calc_pos_CoM( R0, RR);
% CoMv = calc_vel_CoM( v0, vv);
% c2_des = KL*( c0_des - CoM ) + DL*( c1_des - CoMv );
% L1_des = mass*c2_des;
% 
% CoP = Pcon(2,2) + (Pcon(2,1)-Pcon(2,2))*Fcon(1,1)/(Fcon(1,1)+Fcon(1,2));
% if ( time == 0 )
%   CoPv = 0;
% elseif ( time == d_time )
%   CoPv = ( CoP - CoP_old )/d_time;
% else
%   CoPv = (3*CoP - 4*CoP_old + CoP_old_old)/(2*d_time);
% end
% CoP_old_old = CoP_old;
% CoP_old = CoP;
% 
% p2_des = KH*( p0_des - CoP ) + DH*( p1_des - CoPv );
% p1_des = p1_des_old + p2_des*d_time;
% p0_des = p0_des_old + p1_des_old*d_time;
% H1_des = ( p0_des - CoM(2) )*F0(1);
% p0_des_old = p0_des;
% p1_des_old = p1_des;
% 
% if (rem(time,dt_ctrl) ~= 0)
%   return
% end
% 
% param.A0 = A0;
% param.AA = AA;
% param.w0 = w0;
% param.ww = ww;
% param.q = q;
% param.qd = qd;
% param.vd0 = vd0;
% param.wd0 = wd0;
% param.L1_des = L1_des;
% param.qdd = qdd;
% param.H1_des = H1_des;
% param.W = W;
% param.R0 = R0;
% param.RR = RR;
% param.v0 = v0;
% param.vv = vv;
% param.CoM = CoM;
% param.CoMv = CoMv;
% 
% rng(10);
% LB = qdd' - bound;
% UB = qdd' + bound;
% [X,Fpso] = PSO(@func,LB,UB,100,50,1e-4);
% %[X,Fpso,exitflag,output] = particleswarm(@func,num_q,LB,UB);
% qdd = X;
% Force = r_ne(R0,RR,A0,AA,v0,w0,vd0,wd0,q,qd,qdd,Fe,Te);
% tau = Force(7:end);
% tau_save = tau;
% 
% 
% fprintf('%6.3f   %4d    %10.4e',time,length(Fpso),Fpso(end));
% %fprintf('%6.3f   %4d    %10.4e',time,output.iterations,Fpso(end));
% fprintf(' \n');



% Control based on angular momentum
% AA = calc_aa(A0,q);
% RR = calc_pos(R0,A0,AA,q);
% [vv,ww] = calc_vel(A0,AA,v0,w0,q,qd);
% CoM  = calc_pos_CoM(R0,RR);
% [HM0,HMj] = calc_Ang_Mom(R0,A0,RR,AA,v0,w0,vv,ww,CoM);
% Htot = HM0 + sum(HMj,2);

%Kp4 = 1.1;
%Kd4 = 2*0.5*sqrt(m(4)*Kp4);
%tau(4) = - Kp4*( HMj(3) + HMj(3) );
% tau(4) = - Kp4*( H1tot(3) ) ;
 
% Kp3 = 4;
% Kd3 = 2*0.5*sqrt(m(3)*Kp3);
% if (time < 2)
%   q3f = +pi/2;
% else
%   q3f = -pi/2;
% end
% qd3f = 0;
%tau(3) = - Kp3*( q(3) - q3f ) - Kd3*( qd(3) - qd3f );
%tau(3) = tau(3)*sin_ramp(min(time,1),0.1);


end



% EOF