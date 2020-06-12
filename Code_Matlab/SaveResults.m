
% Save the results of the integration

function SaveResults(time,R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau)

global SE m0 m inertia0 inertia Gravity d_time x_plane
global fidPOS fidVEL fidACC fidFOR fidCON fidENE
global Pcon Fcon delta delta1 delta1t tau_lim
global FFs TTs Fjs Tjs FF0s TT0s

Q0 = dc2rpy(A0');          % Main body orientation
AA = calc_aa(A0,q);        % Joints orientations                       
num_q = length(q);
Qi = zeros(3,num_q);      % Joint/centroid RPY frames wrt inertia frame
for i = 1:num_q
  Qi(:,i) = dc2rpy( (AA(:,i*3-2:i*3))' ); 
end
RR = calc_pos(R0,A0,AA,q);                                 % Links centroids position
[vv,ww] = calc_vel(A0,AA,v0,w0,q,qd);                      % Links centroids velocities
[vd0,wd0,qdd] = f_dyn(R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);   % Main body and joint accelerations  
[vd,wd] = calc_acc(A0,AA,w0,ww,vd0,wd0,q,qd,qdd);          % Links centroids accelerations

% End-effector position and linear velocity
Pe = zeros(3,num_q);
vve = zeros(3,num_q);
for i = 1:length(SE)
  if ( SE(i) == 1 )
    joints = j_num( sum(SE(1:i)) );
    [Pe(:,i),~] = f_kin_e(RR,AA,joints);
    vve(:,i) = vv(:,i) - tilde( Pe(:,i)-RR(:,i) )*ww(:,i);
  end
end

% Joint positions
Pj = zeros(3,num_q);
for i = 1:length(SE)
  if ( SE(i) == 1 )
    joints = j_num( sum(SE(1:i)) );
    [Ptmp,~] = f_kin_j(RR,AA,q,joints);
    for k = 1:length(joints)
      Pj(:,joints(k)) = Ptmp(:,k);
    end
  end
end

% CoM
CoM  = calc_pos_CoM( R0, RR);
CoMv = calc_vel_CoM( v0, vv);
CoMa = calc_acc_CoM( vd0, vd);

% Inertial and joint forces (FFs, TTs, Fjs, Tjs, FF0s, TT0s)
Force = r_ne(R0,RR,A0,AA,v0,w0,vd0,wd0,q,qd,qdd,Fe,Te);

% Kinetic energy
[TK0,TKj] = calc_Kin_Ene( A0, AA, v0, w0, vv, ww );
% Potential energy (gravity)
[VG0,VGj] = calc_Pot_Ene( R0, RR );
%Work forces/moments (base, end effector, joints)
[PF0,PFe,Ptau] = calc_Work(v0,w0,vve,ww,qd,F0,T0,Fe,Te,tau);

% Linear and angular momentum (wrt CoM) and their time-derivative
[LM0,LMj] = calc_Lin_Mom(v0,vv);
[HM0,HMj] = calc_Ang_Mom(R0,A0,RR,AA,v0,w0,vv,ww,CoM);

[LM0d,LMjd] = calc_der_Lin_Mom(vd0,vd);
[HM0d,HMjd] = calc_der_Ang_Mom(R0,A0,RR,AA,v0,w0,vv,ww,vd0,wd0,vd,wd, ...
                               CoM,CoMv);

% CoP
CoP(1) = sum( Fcon(1,:).*Pcon(2,:) ) / sum( Fcon(1,:) );
CoP(2) = sum( Fcon(1,:).*Pcon(3,:) ) / sum( Fcon(1,:) );

% ZMP using the 3D equation
xp = [ x_plane  0  0 ]';
Ftmp = (m0+sum(m))*(CoMa - Gravity);
tmp = HM0d + cross( (CoM-xp) , Ftmp );
for i = 1:num_q
  Ftmp = Ftmp - Fe(:,i);
  tmp = tmp - cross( (Pe(:,i)-xp), Fe(:,i) ) + HMjd(:,i);
end
ZMP = cross( [1 0 0]' , tmp ) / Ftmp(1);

% CMP
CMP = zeros(3,1);
for i = 1:size(Fcon,2)
  CMP = CMP + cross( (Pcon(:,i)-CoM), Fcon(:,i) );
end

% Save position data
fprintf(fidPOS,' %15.8e ', time, R0,  Q0,  q,   RR, Qi, CoM,  Pe, Pj, ...
                           CoP, ZMP(2:3), CMP );
fprintf(fidPOS,' \n ');

% Save velocity data
fprintf(fidVEL,' %15.8e ', time, v0,  w0,  qd,  vv, ww, CoMv, vve); 
fprintf(fidVEL,' \n ');

% Save acceleration data
fprintf(fidACC,' %15.8e ', time, vd0, wd0, qdd, vd, wd, CoMa);
fprintf(fidACC,' \n ');

% Save force data
fprintf(fidFOR,' %15.8e ',time, F0, T0, tau, Fe, Te, ...
                                FF0s, TT0s, FFs, TTs, Fjs, Tjs); 
fprintf(fidFOR,' \n ');

% Save contact data
fprintf(fidCON,' %15.8e ',time, Pcon, Fcon, delta, delta1, delta1t, tau_lim);            
fprintf(fidCON,' \n ');

% Save energy data
fprintf(fidENE,' %15.8e ',time, TK0, TKj, VG0, VGj, PF0, PFe, Ptau, ...
                                LM0, LMj, HM0, HMj, LM0d, LMjd, HM0d, HMjd);   
fprintf(fidENE,' \n ');

return
