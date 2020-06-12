

function robot(prefix)

global BB S0 SS SE J_type Ez 
global c0 cc ce Qi Qe
global m0 m inertia0 inertia Gravity
global num_q d_time

global fidPOS fidVEL fidACC fidFOR fidCON fidENE
global Pc x_plane q_ini qd_ini q_lim R0ini
global tau_save CoP_old CoP_old_old p0_des_old p1_des_old 
global Kcon Dcon Ncon mu_c mu_d min_delta1t Kp

tic

% times = [ 1.0  2e-4  10  ;
%           3.0  2e-4  10 ];

times = [ 1  0.0002  4 ];
La = 0.4;
Lu = 0.4;
Ld = 0.4;
d = 0.05;
delta_ini = -8e-4;

% q_ini  = [ 0  0  pi  pi/2 ]';
q_ini = [0.1 0.5 1.0 1.5]';
% qd_ini = [ 0  0  0  0 ]';
qd_ini = [ 0.01  0.02  0.03  0.04 ]';

% Joint limits
% Data: id, lower limit, upper limit, stiffness, damping
% If id=0 the limits are ignored
% K = 5e6, D = 5e6 (0.1 deg. overshooting, no residual)
% K = 1e4, D = 1e4 (1.4 deg. overshooting, 0.3 residual)
q_lim = [ 0  -pi/3  +pi/3  5e6  5e6  ;
          0  -pi/2  +pi/2  5e6  5e6  ;
          0  -pi    +pi    5e6  5e6  ;
          0  -pi    +pi    5e6  5e6 ];

mt3 = 0;
mt4 = 0;

Gravity = [ -9.81  0  0 ]';
Ez = [0 0 1]';

% ===============================================

% Layout    1   2   3   4
J_type = [ 'R' 'R' 'R' 'R' ];
BB     = [  0   1   2   1  ];
num_q = length(BB);

[S0_temp,SS_temp,SE_temp] = build_Inc_Mat(BB);
S0 = S0_temp;
SS = SS_temp;
SE = SE_temp;

% Main body data
R0ini =   d + delta_ini;
Q0 = [ 0.1  0.2  0.3 ]';
R0 = [ 1.0  2.0  3.0 ]'; 
% R0 = [ R0ini  0  0 ]'; 
A0 = rpy2dc(Q0)';                                   
v0 = [ 0  0  0 ]';       % Initial linear and angular velocity
w0 = [ 0  0  0 ]';                                
m0 = 3;
In = mom_inertia(1,[ d  d*0.9  6*d ]);
inertia0 = m0*In;

% Joint initial position and velocity
q  = q_ini;
qd = qd_ini;

% Link mass/inertia
m = [ 7 7 3 3 ]';
inertia = zeros(3,3*num_q);   % Inertia (wrt link/joint frame)
% Link 1
In = mom_inertia(1,[ 0.05  0.9*0.05  Ld ]);
inertia(:,1:3) = m(1)*In;

% Link 2                        
In = mom_inertia(1,[ 0.05  0.9*0.05  Lu ]);
inertia(:,4:6) = m(2)*In;

% Link 3                        
Xcm3 = (m(3)*La/2 + mt3*La)/(m(3) + mt3);
In = mom_inertia(1,[ 0.05  0.9*0.05  La ]);
inertia(1,7) = m(3)*In(1,1);
inertia(2,8) = m(3)*In(2,2) + m(3)*(Xcm3-La/2)^2 + mt3*(La-Xcm3)^2;
inertia(3,9) = m(3)*In(3,3) + m(3)*(Xcm3-La/2)^2 + mt3*(La-Xcm3)^2;
m(3) = m(3) + mt3; 

% Link 4
Xcm4 = (m(4)*La/2 + mt4*La)/(m(4) + mt4);

In = mom_inertia(1,[ 0.05  0.9*0.05  La ]);
inertia(1,10) = m(4)*In(1,1);
inertia(2,11) = m(4)*In(2,2) + m(4)*(Xcm4-La/2)^2 + mt4*(La-Xcm4)^2;
inertia(3,12) = m(4)*In(3,3) + m(4)*(Xcm4-La/2)^2 + mt4*(La-Xcm4)^2;
m(4) = m(4) + mt4; 

% Geometric vectors
c0 = zeros(3,num_q);        % Main body centroid to connected joints
c0(:,1)  = [ +d  0  0 ]';

ce = zeros(3,num_q);        % Link centroid to end point
ce(:,3) = [ (La-Xcm3)  0  0 ]';       
ce(:,4) = [ (La-Xcm4)  0  0 ]';       

cc = zeros(3,num_q,num_q ); % Link centroid to connected joints
cc(:,1,1) = [ -Ld/2  0  0 ]';  
cc(:,1,2) = [ +Ld/2  0  0 ]';  
cc(:,1,4) = [ +Ld/2  0  0 ]';  
cc(:,2,3) = [ +Lu/2  0  0 ]';  
cc(:,2,2) = [ -Lu/2  0  0 ]';       
cc(:,3,3) = [ -Xcm3  0  0 ]';  
cc(:,4,4) = [ -Xcm4  0  0 ]';  

% RPY orientations
Qi = zeros(3,num_q);        % Links/joints (wrt previous one)
Qe = zeros(3,num_q);        % End points (wrt link/joint frame)

% Contact parameters (Pc is given wrt the centroid of idCon)
Pc = [ -d    -d     ;
       +3*d  -3*d   ;
        0     0    ];

x_plane = 0;

Qe(:,1) = [0.1, 0.2, 0.3]';
Qi(:,1) = [0.1, 0.2, 0.3]';
Q0 = [ 0.1  0.2  0.3 ]';
R0 = [ 1.0  2.0  3.0 ]'; 
v0 = R0;
w0 = Q0;
m0 = 2.3;
m = [0.4, 1.4, 0.7, 0.8]';
inertia0 = 1.5*eye(3);
inertia = [0.5*inertia0, 0.3*inertia0, 2.3*inertia0, 1.5*inertia0];
A0 = rpy2dc(Q0)';                                   

AA = calc_aa( A0 , q );
RR = calc_pos( R0, A0, AA, q );
[vv,ww] = calc_vel( A0, AA, v0, w0, q, qd );
vd0 = [0.2, 0.2, 0.2];
wd0 = [0.3, 0.3, 0.3];
qdd = [-0.1, -0.2, -0.3, -0.4];
[ vd , wd ] = calc_acc( A0, AA, w0, ww, vd0, wd0, q, qd, qdd );
joints = j_num( 1 );
Jacobian = calc_je( RR, AA, q, joints );
JJ_t = calc_jt( RR, AA );
JJ_r = calc_jr( AA );
HH = calc_hh( R0, RR, A0, AA );
vd0 = zeros(3,1);
wd0 = zeros(3,1);
qdd = zeros(num_q,1);
Fe = zeros(3,num_q);
Te = zeros(3,num_q);
Force = r_ne( R0, RR, A0, AA, v0, w0, vd0, wd0, q, qd, qdd, Fe, Te );
save('temp')

return

% =============================================

% Output files
fidPOS = fopen(strcat(prefix,'.pos'),'wt');      
fidVEL = fopen(strcat(prefix,'.vel'),'wt');      
fidACC = fopen(strcat(prefix,'.acc'),'wt');      
fidFOR = fopen(strcat(prefix,'.for'),'wt');
fidCON = fopen(strcat(prefix,'.con'),'wt');
fidENE = fopen(strcat(prefix,'.ene'),'wt');

% Basic data
save( strcat(prefix,'.mat'),'BB','S0','SE','SS','J_type','c0','cc','ce', ...
                            'm0','m','inertia0','inertia','Ez','Gravity', ...
                            'Pc','x_plane','prefix','Kcon','Dcon','Ncon', ...
                            'mu_c','mu_d','q_ini','mt3','qd_ini','Kp', ...
                            'min_delta1t','times','Xcm3','Xcm4','mt4', ...
                            'La','Lu','Ld','d','R0ini');


% Save solution at t_start
t_start = 0;
tau_save = zeros(num_q,1);
CoP_old = 0;
CoP_old_old = 0;
p0_des_old = 0;
p1_des_old = 0;
[F0,T0,Fe,Te,tau] = ExtForces(t_start,R0,A0,v0,w0,q,qd);  % Ext. forces
SaveResults(t_start,R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);   % Save


% Loop over time
dt_print = 0.1;
n_times = size(times,1);
for k = 1:n_times
  t_end = times(k,1);
  d_time = times(k,2);
  saveStep = times(k,3);  
  countSave = 0;
  for time = t_start+d_time:d_time:t_end
    
    if (rem(time,dt_print) == 0)
      fprintf(' %6.3f \n',time);
    end

    % Dynamics
    [R0,A0,v0,w0,q,qd] = f_dyn_nb2(R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);
%     Q0 = [ 0  0  0 ]';
%     R0 = [ R0ini  0  0 ]'; 
%     A0 = rpy2dc(Q0)';                                   
%     v0 = [ 0  0  0 ]';       
%     w0 = [ 0  0  0 ]';                                

    [F0,T0,Fe,Te,tau] = ExtForces(time,R0,A0,v0,w0,q,qd);   % Ext. forces
    countSave = countSave + 1;
    if ( countSave == saveStep )
      SaveResults(time,R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);    % Save
      countSave = 0;
    end
    
  end
  t_start = t_end;
  
end

% Close output files
fclose(fidPOS);
fclose(fidVEL);
fclose(fidACC);
fclose(fidFOR);
fclose(fidCON);
fclose(fidENE);

toc

end
