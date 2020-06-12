
% Links

%   0     trunk (TR)
%   1     left thigh (LT)
%   2     left shank (LS)
%   3     left foot (LF)
%   4     right thigh (RT)
%   5     right shank (RS)
%   6     right foot (RF)

% Lengths

%   2d      thigh
%   2d      shank
%   7d/3    trunk (to shoulders)
%    d/3    foot thikness
% d/3,2d/3  foot length

function robotWalk(prefix)

global BB S0 SS SE J_type
global c0 cc ce Qi Qe
global m0 m inertia0 inertia
global num_q Ez Gravity d_time

global fidPOS fidVEL fidACC fidFOR fidCON fidENE
global Pc x_plane Kcon Dcon Ncon mu_c mu_d q_ini qd_ini Kp min_delta1t

tic

times = [ 1  2e-4  10 ];
% times = [ 1.0  0.0002  10  ;
%           4.0  0.0002  10 ];
q_ini  = [ +pi/6  -pi/6  0  -pi/6  -pi/6  0 ]';
qd_ini = [ 0  0  0  0  0  0 ]';

Kp = 1000*[ 1  1  1  1  1  1 ]';
Kcon = 5e6;
Dcon = 5e6;
Ncon = 1.5;
mu_c = 0.0;
mu_d = 0.0;
min_delta1t = 0;

% 1.8m = heigh
height = 1.80;
BMI = 25;
d = height/8;
mass = BMI*height^2;

Gravity = [ -9.81  0  0 ]';
Ez = [0 0 1]';

R0 = [ 167*d/30  0  0 ]'; 
Q0 = [ 0  0  0 ]';
v0 = [ 0  0  0 ]';      
w0 = [ 0  0  0 ]';                                

% ===============================================

% Layout    1   2   3   4   5   6 
J_type = [ 'R' 'R' 'R' 'R' 'R' 'R' ];
BB     = [  0   1   2   0   4   5  ];
num_q = length(BB);

[S0_temp,SS_temp,SE_temp] = BuildIncidence(BB);
S0 = S0_temp;
SS = SS_temp;
SE = SE_temp;

% Mass and inertia

%   50%   Trunk
%   7%    Neck (20%) + Head (80%)
%   10%   Thigh
%   5%    Shank
%   2.7%  Upper arm
%   2.3%  Forearm (80%) + hand (20%)
%   1.5%  Foot
m0 = mass*0.5;    % TR
m = mass*[ 0.1 0.05 0.015 0.1 0.05 0.015 ]';
%          LT  LS   LF    RT  RS   RF    LUA   LFA   RUA   RFA   NH

% Main body inertia
inertia0 = m0*[ 0.01  0     0     ;
                0     0.01  0     ;
                0     0     0.01 ];
              
% Inertia (wrt link/joint frame)
inertia = zeros(3,3*num_q);   

% Left leg (thigh, shank, foot)
In = mom_inertia(1,[ 0.05  0.9*0.05  2*d ]);
inertia(:,1:3) = m(1)*In;

In = mom_inertia(1,[ 0.05  0.9*0.05  2*d ]);
inertia(:,4:6) = m(2)*In;

In = mom_inertia(1,[ d/6  0.9*d/6  d ]);
inertia(:,7:9) = m(3)*In;

                 
% Right leg (thigh, shank, foot)

In = mom_inertia(1,[ 0.05  0.9*0.05  2*d ]);
inertia(:,10:12) = m(4)*In;

In = mom_inertia(1,[ 0.05  0.9*0.05  2*d ]);
inertia(:,13:15) = m(5)*In;

In = mom_inertia(1,[ d/6  0.9*d/6  d ]);
inertia(:,16:18) = m(6)*In;

% System geometry 

% Main body centroid to connected joints (S0 = +1)
c0 = zeros(3,num_q);        
c0(:, 1)  = [ -7*d/6  0  0 ]';       % Left thigh
c0(:, 4)  = [ -7*d/6  0  0 ]';       % Right thigh

% Links centroid to end points (SE = +1)
ce = zeros(3,num_q);
ce(:, 3) = [ 0  0  0 ]';       
ce(:, 6) = [ 0  0  0 ]';       

% Links centroid to own joints (SS = -1)
cc = zeros(3,num_q,num_q ); 
cc(:, 1, 1) = [ +d  0  0 ]';  
cc(:, 2, 2) = [ +d  0  0 ]';       
cc(:, 3, 3) = [ +d/5  0  0 ]';  
cc(:, 4, 4) = [ +d  0  0 ]';  
cc(:, 5, 5) = [ +d  0  0 ]';       
cc(:, 6, 6) = [ +d/5  0  0 ]';  
% Links centroid to connected joints (SS = +1)
cc(:, 1, 2) = [ -d  0  0 ]';  
cc(:, 2, 3) = [ -d  0  0 ]';  
cc(:, 4, 5) = [ -d  0  0 ]';  
cc(:, 5, 6) = [ -d  0  0 ]';  

% RPY orientations
Qi = zeros(3,num_q);        % Links/joints (wrt previous one)
Qe = zeros(3,num_q);        % End points (wrt link/joint frame)

% Contact points position (wrt centroid link 3 and 4)
Pc = [ -d/5    -d/5    -d/5    -d/5    ;
       +d/3    -2*d/3  +d/3    -2*d/3  ;
        0       0       0       0     ];
x_plane = 0;
      
% =============================================

% Output files
fidPOS = fopen(strcat(prefix,'.pos'),'wt');      
fidVEL = fopen(strcat(prefix,'.vel'),'wt');      
fidACC = fopen(strcat(prefix,'.acc'),'wt');      
fidFOR = fopen(strcat(prefix,'.for'),'wt');
fidCON = fopen(strcat(prefix,'.con'),'wt');
fidENE = fopen(strcat(prefix,'.ene'),'wt');

save( strcat(prefix,'.mat'),'BB','S0','SE','SS','J_type','c0','cc','ce', ...
                            'm0','m','inertia0','inertia','Ez','Gravity', ...
                            'Pc','x_plane','prefix','Kcon','Dcon','Ncon', ...
                            'mu_c','mu_d','q_ini','qd_ini','Kp','d', ...
                            'height','min_delta1t','times','XcmLFA', ...
                            'XcmRFA','XcmNH','mass','BMI');

% Save solution at t_start
A0 = rpy2dc(Q0)';                                   
q  = q_ini;
qd = qd_ini;
t_start = 0;
[F0,T0,Fe,Te,tau] = ExtForcesFull(t_start,R0,A0,v0,w0,q,qd);  % Ext. forces 
SaveResults(t_start,R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);   % Save

% Loop over time
n_times = size(times,1);
for k = 1:n_times
  t_end = times(k,1);
  d_time = times(k,2);
  saveStep = times(k,3);  
  countSave = 0;
  for time = t_start+d_time:d_time:t_end
    
    [R0,A0,v0,w0,q,qd] = f_dyn_nb2_fixed(R0,A0,v0,w0,q,qd,F0,T0,Fe,Te,tau);  
    [F0,T0,Fe,Te,tau] = ExtForcesFull(time,R0,A0,v0,w0,q,qd);   % Ext. forces
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

