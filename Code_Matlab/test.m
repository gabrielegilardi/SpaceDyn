
% Testing script

fprintf('\n==========  Test cx, cy, cz  ==========');
theta = 30*pi/180;
Cx = cx(theta)
Cy = cy(theta)
Cz = cz(theta)

fprintf('==========  Test tilde  ==========');
a = [1.3, -2.3, 0.7];
b = tilde(a)

fprintf('==========  Test rpy2dc  ==========');
rpy = [1.3, -2.3, 0.7];
roll = rpy(1);
pitch = rpy(2);
yaw = rpy(3);
C = rpy2dc(rpy)
C = rpy2dc(roll,pitch,yaw)

fprintf('==========  Test eul2dc  ==========');
eul = [1.3, -2.3, 0.7];
phi = eul(1);
theta = eul(2);
psi = eul(3);
C = eul2dc(rpy)
C = eul2dc(phi,theta,psi)

fprintf('==========  Test cross  ==========');
u = [1.3, -2.3,  0.7];
v = [0.5,  1.2, -3.2];
n=cross(u,v)

fprintf('==========  Test dc2rpy  ==========');
C = rpy2dc([1.3, -2.3, 0.7])
rpy = dc2rpy(C)

fprintf('==========  Test dc2eul  ==========');
C = eul2dc([1.3, -2.3, 0.7])
eul = dc2eul(C)

fprintf('==========  Test tr2diff  ==========');
AA1 = eul2dc([1.3, -2.3, 0.7]);
AA2 = eul2dc([2.1,  0.5, 4.7]);
diff = tr2diff(AA1,AA2)

fprintf('==========  Test aw  ==========');
w0 = [1.3, -2.3,  0.7];
E0 = aw(w0)

fprintf('==========  Test build_Inc_Mat  ==========');
BB = [0, 1, 2, 1];
[S0,SS,SE] = build_Inc_Mat(BB)
