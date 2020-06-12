
num_q = length(BB);
npc = size(Pc,2);

% Position data
fid = fopen(strcat(prefix,'.pos'),'r');
ncol = 12 + 13*num_q;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
t = Dmat(:,1);    
R0 = Dmat(:,2:4); 
Q0 = Dmat(:,5:7); 
is = 8;
ie = is+num_q-1;
q = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
RR = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Qi = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
CoM = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Pe = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Pj = Dmat(:,is:ie);
is = ie+1;
ie = is + 1;
CoP = Dmat(:,is:ie);

% Velocity data
fid = fopen(strcat(prefix,'.vel'),'r');
ncol = 10 + 10*num_q;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
v0 = Dmat(:,2:4); 
w0 = Dmat(:,5:7); 
is = 8;
ie = is+num_q-1;
qd = Dmat(:,is:ie); 
is = ie+1;
ie = is+3*num_q-1;
vv = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
ww = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
CoMv = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
vve = Dmat(:,is:ie);

% Acceleration data
fid = fopen(strcat(prefix,'.acc'),'r');
ncol = 10 + 7*num_q;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
vd0 = Dmat(:,2:4); 
wd0 = Dmat(:,5:7); 
is = 8;
ie = is+num_q-1;
qdd = Dmat(:,is:ie); 
is = ie+1;
ie = is+3*num_q-1;
vd = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
wd = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
CoMa = Dmat(:,is:ie);

% Force data
fid = fopen(strcat(prefix,'.for'),'r');
ncol = 13 + 19*num_q;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
F0 = Dmat(:,2:4); 
T0 = Dmat(:,5:7); 
is = 8;
ie = is+num_q-1;
tau = Dmat(:,is:ie); 
is = ie+1;
ie = is+3*num_q-1;
Fe = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Te = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
FF0 = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
TT0 = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
FF = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
TT = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Fj = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
Tj = Dmat(:,is:ie);

% Contact data
fid = fopen(strcat(prefix,'.con'),'r');
ncol = 1 + 9*npc;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
is = 2;
ie = is+3*npc-1;
Pcon = Dmat(:,is:ie); 
is = ie+1;
ie = is+3*npc-1;
Fcon = Dmat(:,is:ie); 
is = ie+1;
ie = is+npc-1;
delta = Dmat(:,is:ie); 
is = ie+1;
ie = is+npc-1;
delta1 = Dmat(:,is:ie); 
is = ie+1;
ie = is+npc-1;
delta1t = Dmat(:,is:ie); 

% Energy data
fid = fopen(strcat(prefix,'.ene'),'r');
ncol = 10 + 10*num_q;
Dmat = fscanf(fid,'%g',[ncol inf]);
Dmat = Dmat';
fclose(fid);
TK0 = Dmat(:,2);
is = 3;
ie = is+num_q-1;
TKj = Dmat(:,is:ie);
is = ie+1;
ie = is;
VG0 = Dmat(:,is:ie);
is = ie+1;
ie = is+num_q-1;
VGj = Dmat(:,is:ie);
is = ie+1;
ie = is;
PF0 = Dmat(:,is:ie);
WF0 = cumtrapz(t,PF0);
is = ie+1;
ie = is+num_q-1;
PFe = Dmat(:,is:ie);
WFe = cumtrapz(t,PFe);
is = ie+1;
ie = is+num_q-1;
Ptau = Dmat(:,is:ie);
Wtau = cumtrapz(t,Ptau);
is = ie+1;
ie = is+2;
LM0 = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
LMj = Dmat(:,is:ie);
is = ie+1;
ie = is+2;
HM0 = Dmat(:,is:ie);
is = ie+1;
ie = is+3*num_q-1;
HMj = Dmat(:,is:ie);

% Total energy and work
Etot = TK0 + sum(TKj,2) + VG0 + sum(VGj,2);
Wtot = WF0 + sum(WFe,2) + sum(Wtau,2);
error = 100*max(abs(Etot-Etot(1)-Wtot))/max(abs(Etot));
fprintf('\n Error: %4.2g%% \n\n',error);

% Total linear/angular momentum and linear impulse
% Angular momentum is wrt to CoM
nt = length(t);
IF0 = cumtrapz(t,F0);
IFE = cumtrapz(t,Fe);
Fg = repmat( (m0+sum(m))*Gravity', nt, 1 );
IFG = cumtrapz(t,Fg);
Ltot = zeros(nt,3);
Htot = zeros(nt,3);
IFtot = zeros(nt,3);
for i = 1:3
  Ltot(:,i) = LM0(:,i) + sum(LMj(:,i:3:3*num_q),2);
  Htot(:,i) = HM0(:,i) + sum(HMj(:,i:3:3*num_q),2);
  IFtot(:,i) = IF0(:,i) + sum(IFE(:,i:3:3*num_q),2) + IFG(:,i);
end

clear Dmat fid is ie nt i Fg
