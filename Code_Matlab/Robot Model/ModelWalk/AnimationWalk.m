
% Initialize plot properties
hf = figure(1);
clf
hold on
hf.Renderer = 'painters';
ax = gca;
ax.XDir = 'reverse';
ylim([-0.2 2.0])
axis equal

% Parameters
nt = size(t,1);   % Number of instants
nq = 3*num_q;     % Number of joints
nc = 3*npc;       % Number of contact points

% If not defined use step = 1
flag_step = exist('step','var');
if ( flag_step == 0 )
  step = 1;
end
% If not defined don't record frames
flag_rec = exist('rec','var');
if ( flag_rec == 1 )
  clear anim_movie
end

% Position CoM, XCoM, and CoP
Y = CoM(1,1);
X = CoM(1,2);
hcm = plot(X,Y,'hk','MarkerFaceColor','k','MarkerSize',8);
Y = -0.03;
X = CoP(1,1);
hcp1 = plot(X,Y,'hr','MarkerFaceColor','r','MarkerSize',8);
Y = -0.03;
X = CoP(1,2);
hcp2 = plot(X,Y,'hr','MarkerFaceColor','r','MarkerSize',8);
% Link centroids
Y = [ R0(1,1) RR(1,1:3:nq) ]; 
X = [ R0(1,2) RR(1,2:3:nq) ]; 
hce = plot(X,Y,'sr','MarkerSize',7);
% Links 3-2-1-0-11
Y = [ RR(1,7) Pj(1,7) Pj(1,4) Pj(1,1) Pj(1,31) Pe(1,31) ]; 
X = [ RR(1,8) Pj(1,8) Pj(1,5) Pj(1,2) Pj(1,32) Pe(1,32) ]; 
hje1 = plot(X,Y,'ob-','MarkerSize',5,'MarkerFaceColor','b');
% Links 6-5-4
Y = [ RR(1,16) Pj(1,16) Pj(1,13) Pj(1,10) ]; 
X = [ RR(1,17) Pj(1,17) Pj(1,14) Pj(1,11) ]; 
hje2 = plot(X,Y,'ob-','MarkerSize',5,'MarkerFaceColor','b');
% Links 8-7-9-10
Y = [ Pe(1,22) Pj(1,22) Pj(1,19) Pj(1,28) Pe(1,28) ]; 
X = [ Pe(1,23) Pj(1,23) Pj(1,20) Pj(1,29) Pe(1,29) ]; 
hje3 = plot(X,Y,'ob-','MarkerSize',5,'MarkerFaceColor','b');
% Joints
Y = Pj(1,1:3:nq); 
X = Pj(1,2:3:nq); 
hjo = plot(X,Y,'ob','MarkerSize',10);
% Contact points
Y = Pcon(1,1:3:nc);     % Markers
X = Pcon(1,2:3:nc); 
hco = plot(X,Y,'ob','MarkerSize',5,'MarkerFaceColor','b');
Y = [ Pcon(1,1) RR(1,7) Pcon(1,4) ]; 
X = [ Pcon(1,2) RR(1,8) Pcon(1,5) ]; 
hec1 = plot(X,Y,'b--','LineWidth',0.5);
Y = [ Pcon(1,7) RR(1,16) Pcon(1,10) ]; 
X = [ Pcon(1,8) RR(1,17) Pcon(1,11) ]; 
hec2 = plot(X,Y,'b--','LineWidth',0.5);
% Contact forces
Fmax = max( [ Fcon(:,1); Fcon(:,4); Fcon(:,7); Fcon(:,10) ] );
Y = Pcon(1,1:3:nc)-0.03;    % Arrows tip
X = Pcon(1,2:3:nc); 
hfc0 = plot(X,Y,'^k','MarkerSize',4,'MarkerFaceColor','k');
Y = [ -0.03  -0.16*Fcon(1,1)/Fmax-0.03 ];    % Left arrow P1
X = [ Pcon(1,2) Pcon(1,2) ];
hfc1 = plot(X,Y,'-k','LineWidth',2);
Y = [ -0.03  -0.16*Fcon(1,4)/Fmax-0.03 ];    % Right arrow P1
X = [ Pcon(1,5) Pcon(1,5) ];
hfc2 = plot(X,Y,'-k','LineWidth',2);
Y = [ -0.03  -0.16*Fcon(1,7)/Fmax-0.03 ];    % Left arrow P2
X = [ Pcon(1,8) Pcon(1,8) ];
hfc3 = plot(X,Y,'-k','LineWidth',2);
Y = [ -0.03  -0.16*Fcon(1,10)/Fmax-0.03 ];    % Right arrow P2
X = [ Pcon(1,11) Pcon(1,11) ];
hfc4 = plot(X,Y,'-k','LineWidth',2);

% Plot impact plane
X = xlim;
Y = [ x_plane x_plane ];
plot(X,Y,':k','LineWidth',0.5);
% Title
title( [ 't = ' num2str( t(1),'%5.1f' ) ' s' ] );
% Text
fact = height/40;
tLL = text(Pj(1,5)-fact,Pj(1,4),'Left','FontSize',12,'Color','r');
tRL = text(Pj(1,14)-fact,Pj(1,13),'Right','FontSize',12,'Color','r');
tLA = text(Pj(1,23)-fact,Pj(1,22),'Left','FontSize',12,'Color','r');
tRA = text(Pj(1,29)-fact,Pj(1,28),'Right','FontSize',12,'Color','r');

% Record frames for movie
if ( flag_rec == 1 )
  idx = 1;
  anim_movie(idx) = getframe(hf);
end

for i = 2:step:nt
  pause(0.01)
  % Update CoM, XCoM, and CoP
  Y = CoM(i,1);
  X = CoM(i,2);
  hcm.XData = X;
  hcm.YData = Y;
  X = CoP(i,1);
  hcp1.XData = X;
  X = CoP(i,2);
  hcp2.XData = X;
  % Update link centroids
  Y = [ R0(i,1) RR(i,1:3:nq) ]; 
  X = [ R0(i,2) RR(i,2:3:nq) ]; 
  hce.XData = X;
  hce.YData = Y;
  % Links 3-2-1-0-11
  Y = [ RR(i,7) Pj(i,7) Pj(i,4) Pj(i,1) Pj(i,31) Pe(i,31) ]; 
  X = [ RR(i,8) Pj(i,8) Pj(i,5) Pj(i,2) Pj(i,32) Pe(i,32) ]; 
  hje1.XData = X;
  hje1.YData = Y;
  % Links 6-5-4
  Y = [ RR(i,16) Pj(i,16) Pj(i,13) Pj(i,10) ]; 
  X = [ RR(i,17) Pj(i,17) Pj(i,14) Pj(i,11) ]; 
  hje2.XData = X;
  hje2.YData = Y;
  % Links 8-7-9-10
  Y = [ Pe(i,22) Pj(i,22) Pj(i,19) Pj(i,28) Pe(i,28) ]; 
  X = [ Pe(i,23) Pj(i,23) Pj(i,20) Pj(i,29) Pe(i,29) ]; 
  hje3.XData = X;
  hje3.YData = Y;
  % Joints
  Y = Pj(i,1:3:nq); 
  X = Pj(i,2:3:nq); 
  hjo.XData = X;
  hjo.YData = Y;
  % Update contact points
  Y = Pcon(i,1:3:nc); 
  X = Pcon(i,2:3:nc); 
  hco.XData = X;
  hco.YData = Y;
  Y = [ Pcon(i,1) RR(i,7) Pcon(i,4) ]; 
  X = [ Pcon(i,2) RR(i,8) Pcon(i,5) ]; 
  hec1.XData = X;
  hec1.YData = Y;
  Y = [ Pcon(i,7) RR(i,16) Pcon(i,10) ]; 
  X = [ Pcon(i,8) RR(i,17) Pcon(i,11) ]; 
  hec2.XData = X;
  hec2.YData = Y;
  % Update contact forces
  Y = Pcon(i,1:3:nc)-0.03; 
  X = Pcon(i,2:3:nc); 
  hfc0.XData = X;
  hfc0.YData = Y;
  Y = [ -0.03  -0.16*Fcon(i,1)/Fmax-0.03 ];
  X = [ Pcon(i,2) Pcon(i,2) ];
  hfc1.XData = X;
  hfc1.YData = Y;
  Y = [ -0.03  -0.16*Fcon(i,4)/Fmax-0.03 ];
  X = [ Pcon(i,5) Pcon(i,5) ];
  hfc2.XData = X;
  hfc2.YData = Y;
  Y = [ -0.03  -0.16*Fcon(i,7)/Fmax-0.03 ];
  X = [ Pcon(i,8) Pcon(i,8) ];
  hfc3.XData = X;
  hfc3.YData = Y;
  Y = [ -0.03  -0.16*Fcon(i,10)/Fmax-0.03 ];
  X = [ Pcon(i,11) Pcon(i,11) ];
  hfc4.XData = X;
  hfc4.YData = Y;
 
  % Update time in title and text positions
  title( [ 't = ' num2str( t(i),'%5.1f' ) ' s' ] );
  tLL.Position = [ Pj(i,5)-fact  Pj(i,4) ];
  tRL.Position = [ Pj(i,14)-fact  Pj(i,13) ];
  tLA.Position = [ Pj(i,23)-fact  Pj(i,22) ];
  tRA.Position = [ Pj(i,29)-fact  Pj(i,28) ];
 
  % Record frames for movie
  if ( flag_rec == 1 )
    idx = idx + 1;
    anim_movie(idx) = getframe(hf);
  end
  
end

clear hf ax nt nq nc flag_step flag_rec X Y hcm hxcm hcp hce hje1 hje2 ...
      hje3 hjo hco hec1 hec2 Fmax hfc0 hfc1 hfc2 hfc2 hfc3 i idx hcp1 ...
      hcp2 hfc4

% End of script

