
% Initialize plot properties
hf = figure(1);
clf
hold on
set(hf,'Renderer','painters');
ax = gca;
set(ax,'XDir','reverse');
ylim([-0.2 +1.3])
xlim([-0.9 +0.6])
%axis equal

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

% Position CoM, ZMP, and CoP
Y = CoM(1,1);
X = CoM(1,2);
hcm = plot(X,Y,'hk','MarkerFaceColor','k','MarkerSize',8);
%Y = -0.1;
%X = ZMP(1);
%hzmp = plot(X,Y,'hm','MarkerFaceColor','m','MarkerSize',8);
Y = -0.03;
X = CoP(1);
hcp = plot(X,Y,'hr','MarkerFaceColor','r','MarkerSize',8);
% Link centroids
Y = [ R0(1,1) RR(1,1:3:nq) ]; 
X = [ R0(1,2) RR(1,2:3:nq) ]; 
hce = plot(X,Y,'sr','MarkerSize',7);
% Joints/End effectors/Links
Y = [ R0(1,1) Pj(1,1) Pj(1,4) Pj(1,7) Pe(1,7) ]; 
X = [ R0(1,2) Pj(1,2) Pj(1,5) Pj(1,8) Pe(1,8) ]; 
hje1 = plot(X,Y,'ob-','MarkerSize',5,'MarkerFaceColor','b');
Y = [ Pj(1,10) Pe(1,10) ]; 
X = [ Pj(1,11) Pe(1,11) ]; 
hje2 = plot(X,Y,'ob-','MarkerSize',5,'MarkerFaceColor','b');
Y = Pj(1,1:3:nq); 
X = Pj(1,2:3:nq); 
hjo = plot(X,Y,'ob','MarkerSize',10);
% Contact points
Y = Pcon(1,1:3:nc);     % Markers
X = Pcon(1,2:3:nc); 
hco = plot(X,Y,'ob','MarkerSize',5,'MarkerFaceColor','b');
% Base to contact points
Y = [ Pcon(1,1) R0(1,1) Pcon(1,4) ]; 
X = [ Pcon(1,2) R0(1,2) Pcon(1,5) ]; 
hec = plot(X,Y,'b--','LineWidth',0.5);
% Contact forces
Fmax = max( [ Fcon(:,1); Fcon(:,4) ] );
Y = Pcon(1,1:3:nc)-0.03;    % Arrows tip
X = Pcon(1,2:3:nc); 
hfc0 = plot(X,Y,'^k','MarkerSize',4,'MarkerFaceColor','k');
Y = [ -0.03  -0.16*Fcon(1,1)/Fmax-0.03 ];    % Left arrow
X = [ Pcon(1,2) Pcon(1,2) ];
hfc1 = plot(X,Y,'-k','LineWidth',2);
Y = [ -0.03  -0.16*Fcon(1,4)/Fmax-0.03 ];    % Right arrow
X = [ Pcon(1,5) Pcon(1,5) ];
hfc2 = plot(X,Y,'-k','LineWidth',2);
% Plot impact plane
X = xlim;
Y = [ x_plane x_plane ];
plot(X,Y,':k','LineWidth',0.5);
title( [ 't = ' num2str( t(1),'%5.1f' ) ' s' ] );

% Plot any other elements
X = [ 0.0261  -0.1239  -0.1239  0.0261 ];
%X = [ 0.0261  -0.1239  -0.1239 -0.2739 ];
Y = [ 0.4666   0.3166   0.4666  0.4666 ];
plot(X,Y,':r','LineWidth',1);

% Plot any other text
X = Pe(1,11)+0.1;
Y = Pe(1,10);
htxt1 = text(X,Y,'Link 4');

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
  set(hcm,'XData',X);
  set(hcm,'YData',Y);
  %X = ZMP(i);
 % hzmp.XData = X;
  X = CoP(i);
  set(hcp,'XData',X);
  % Update link centroids
  Y = [ R0(i,1) RR(i,1:3:nq) ]; 
  X = [ R0(i,2) RR(i,2:3:nq) ]; 
  set(hce,'XData',X);
  set(hce,'YData',Y);
  % Update joints/end effectors/links
  Y = [ R0(i,1) Pj(i,1) Pj(i,4) Pj(i,7) Pe(i,7) ]; 
  X = [ R0(i,2) Pj(i,2) Pj(i,5) Pj(i,8) Pe(i,8) ]; 
  set(hje1,'XData',X);
  set(hje1,'YData',Y);
  Y = [ Pj(i,10) Pe(i,10) ]; 
  X = [ Pj(i,11) Pe(i,11) ]; 
  set(hje2,'XData',X);
  set(hje2,'YData',Y);
  Y = Pj(i,1:3:nq); 
  X = Pj(i,2:3:nq); 
  set(hjo,'XData',X);
  set(hjo,'YData',Y);
  % Update contact points
  Y = Pcon(i,1:3:nc); 
  X = Pcon(i,2:3:nc); 
  set(hco,'XData',X);
  set(hco,'YData',Y);
  Y = [ Pcon(i,1) R0(i,1) Pcon(i,4) ]; 
  X = [ Pcon(i,2) R0(i,2) Pcon(i,5) ]; 
  set(hec,'XData',X);
  set(hec,'YData',Y);
  % Update contact forces
  Y = Pcon(i,1:3:nc)-0.03; 
  X = Pcon(i,2:3:nc); 
  set(hfc0,'XData',X);
  set(hfc0,'YData',Y);
  Y = [ -0.03  -0.16*Fcon(i,1)/Fmax-0.03 ];
  X = [ Pcon(i,2) Pcon(i,2) ];
  set(hfc1,'XData',X);
  set(hfc1,'YData',Y);
  Y = [ -0.03  -0.16*Fcon(i,4)/Fmax-0.03 ];
  X = [ Pcon(i,5) Pcon(i,5) ];
  set(hfc2,'XData',X);
  set(hfc2,'YData',Y);
  % Update time in title
  title( [ 't = ' num2str( t(i),'%5.1f' ) ' s' ] );
  % Update any text
  X = Pe(i,11)+0.1;
  Y = Pe(i,10);
  set(htxt1,'Position',[X Y]);

  % Record frames for movie
  if ( flag_rec == 1 )
    idx = idx + 1;
    anim_movie(idx) = getframe(hf);
  end
  
end

clear hf ax nt nq nc flag_step flag_rec X Y hcm hzmp hcp hce ...
      hje hjo hco hec Fmax hfc0 hfc1 hfc2 i idx hje1 hje2

% End of script

