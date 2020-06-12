
plot(t,CoM(:,2),'r',t,CoP(:,1),'b',t,Pcon(:,2),'k--',t,Pcon(:, 5),'k--', ...
                    t,CoP(:,2),'b',t,Pcon(:,8),'k--',t,Pcon(:,11),'k--')
legend( { ' CoM',' CoP',' P_{con}' },'FontSize',14)
text( 0.01, (Pcon(1,2)+Pcon(1, 5))/2, 'Left foot',  'FontSize',12, 'Color','b' )
text( 0.01, (Pcon(1,8)+Pcon(1,11))/2, 'Right foot', 'FontSize',12, 'Color','b' )
grid on
xlabel('Time [s]','FontSize',18)
ylabel('Horizontal component [m]','FontSize',18)
set(gcf,'color','w')