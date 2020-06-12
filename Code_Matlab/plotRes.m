
plot(t,CoM(:,2),'r',t,CoP(:,1),'b',t,Pcon(:,2),'k--',t,Pcon(:,5),'k--')
legend({' CoM',' CoP',' P_C_1 & P_C_2'},'FontSize',14)
grid on
xlabel('Time [s]','FontSize',18)
ylabel('Horizontal component [m]','FontSize',18)
set(gcf,'color','w')