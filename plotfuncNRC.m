function plotfuncNRC(FOM1,FOM2,FOM3,FOM4,FOM5,FOM6,ITER,PhanNum)
lw = 3;
ms = 12;
%--------------------------------------------------------------------------
figure;
plot(1:ITER,FOM1.NRC1,'-cs','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM2.NRC1,'-bv','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM3.NRC1,'-k+','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM4.NRC1,'-m*','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM5.NRC1,'-gd','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM6.NRC1,'-ro','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);
hold off
h=legend('PKMA','PPGA','APPGA(\omega=1/4)','APPGA(\omega=1/2)','APPGA(\omega=3/4)','APPGA(\omega=1)');
xlabel('Iteration number','FontWeight','bold');ylabel('NRC','FontWeight','bold');
set(gca,'FontSize',35,'FontWeight','bold');
set(gcf,'position',[0 50 910 710],'Units','normalized'); 
set(h,'FontSize',28);
yticks([0:0.1:0.6]);
xticks([0:100:500]);
xlim([0, 500]);
title('Smallest hot sphere','FontWeight','bold');
%--------------------------------------------------------------------------
figure;
plot(1:ITER,FOM1.NRC2,'-cs','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM2.NRC2,'-bv','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM3.NRC2,'-k+','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM4.NRC2,'-m*','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM5.NRC2,'-gd','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM6.NRC2,'-ro','LineWidth',lw,'MarkerIndices',1:50:ITER,'MarkerSize',ms+2);
hold off
h=legend('PKMA','PPGA','APPGA(\omega=1/4)','APPGA(\omega=1/2)','APPGA(\omega=3/4)','APPGA(\omega=1)');
xlabel('Iteration number','FontWeight','bold');ylabel('NRC','FontWeight','bold');
set(gca,'FontSize',35,'FontWeight','bold');
set(gcf,'position',[0 50 910 710],'Units','normalized'); 
set(h,'FontSize',28);
yticks([0:0.2:1]);
xticks([0:100:500]);
xlim([0, 500]);
title('Largest hot sphere','FontWeight','bold');
end