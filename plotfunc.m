function plotfunc(FOM1,FOM2,FOM3,FOM4,FOM5,FOM6,ITER,PhanNum)
lw = 3;
ms = 10;
%% 
%--------------------------------------------------------------------------
figure;
semilogy(1:ITER,FOM1.NOFV,'-cs','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
semilogy(1:ITER,FOM2.NOFV,'-bv','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
semilogy(1:ITER,FOM3.NOFV,'-k+','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
semilogy(1:ITER,FOM4.NOFV,'-m*','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
semilogy(1:ITER,FOM5.NOFV,'-gd','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
semilogy(1:ITER,FOM6.NOFV,'-ro','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);
hold off
h=legend('PKMA','PPGA','APPGA(\omega=1/4)','APPGA(\omega=1/2)','APPGA(\omega=3/4)','APPGA(\omega=1)');
xlabel('Iteration number','FontWeight','bold');ylabel('NOFV','FontWeight','bold');
set(gca,'FontSize',35,'FontWeight','bold');
set(gcf,'position',[0 0 910 710],'Units','normalized'); 
set(h,'FontSize',24);
xticks([0:20:100]);
yticks([1e-6 1e-4 1e-2 1]);
%--------------------------------------------------------------------------
figure;
start_iter = 3; 
plot(start_iter:ITER,FOM1.NRMSE(start_iter:ITER),'-cs','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
plot(start_iter:ITER,FOM2.NRMSE(start_iter:ITER),'-bv','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms,'MarkerFaceColor','none');hold on
plot(start_iter:ITER,FOM3.NRMSE(start_iter:ITER),'-k+','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
plot(start_iter:ITER,FOM4.NRMSE(start_iter:ITER),'-m*','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
plot(start_iter:ITER,FOM5.NRMSE(start_iter:ITER),'-gd','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
plot(start_iter:ITER,FOM6.NRMSE(start_iter:ITER),'-ro','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);
hold off
h=legend({'PKMA','PPGA','APPGA(\omega=1/4)','APPGA(\omega=1/2)','APPGA(\omega=3/4)','APPGA(\omega=1)'});
xlabel('Iteration number','FontWeight','bold');ylabel('NRMSE','FontWeight','bold');
set(gca,'FontSize',35,'FontWeight','bold');
set(gcf,'position',[0 0 910 710],'Units','normalized'); 
set(h,'FontSize',28);
xticks([0:20:100]);
%--------------------------------------------------------------------------
figure;
plot(1:ITER,FOM1.PSNR,'-cs','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM2.PSNR,'-bv','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms,'MarkerFaceColor','none');hold on
plot(1:ITER,FOM3.PSNR,'-k+','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM4.PSNR,'-m*','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);hold on
plot(1:ITER,FOM5.PSNR,'-gd','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms);hold on
plot(1:ITER,FOM6.PSNR,'-ro','LineWidth',lw,'MarkerIndices',1:10:ITER,'MarkerSize',ms+2);
hold off
h=legend({'PKMA','PPGA','APPGA(\omega=1/4)','APPGA(\omega=1/2)','APPGA(\omega=3/4)','APPGA(\omega=1)'});
xlabel('Iteration number','FontWeight','bold');ylabel('PSNR','FontWeight','bold');
set(gca,'FontSize',35,'FontWeight','bold');
set(gcf,'position',[0 0 910 710],'Units','normalized'); 
set(h,'FontSize',28);
xticks([0:20:100]);
%--------------------------------------------------------------------------

end