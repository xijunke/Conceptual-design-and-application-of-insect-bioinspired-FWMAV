% 手动输入线标
clear all;clc;
x=0.01:0.001:1.999; % 更高精度显示
y=0:0.01:4;
[wr lr]=meshgrid(x,y);
gc=8*(1-wr).^3;
gd=6*(wr-1).*(3+4*lr-2*wr-4*lr.*wr);
ge=3*(-2-2*lr+wr+2*lr.*wr).^2.*log((2-wr)./wr);
U=(1+2*lr).^2.*gc./(gd+ge);

% Energy_improvement_geometry_factor_3D
figure(3)
mesh(wr,lr,U);
% xlabel('\itw_r');
% ylabel('\itl_r');
% zlabel('\itG_U');
xlabel('宽度比 \itw_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('长度比 \itl_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
zlabel('能量改善几何因子 \itG_U','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
grid on
box on
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')

% Energy_improvement_geometry_factor_2D_contour
figure(4)
%%%%%%%%%%%%%%%%%%%%%%%%%
[C,h]=contour(wr,lr,U,'LineWidth',2.0);
clabel(C,h)
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contour(wr,lr,U,[0.4 0.6 0.8 1.0 1.2 1.3 1.33])
% clabel(C,h,'manual')
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contourf(wr,lr,U,'c-','LineWidth',2.0);
% colormap autumn
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contourf(wr,lr,U,9,'c-','LineWidth',2);
% clabel(C,h,'Color','k','FontSize',12,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%
% xlabel('\itw_r');
% ylabel('\itl_r');
xlabel('宽度比 \itw_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('长度比 \itl_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
grid on
box on
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')