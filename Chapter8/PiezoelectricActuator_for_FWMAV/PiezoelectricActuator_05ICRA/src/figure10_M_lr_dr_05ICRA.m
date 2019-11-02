%x=0.04:0.05:1.99;y=0.09:0.1:3.99;[w_r,l_r]=meshgrid(x,y)； g_delta=(1+2*l_r); g_c=8*(1-w_r).^3;g_d=6*(w_r-1).*(3+4*l_r-2*w_r-4*l_r.*w_r);
% g_e=3*(-2-2*l_r+w_r+2*l_r.*w_r).^2.*log((2-w_r)./w_r); G_le=(g_d+g_e)./g_c;G_F=g_delta./G_le; G_U=g_delta.*G_F; surf(w_r,l_r,G_U)
clear all; clc;
x=linspace(0,2,100);
y=linspace(0,2,100);
[lr,dr]=meshgrid(x,y);
M=(33+7*lr.*(13+9*lr+5*dr.*(4+3*lr.*(6+lr.*(11+2*lr.*(5+2*lr))))))./(140.*(1+dr.*lr).*(1+3*lr.*(1+lr)).^2);

% Effective_mass_for_rectangular_cantilever_beam_3D
figure(1)
mesh(lr,dr,M);
xlabel('长度比 \itl_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('厚度比 \itd_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
zlabel('有效质量系数 \itM(l_r,d_r)','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
grid on
box on
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')

% Effective_mass_for_rectangular_cantilever_beam_2D_contour
figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contour(lr,dr,M,'LineWidth',1.5);
% clabel(C,h)
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contourf(lr,dr,M,'c-','LineWidth',1.5);
% colormap autumn
%%%%%%%%%%%%%%%%%%%%%%%%%
[C,h]=contourf(lr,dr,M,9,'c-','LineWidth',1.5);
clabel(C,h,'Color','k','FontSize',12,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('长度比 \itl_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('厚度比 \itd_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
grid on
box on
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')