clear all; clc;
%%%%%%%%%%%%%%%
% dr=0.15;
% dr=0.3;
% dr=0.6;
dr=1.2;
%%%%%%%%%%%%%%%
% dr=0.25;
% dr=0.5;
% dr=0.75;
% dr=1;
%%%%%%%%%%%%%%%
% L=1e-3;
% dr=0.229;
% dr=0.5;
% Wr=1.5;
% Lr=0.5;
%%%%%%%%%%%%%%%
 x=0:0.01:2;
%  y=[0.01:0.01:0.9,1.1:0.01:1.9,1.9:0.01:2];
 y=[0.01:0.01:0.9,1.1:0.01:2];
 %%%%%%%%%%%%%%%
[Lr Wr]=meshgrid(x,y);
a=Wr-1;
b=Wr-2;
c=1-dr.*Lr.*Wr+2*dr.*Lr;
d=log(-(-2+Wr)./Wr);
f=2*Lr-Wr-2*Lr.*Wr+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g5=-9*(Wr - 2)^2*(2*Lr - Wr - 2*Lr*Wr + 2)^2*log(-(-2+Wr)/Wr)^2;
% g6=3*(Wr - 2)^2*(2*Lr - Wr - 2*Lr*Wr + 2)*(6*Lr - 13*Wr - 6*Lr*Wr + 14)*log(-(-2+Wr)/Wr);
% g7=2*(Wr - 1)*(36*Lr^2*Wr^2 - 72*Lr^2*Wr + 36*Lr^2 - 8*Lr*Wr^3 + 96*Lr*Wr^2 - 192*Lr*Wr + 104*Lr - 16*Wr^3 + 103*Wr^2 - 188*Wr + 104);
% g8=((2*Lr - Wr - 2*Lr*Wr + 2)*log(-(-2+Wr)./Wr)+2*(Wr - 1))^2*36*(Wr-1)*(2-Wr)*dr*Lr;
% g9=((2*Lr - Wr - 2*Lr*Wr + 2)*log(-(-2+Wr)./Wr)+2*(Wr - 1))^2*36*(Wr-1)*((2-Wr)*dr*Lr+1);
N=a.*c.*(8*Lr-10*Wr-16*Lr.*Wr+8*Lr.*Wr.^2+4*Wr.^2-d.*f.^2+6).^2;
A=8*Lr.*Wr.*(140680*Wr-29800*Lr+62225*Lr.*Wr-67975*Lr.*Wr.^2+40825*Wr.^3.*Lr-...
    12725*Wr.^4.*Lr+1600*Wr.^5.*Lr-143585*Wr.^2+80750*Wr.^3-23671*Wr.^4+2816*Wr.^5-72146);
B=4800*Lr.*Wr.*dr.*(28*Wr.^5.*Lr.^2-196*Wr.^4.*Lr.^2+560*Wr.^3.*Lr.^2-840.*Lr.^2.*Wr.^2+700*Lr.^2.*Wr-308*Lr.^2+...
    36*Wr.^5.*Lr-270*Wr.^4.*Lr+828*Wr.^3.*Lr-1332*Lr.*Wr.^2+1188*Lr.*Wr-558*Lr+12*Wr.^5-96*Wr.^4+315*Wr.^3-543*Wr.^2+519*Wr-261);
C=2*(60624*Lr-175432*Wr+129600*dr.*Lr+259200*dr*Lr.^2+134400*dr*Lr.^3+23400*Lr.^2+...
    318590*Wr.^2-303660*Wr.^3+160135*Wr.^4-44249*Wr.^5+4992*Wr.^6+39624);
D=150.*d.^2.*b.*f.^2.*(8*Lr.*Wr.*dr.*(4*Lr.^2.*Wr.^2-12*Lr.^2.*Wr+12*Lr.^2+6*Lr.*Wr.^2-24*Lr.*Wr+...
    30*Lr+3*Wr.^2-15*Wr+24)+32*(-Lr.^3-3*Lr.^2-3*Lr).*dr+3*Wr.^3-18*Wr.^2+36*Wr-24);
E=15*d.*b.*f.*(320*Lr.*Wr.*dr.*(10*Wr.^3.*Lr.^2-40*Lr.^2.*Wr.^2+60*Lr.^2.*Wr-40*Lr.^2+15*Wr.^3.*Lr-72*Lr.*Wr.^2+126*Lr.*Wr-96*Lr+6*Wr.^3-...
    33*Wr.^2+66*Wr-57)+(3200*Lr.^3+8640*Lr.^2+5760*Lr).*dr+1200*Lr-4288*Wr-3640*Lr.*Wr+3980*Lr.*Wr.^2-1850*Wr.^3.*Lr+...
    310*Wr.^4.*Lr+4216*Wr.^2-1792*Wr.^3+279*Wr.^4+1584);
M=-(A+B+C+D+E)./(3600*N);
%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective_mass_for_tapering_cantilever_beam_3D
% figure(3)
% mesh(Lr,Wr,M);
% xlabel('\itw_r');
% ylabel('\itl_r');
% zlabel('\itM_{coeff}');
%%%%%%%%%%%%%%%%%%%%%%%%%
% Effective_mass_for_tapering_cantilever_beam_2D_contour
figure(4)
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(224)
%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,h]=contour(Lr,Wr,M);
% clabel(C,h)
%%%%%%%%%%%%%%%%%%%%%%%%%
[C,h]=contourf(Lr,Wr,M,9,'c-','LineWidth',1.5);
clabel(C,h,'Color','k','FontSize',12,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('长度比 \itl_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('宽度比 \itw_r','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
% title('厚度比 \itd_r=0.15','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
% title('厚度比 \itd_r=0.3','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
% title('厚度比 \itd_r=0.6','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
title('厚度比 \itd_r=1.2','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
grid on
box on
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')
hold on