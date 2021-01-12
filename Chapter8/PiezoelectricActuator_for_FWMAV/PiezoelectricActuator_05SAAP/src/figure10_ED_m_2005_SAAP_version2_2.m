clear all; clc;
syms E v E1 E2 v12 v21 G12  z1 z2 z3 z0 E3 L W d31 d32 y p1 p2 tc;
% t=0:0.0001:0.025;
% E3=100*sin(2*pi*120*t)+100;     
% plot(t,E3)
% PZT
p1=7800;   %kg/m
t=127e-6;
E3=300/t;            %v
E=6.2e10;   v=0.31; %Gpa &G12=24e9
d31=-320e-12; d32=-320e-12; %pm/v
d=[d31;d32;0];
% CF
p2=1500;   %kg/m
E1=350e9;  E2=7e9; G12=5e9;   %Gpa
v12=0.33; v21=E2*v12/E1;         %
%  tr=tc/t;  CF
l_r=1;w_r=1.5;

%% PZT
Q112=E/(1-v^2);  Q222=Q112;
Q122=E*v/(1-v^2);
Q662=E/(2*(1+v)); G=Q662;
%CF
Q11=E1/(1-v12*v21);
Q12=v12*E2/(1-v12*v21);
Q22=E2/(1-v12*v21);
Q66=G12;
%
%%  
Du=[];%
TD=[];%tip displacement
BF=[];%block force

% PZT layer
Q0=[Q112 Q122 0
      Q122 Q222 0
      0      0    Q662];
% CF layer
Q1=[Q11 Q12 0
      Q12 Q22 0
      0      0    Q66];
% 
Qmid=Q1;
Qfirst=Q0;
Qthird=Q0;
%%
for tr=0:0.01:2;
    tc=tr*t;
   % 
   wnom=1e-3; l=5e-3;    %
   % m
   z0=-tc/2-t; 
   z1=-tc/2; 
   z2=tc/2; 
   z3=tc/2+t;     
   
   % 
   g_delta=(1+2*l_r);
   g_c=8*(1-w_r).^3;
   g_d=-6*(w_r-1).*(3+4*l_r-2*w_r-4*l_r.*w_r);
   g_e=3*(-2-2*l_r+w_r+2*l_r.*w_r).^2.*log((2-w_r)./w_r);
   G_le=(g_d+g_e)./g_c;
   G_F=g_delta/G_le;
   G_U=g_delta*G_F;
   
   % 
   A=Qfirst*(z1-z0)+Qmid*(z2-z1)-Qthird*(z3-z2);
   B=(Qfirst*(z1^2-z0^2)+Qmid*(z2^2-z1^2)+Qthird*(z3^2-z2^2))*(1/2);
   D=(Qfirst*(z1^3-z0^3)+Qmid*(z2^3-z1^3)+Qthird*(z3^3-z2^3))*(1/3);
   F=[A,B;B,D];
   C=pinv(F); % pinv_Moore-Penrose pseudoinverse of matrix
   
  % 
   Fp=E3*(z1-z0)*Qfirst*d+E3*(z3-z2)*Qthird*d;  
   Mp=E3*1/2*(z1^2-z0^2)*Qfirst*d-E3*1/2*(z3^2-z2^2)*Qthird*d;
   Np=[Fp;Mp];
   PE3=C(4,1)*Np(1)+C(4,2)*Np(2)+C(4,4)*Np(4)+C(4,5)*Np(5);
  % du=3/8*PE3*PE3/F(4,4)/(2*p1*t+p2*tc)*G_U;
   du=-3/8*PE3*PE3/C(4,4)/(2*p1*t+p2*tc)*G_U;
   Du=[Du du];
   
   % displacement output
   td=PE3*l^2/2*g_delta;
   td=td*10^6;%um
   TD=[TD td];
   
   % blocked force output
   bf=3*PE3*wnom/(2*C(4,4)*l)*G_F;
   bf=bf*10^3;%mN
   BF=[BF bf];
end
%%
tr=0:0.01:2;
% Energy density output
figure(1);
plot(tr,Du,'k-','LineWidth',2);
% xlabel('Passive Layer Thickness Ratio ( \itt_r )');
% ylabel('Energy density \itED_m (J/kg)');
xlabel('中间被动结构碳纤维层与压电陶瓷层的厚度比 ( \itt_r )','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
ylabel('能量密度 \itED_m (J/kg)','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold');
grid on
box on
set(gca,'XTick',(0:0.1:2),'LineStyle','-','LineWidth',1.5)
set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')
axis(gca,[0,2,0,2.6])


% % displacement output
% figure(2);
% plot(tr,TD,'r-','LineWidth',2)
% % xlabel('Passive Layer Thickness Ratio ( \itt_r )');
% % ylabel('Tip Displacement \it\delta_{pp} (um)');
% xlabel('中间被动结构碳纤维层与压电陶瓷层的厚度比 ( \itt_r )');
% ylabel('压电驱动器远端尖部位移输出 \it\delta_{pp} (um)');
% grid on

% % blocked force output
% figure(3);
% plot(tr,BF,'b-','LineWidth',2)
% % xlabel('Passive Layer Thickness Ratio ( \itt_r )');
% % ylabel('Block Force \itF_b (mN)');
% xlabel('中间被动结构碳纤维层与压电陶瓷层的厚度比 ( \itt_r )');
% ylabel('压电驱动器远端尖部堵死力输出 \itF_b (mN)');
% grid on