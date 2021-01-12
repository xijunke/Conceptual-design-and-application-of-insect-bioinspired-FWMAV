% Linearity_analysis_transmission_drive_mechanism
% Transmission_ratio_four_link_size
clear all;clc;
% Refer to��2011-ICIRS-System identification, modeling, and optimization of ...
% an insect-sized flapping-wing micro air vehicle-Finio_IROS11_c-8p
%% 
% L_c=300e-6;   % um����to����m:  e-6
% L_s1=500e-6;   % um����to����m:  e-6
% L_s2=300e-6;   % um����to����m:  e-6
% L_s3=630e-6;   % um����to����m:  e-6
% %  L_s3=L_s1;
%% 
L_c=312e-6;   % um����to����m:  e-6
L_s1=400e-6;   % um����to����m:  e-6
L_s2=291e-6;   % um����to����m:  e-6
L_s3=498e-6;   % um����to����m:  e-6
% L_s3=L_s1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% % lb = [50e-6,50e-6,50e-6,50e-6];          % Set lower bounds 
% % ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
% % x0= [L_c,L_s1,L_s2,L_s3];  
% optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; % objTdiff_linear % delta_pp=250e-6;
% L_c=optimal_four_link_size(1)*10^-6;   % m
% L_s1=optimal_four_link_size(2)*10^-6;   % m 
% L_s2=optimal_four_link_size(3)*10^-6;   % m
% L_s3=optimal_four_link_size(4)*10^-6;   % m
% % L_c=linspace(50e-6,1000e-6,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
delta=linspace(-250e-6,250e-6,1000);
% delta=linspace(-300e-6,300e-6,1000);
% delta=linspace(-350e-6,350e-6,1000);
% 
% delta_pp=250e-6;   % um��
% delta=delta_pp;     % 
%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% phi_l=(1/2)*pi-acos((1/2)*((L_c+L_s1-L_s3-delta).^2+2*L_s2^2+(L_s1-L_s3)^2-L_c^2)...
%           ./(sqrt(L_s2^2+(L_s1-L_s3)^2)*sqrt(L_s2^2+(L_c+L_s1-L_s3-delta).^2)))-...
%          atan(L_s2./(L_c+L_s1-L_s3-delta))-atan((L_s1-L_s3)./L_s2);  
phi_l=(1/2)*pi-acos((1/2)*((L_c+L_s1-L_s3-delta).^2+2*L_s2^2+(L_s1-L_s3)^2-L_c^2)...
          ./(sqrt(L_s2^2+(L_s1-L_s3)^2)*sqrt(L_s2^2+(L_c+L_s1-L_s3-delta).^2)))-...
         atan2(L_s2,(L_c+L_s1-L_s3-delta))-atan2((L_s1-L_s3),L_s2);
%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ��
% dphi_l_to_ddelta=(2*((1/2)*(-2*L_c-2*L_s1+2*L_s3+2*delta)./...
%     (sqrt(L_s2^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
%     sqrt(L_s2^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2))-...
%     (1/4)*((L_c+L_s1-L_s3-delta).^2+2*L_s2^2+(L_s1-L_s3)^2-L_c^2).*(-2*L_c-2*L_s1+2*L_s3+2*delta)...
%     ./(sqrt(L_s2^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
%     (L_s2^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2).^(3/2))))./...
%     sqrt(4-((L_c+L_s1-L_s3-delta).^2+2*L_s2^2+(L_s1-L_s3)^2-L_c^2).^2./...
%     ((L_s2^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
%     (L_s2^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2)))...
%     -L_s2./((L_c+L_s1-L_s3-delta).^2.*(1+L_s2^2./(L_c+L_s1-L_s3-delta).^2));
%%%%%%%%%%%%%%%%%%
% A0=L_c+L_s1-L_s3-delta;
% % A=-2*L_c-2*L_s1+2*L_s3+2*delta;
% % A=-2*(L_c+L_s1-L_s3-delta);
% A=-2*A0;  % size(A)
% B=L_s2^2+L_s1^2-2*L_s1*L_s3+L_s3^2; % size(B) % (1*1)
% C=L_s2^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2;
% size(C)
% D=(A0).^2+2*L_s2^2+(L_s1-L_s3)^2-L_c^2;
% size(D)
% E=0.5*(A)./(sqrt(B).*sqrt(C));
% size(E)
% F=0.25*D.*A./(sqrt(B).*C.^(3/2));
% size(F)
% G=sqrt(4-D.^2./(B.*C));
% size(G)
% H=(A0.^2.*(1+L_s2^2./(A0).^2));
% size(H)
%%%%%%%%%%%%%%%%%
% dphi_l_to_ddelta=(2*(E-F))./G-L_s2./H;
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
figure(1)
T_est=1/L_s2;
phi_l_est=T_est*delta;
plot(delta*10^6,-phi_l*180/pi,'g-',delta*10^6,phi_l_est*180/pi,'b-.','LineWidth',3)
hold on
plot([0,0],[-60,60],'k-','LineWidth',1.5)
hold on
plot([-250,250],[0,0],'k-','LineWidth',1.5)
xlabel('ѹ�������������λ�� \rm(\it\delta \rm(t) (um))','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
ylabel('�����Ĵ�� (\it\phi_{l,w} \rm(��))','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
legend('���ŵ����Զ�','�������Զ�')
% axis([-260,260,-inf,inf])
axis([-250,250,-60,60])
% axis([-360,360,-70,70])
grid on
box on
hold on
set(gca,'XTick',(-250:50:250),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')

% 
figure(2)
T_nonlinear=-phi_l./delta;  % transmission ratio
plot(delta*10^6,T_nonlinear,'g-','LineWidth',3)
xlabel('ѹ��������������λ�� (\it\delta \rm(t) (um))','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')
ylabel('����ѧ���������Ĵ����� (\itT_{ratio} \rm(rad/m))','Color','k','FontSize',16,'FontName','Times','FontWeight','Bold')  
legend('������ vs ����������λ��')
grid on
box on
hold on
set(gca,'XTick',(-250:50:250),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
set(gca,'Fontsize',14,'FontName','Times','FontWeight','Bold','Ycolor','k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%