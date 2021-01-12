% Sensitivity_analysis_transmission_ratio_difference_link_length
clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
delta_pp=250e-6;   % um——to——m:  e-6——delta_pp=(-250,250)um
% delta_pp=300e-6;
delta=delta_pp;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% lb = [50e-6,50e-6,50e-6,50e-6];          % Set lower bounds 
% ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
% x0= [L_c,L_s1,L_s2,L_s3];  
optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; % objTdiff_linear % delta_pp=250e-6;
% L_c=optimal_four_link_size(1)*10^-6;   % m
L_s1=optimal_four_link_size(2)*10^-6;   % m 
L_s2=optimal_four_link_size(3)*10^-6;   % m
L_s3=optimal_four_link_size(4)*10^-6;   % m
L_c=linspace(50e-6,1000e-6,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%：
dphi_l_to_ddelta=(2*((1/2)*(-2*L_c-2*L_s1+2*L_s3+2*delta)./...
    (sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    sqrt(L_s2.^2+L_c.^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2))-...
    (1/4)*((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3)^2-L_c.^2).*(-2*L_c-2*L_s1+2*L_s3+2*delta)...
    ./(sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c.^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2).^(3/2))))./...
    sqrt(4-((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3)^2-L_c.^2).^2./...
    ((L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c.^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2)))...
    -L_s2./((L_c+L_s1-L_s3-delta).^2.*(1+L_s2.^2./(L_c+L_s1-L_s3-delta).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%
T_est=1./L_s2;
obj_Tdiff=(abs(dphi_l_to_ddelta)-T_est).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity_analysis_transmission_ratio_difference_L_c
figure(3)
plot(L_c*10^6,obj_Tdiff,'k-','LineWidth',2)
hold on
plot([445.8109,445.8109],[min(obj_Tdiff)-0.25*10^-23,max(obj_Tdiff)+0.25*10^-23],'r-.','LineWidth',3)
xlabel('连杆(\itL_{c}\rm) (um)','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
ylabel('\itf(L_{s2})=\rm(\itd\phi_{theor}/d\delta|_{\delta_{pp}\rm=250um} -\itT_{est}\rm)^2 ','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
% title('在压电驱动器峰值位移下的理论传动比和线性传动比的差的平方和','Color','k','FontSize',20,'FontName','Times','FontWeight','Bold')
grid on
box on
set(gca,'XTick',(45:50:1005),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
axis(gca,[45,1005,min(obj_Tdiff)-0.25*10^-23,max(obj_Tdiff)+0.25*10^-23])
set(gca,'Fontsize',19,'FontName','Times','FontWeight','Bold','Ycolor','k')

%% (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% lb = [50e-6,50e-6,50e-6,50e-6];          % Set lower bounds 
% ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
% x0= [L_c,L_s1,L_s2,L_s3];  
optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; % objTdiff_linear % delta_pp=250e-6;
L_c=optimal_four_link_size(1)*10^-6;   % m
% L_s1=optimal_four_link_size(2)*10^-6;   % m 
L_s2=optimal_four_link_size(3)*10^-6;   % m
L_s3=optimal_four_link_size(4)*10^-6;   % m
L_s1=linspace(50e-6,1000e-6,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%：
dphi_l_to_ddelta=(2*((1/2)*(-2*L_c-2*L_s1+2*L_s3+2*delta)./...
    (sqrt(L_s2.^2+L_s1.^2-2*L_s1*L_s3+L_s3^2).*...
    sqrt(L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1.^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2))-...
    (1/4)*((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3).^2-L_c^2).*(-2*L_c-2*L_s1+2*L_s3+2*delta)...
    ./(sqrt(L_s2.^2+L_s1.^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1.^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2).^(3/2))))./...
    sqrt(4-((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3).^2-L_c^2).^2./...
    ((L_s2.^2+L_s1.^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1.^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2)))...
    -L_s2./((L_c+L_s1-L_s3-delta).^2.*(1+L_s2.^2./(L_c+L_s1-L_s3-delta).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%
T_est=1./L_s2;
obj_Tdiff=(abs(dphi_l_to_ddelta)-T_est).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity_analysis_transmission_ratio_difference_L_s1
figure(4)
plot(L_s1*10^6,obj_Tdiff,'k-','LineWidth',2)
hold on
plot([564.1168,564.1168],[min(obj_Tdiff)-0.25*10^6,max(obj_Tdiff)+0.25*10^6],'r-.','LineWidth',3)
xlabel('翅前缘杆根部前段杆(\itL_{s1}\rm) (um)','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
ylabel('\itf(L_{s2})=\rm(\itd\phi_{theor}/d\delta|_{\delta_{pp}\rm=250um} -\itT_{est}\rm)^2 ','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
% title('在压电驱动器峰值位移下的理论传动比和线性传动比的差的平方和','Color','k','FontSize',20,'FontName','Times','FontWeight','Bold')
grid on
box on
set(gca,'XTick',(45:50:1005),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
axis(gca,[45,1005,min(obj_Tdiff)-0.25*10^6,max(obj_Tdiff)+0.25*10^6])
set(gca,'Fontsize',19,'FontName','Times','FontWeight','Bold','Ycolor','k')

%% (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% x0= [L_c,L_s1,L_s2,L_s3];  
optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; % objTdiff_linear % delta_pp=250e-6;
L_c=optimal_four_link_size(1)*10^-6;   % m
L_s1=optimal_four_link_size(2)*10^-6;   % m 
% L_s2=optimal_four_link_size(3)*10^-6;   % m
L_s3=optimal_four_link_size(4)*10^-6;   % m
L_s2=linspace(50e-6,300e-6,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%：
dphi_l_to_ddelta=(2*((1/2)*(-2*L_c-2*L_s1+2*L_s3+2*delta)./...
    (sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    sqrt(L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2))-...
    (1/4)*((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3)^2-L_c^2).*(-2*L_c-2*L_s1+2*L_s3+2*delta)...
    ./(sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2).^(3/2))))./...
    sqrt(4-((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3)^2-L_c^2).^2./...
    ((L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3^2+2*L_s3.*delta+delta.^2)))...
    -L_s2./((L_c+L_s1-L_s3-delta).^2.*(1+L_s2.^2./(L_c+L_s1-L_s3-delta).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%
T_est=1./L_s2;
obj_Tdiff=(abs(dphi_l_to_ddelta)-T_est).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity_analysis_transmission_ratio_difference_L_s2
figure(5)
plot(L_s2*10^6,obj_Tdiff,'k-','LineWidth',2)
hold on
plot([297.1740,297.1740],[min(obj_Tdiff)-0.5*10^-21,max(obj_Tdiff)+0.25*10^-21],'r-.','LineWidth',3)
xlabel('翅前缘杆根部短杆(\itL_{s2}\rm) (um)','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
ylabel('\itf(L_{s2})=\rm(\itd\phi_{theor}/d\delta|_{\delta_{pp}\rm=250um} -\itT_{est}\rm)^2 ','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
% title('在压电驱动器峰值位移下的理论传动比和线性传动比的差的平方和','Color','k','FontSize',20,'FontName','Times','FontWeight','Bold')
grid on
box on
set(gca,'XTick',(45:25:305),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
axis(gca,[45,305,min(obj_Tdiff)-0.5*10^-21,max(obj_Tdiff)+0.25*10^-21])
set(gca,'Fontsize',19,'FontName','Times','FontWeight','Bold','Ycolor','k')

%% (4) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% lb = [50e-6,50e-6,50e-6,50e-6];          % Set lower bounds 
% ub = [1000e-6, 1000e-6,300e-6,1000e-6];  % Set upper bounds
% x0= [L_c,L_s1,L_s2,L_s3];  
optimal_four_link_size =[445.8109  564.1168  297.1740  439.1168]; % objTdiff_linear % delta_pp=250e-6;
L_c=optimal_four_link_size(1)*10^-6;   % m
L_s1=optimal_four_link_size(2)*10^-6;   % m 
L_s2=optimal_four_link_size(3)*10^-6;   % m
% L_s3=optimal_four_link_size(4)*10^-6;   % m
L_s3=linspace(50e-6,1000e-6,1000);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
dphi_l_to_ddelta=(2*((1/2)*(-2*L_c-2*L_s1+2*L_s3+2*delta)./...
    (sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3.^2).*...
    sqrt(L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3.^2+2*L_s3.*delta+delta.^2))-...
    (1/4)*((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3).^2-L_c^2).*(-2*L_c-2*L_s1+2*L_s3+2*delta)...
    ./(sqrt(L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3.^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3.^2+2*L_s3.*delta+delta.^2).^(3/2))))./...
    sqrt(4-((L_c+L_s1-L_s3-delta).^2+2*L_s2.^2+(L_s1-L_s3).^2-L_c^2).^2./...
    ((L_s2.^2+L_s1^2-2*L_s1*L_s3+L_s3.^2).*...
    (L_s2.^2+L_c^2+2*L_c*L_s1-2*L_c*L_s3-2*L_c.*delta+L_s1^2-2*L_s1*L_s3-2*L_s1.*delta+L_s3.^2+2*L_s3.*delta+delta.^2)))...
    -L_s2./((L_c+L_s1-L_s3-delta).^2.*(1+L_s2.^2./(L_c+L_s1-L_s3-delta).^2));
%%%%%%%%%%%%%%%%%%%%%%%%%
T_est=1./L_s2;
obj_Tdiff=(abs(dphi_l_to_ddelta)-T_est).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivity_analysis_transmission_ratio_difference_L_s3
figure(6)
plot(L_s3*10^6,obj_Tdiff,'k-','LineWidth',2)
hold on
plot([439.1168,439.1168],[min(obj_Tdiff)-0.25*10^6,max(obj_Tdiff)+0.25*10^6],'r-.','LineWidth',3)
xlabel('翅前缘杆根部中段杆(\itL_{s3}\rm) (um)','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
ylabel('\itf(L_{s2})=\rm(\itd\phi_{theor}/d\delta|_{\delta_{pp}\rm=250um} -\itT_{est}\rm)^2 ','Color','k','FontSize',24,'FontName','Times','FontWeight','Bold')
% title('在压电驱动器峰值位移下的理论传动比和线性传动比的差的平方和','Color','k','FontSize',20,'FontName','Times','FontWeight','Bold')
grid on
box on
set(gca,'XTick',(45:50:1005),'LineStyle','-','LineWidth',1.5)
% set(gca,'YTick',(0:0.4:2.8),'LineStyle','-','LineWidth',1.5)
set(gca,'LineStyle','-','LineWidth',1.5)
axis(gca,[45,1005,min(obj_Tdiff)-0.25*10^6,max(obj_Tdiff)+0.25*10^6])
set(gca,'Fontsize',19,'FontName','Times','FontWeight','Bold','Ycolor','k')



