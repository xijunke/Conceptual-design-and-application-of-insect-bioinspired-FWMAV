% Nonlinear_resonance_modeling
% standard harmonic oscillator motion of the actuator tip (x(t))
% x=0;    % 压电驱动器的尖端位移输出
% psi=0; % 扭转角
clear all;clc;
%%  翅膀惯性张量和虚质量惯性矩
I_xx=0.95*1.29;     % mg*mm^2
I_xz=0.29*2.8;       % mg*mm^2
I_zz=51.1;             % mg*mm^2
I_am=10.22;         % mg*mm^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 驱动器的参数
% t_pzt=135e-6;    % um % thickness of each PZT plate in the actuator
% t_cf=50e-6;        % um % the thickness of the actuator's central carbon fiber layer
% E_cf=340e9;     % Gpa % modulus of the carbon fiber
% L_pzt=9e-3;       % mm % the length of the active portion of the actuator
% L_r=0.25;           % the ratio of the length of the actuator's rigid extension to the length of the PZT
% w_nom=1.125e-3;  % mm % the mean width of the actuator  
% w_r=1.556;             % w_r is the width ratio (the width of the PZT's base divided by w_nom)
%%%%%%%%%%%%%
t_pzt=127e-6;
t_cf=40e-6; 
E_cf=350e9;
L_pzt=6.021e-3; 
L_r=5.979e-3/L_pzt;      % L_ext=5.979e-3;
w_nom=1.1755e-3;       % w_nom=(w0+w1)/2; % w0=1.569e-3; w1=0.782e-3;
w_r=1.569e-3/w_nom;  % w_r =1.3348;
%%%%%%%%%%%%%%%%%%%%%%%%%
L_act=L_pzt;       % L_act is the length of the PZT
L_ext=L_r*L_act;  % L_ext is the length of the alumina extension 
d_31=320e-12;          % pm/V  % d_31 is the piezo coefficient
E_pzt=62e9;  % pa  % E_pzt is the modulus of the PZT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 传动机构
T=3333;        % rad/m % transmission ratio
k_t=197.7;    % N/m % transmission stiffness
%%  the strain for actuator——应变
f=155; % Hz
omega=2*pi*f;
T_period=1/f;
%%%%%%%%%%%%%%%%%%%%%%%%%
Phi_amplitude=60*pi/180;
t=linspace(0,2*T_period,1000);
phi=Phi_amplitude*cos(omega*t); % flapping angle
x=phi/T;  % 压电驱动器的尖端位移输出
%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_amplitude=pi/4;
% Psi_amplitude=70*pi/180;
psi=Psi_amplitude*cos(omega*t-pi/2); % pitch angle
%%%%
figure(1)
plot(t,phi*180/pi,'r-',t,psi*180/pi,'b-','LineWidth',2)
xlabel('time (second)')
ylabel('Wingbeat angle (degree)')
legend('\phi(t)','\psi(t)')
title('悬飞FWMAV的翅拍运动角度')  % Wingbeat angle of hovering FWMAV
%%%%
figure(2)
plot(t,x*10^6,'r-','LineWidth',2)
xlabel('time (second)')
ylabel('The free deflection of actuator (um)')
legend('x(t)')
title('压电驱动器的尖端位移输出')
%%%%%%%%%%%%%%%%%%%%%%%%%
% The compressive strain in each of the two PZT plates——应变
epsilon_1=-(t_pzt+t_cf*x)/(L_pzt^2*(1+2*L_r));
epsilon_2=(t_pzt+t_cf*x)/(L_pzt^2*(1+2*L_r));
epsilon_0=-0.00047;  % Strain at which the transition occurs——2015-SMS
% epsilon_0=-0.047; % unit (%)  % the strain at which the modulus transitions from E_min to E_max
epsilon=-abs(epsilon_1); % the magnitude of the compressive strain in the actuator
E_min=38.5e9;  % Gpa % the minimum moduli of the PZT under varying strain
E_max=81e9;    % Gpa % the maximum moduli of the PZT under varying strain
a_1=8000;    % 无量纲 % the steepness of this transition
E_e=E_min-((E_max-E_min)/2)*log((1+exp(a_1*(epsilon_0-epsilon)))/(1+exp(a_1*epsilon_0)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The time-varying electric field applied to each of the PZT plates——场强
xi_0=0.4e-6;       % V/um % the field at which this transition occurs
V_pp=200;         % V % applied peak-to-peak voltage
xi_1_t=(0.5*V_pp/t_pzt)*(1-sin(2*pi*(1+15*t/2).*t));   % The compressive strain in each of the two PZT plates
xi_2_t=(0.5*V_pp/t_pzt)*(1+sin(2*pi*(1+15*t/2).*t));  % The compressive strain in each of the two PZT plates
figure(3)
plot(t,xi_1_t/10^6,'r-',t,xi_2_t/10^6,'b-','LineWidth',2)
xlabel('time (second)')
ylabel('External actuating electiric field (V/um)')
legend('\xi_{1}(t)','\xi_{2}(t)')
title('外加驱动电场的场强')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 形状因子:
GF=8*(1-w_r)^3*(1+2*(L_ext/L_act))/(-6*(w_r-1)*(-3+4*L_r*(w_r-1)+2*w_r)+3*(-2+2*L_r*(w_r-1)+w_r)^2*abs(log((2-w_r)/w_r)));
G=GF;       % a gain factor (defined in [9]) determined by the shape of the actuator
%%%%%%%%%%%%%%%%%%%%%%%%%   
% actuator stiffness——驱动器的刚度
% p_4 is a fitting coefficient to account for an imperfect attachment between the actuator and the airframe
P_4=0.79; % 拟合系数
% p_5 is set to 0.8 to bring the predicted free deflection of the actuator closer to its measured values 
% (in [9], the actuator deflections are slightly under-predicted)
P_5=0.8; % 拟合系数
K_act=P_4*P_5*((w_nom*G*(E_e*t_pzt*(1.5*t_cf^2+3*t_cf*t_pzt+2*t_pzt^2)+E_cf*t_cf^3/4))/(L_pzt^3*(1+2*L_r)));  
%%%%%%%%%%%%%%%%%%%%%%%%%
% The effective mass felt by the actuator——驱动器感受到的质量
M_e=(I_zz+I_am)*T^2;
% The stiffness felt by the actuator tip——驱动器尖端感受到的刚度
K=K_act+T^2*k_t;
%% The driving force provided by the actuator——驱动力的计算
% (derived from [9]) is given by the difference in forces provided by each
% PZT plate (F_act = F_1-F_2). F_1 and F_2 are computed according to (6).
%%%%%%%%%%%%%%%%%%%%%%%%%
% f_31min and f_31max are the minimum and maximum values of the piezoelectric stress per field coupling coefficient 
% in the field dependency of this piezoelectric coupling coefficient
f_31min=14;   % N/Vm
f_31max=29;  % N/Vm
a_2=-230;     % 无量纲  % a slight linear increase in the force with increasing compressive strain
% a_3=10e-5; % ——2015-SMS
a_3=10e-6;   % um/V the steepness of this transition 
a_4=69e-9;   % nm/V a slight linear decline in the piezo stress coefficient at high fields
%%%%%%%%%%%%%%%%%%%%%%%%%
% p_3 is a fitting coefficient to account  for an imperfect attachment between the actuator and the air frame
P_3=0.7;   %拟合系数
% P_3=0.75; % ——2015-SMS
F_1=P_3*((3*xi_1_t.*G*w_nom*t_pzt*(t_pzt+t_cf).*(1+a_2*epsilon_1).*...
          (f_31min+(f_31max*(1-a_4*xi_1_t)-f_31min).*exp(a_3*(xi_1_t-xi_0))))./...
          (4*L_pzt*(1+exp(a_3*(xi_1_t-xi_0)))));
F_2=P_3*((3*xi_2_t.*G*w_nom*t_pzt*(t_pzt+t_cf).*(1+a_2*epsilon_2).*...
          (f_31min+(f_31max*(1-a_4*xi_2_t)-f_31min).*exp(a_3*(xi_2_t-xi_0))))./...
          (4*L_pzt*(1+exp(a_3*(xi_2_t-xi_0)))));
F_act=F_1-F_2;
figure(4)
plot(t,F_1*10^3,'r-',t,F_2*10^3,'b-',t,F_act*10^3,'k-','LineWidth',2)
xlabel('time (second)')
ylabel('Actuating force (mN)')
legend('F_1(t)','F_2(t)','F_{act}(t)')
title('压电驱动器产生的驱动力') % The driving force provided by the actuator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀的气动参数
rou=1.225e-3;                  % density of air %单位是Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3
A_w=54.59e-6;                 % mm^2  % wing area
R_cp=(1.42*9.56)*10^-3; % mm % radial center of pressure
r_nd2=0.564;                    % 无量纲 % standard wing-shape parameter
R=17e-3;                          % mm % distance from the transmission flexure to the wing tip
W_av=A_w/R;                  % mean width of the wing
w_nomd2=0.43;     % 弦向面积矩
W=4.5e-3; % mm   % 类似于弦长
W_cp1=(0.39*1.93)*10^-3;  % mm % 弦向压心
p_min=0.132;         % 无量纲 % fitting parameters
p_max=0.151;        % 无量纲 % fitting parameters
k_wh=(0.985*1.52)*10^-6;  % uNm/rad=*10^-6Nm/rad  % the stiffness of the passive wing hinge
b_min=0.39*0.4;               % 无量纲 % the drag at  psi= 0 degree.
b_max=1.25*3.4;              % 无量纲 % the drag at  psi=90 degree
b_1=(b_max+b_min)/2;
b_2=(b_max-b_min)/2;
b_L=1.8;                  % 无量纲
%% 压心和气动阻尼系数等
C_A=0.5*rou*A_w*R_cp*r_nd2^2*T^3*R^2*b_1;  % standard drag coefficients
C_B=0.5*rou*A_w*R_cp*r_nd2^2*T^3*R^2*b_2;  % standard drag coefficients
C_C=0.5*rou*A_w*W_cp1*w_nomd2^2*W^2*b_max;  % damping coefficient for the wing hinge rotation
W_cp2=W_av*(p_min+(p_max-p_min)*abs(1-psi*(2/pi))); %  the effective chord-wise center of pressure as a function of wing hinge angle
C_D=0.5*rou*A_w*W_cp2*r_nd2^2*T^2*R^2.*((b_1+b_2*cos(2*psi)).*cos(psi)+b_L*sin(2*psi).*sin(psi));
figure(5)
plot(t,C_D*10^0,'r-','LineWidth',2)
xlabel('time (second)')
ylabel('C_D (-)')
legend('C_D(t)')
title('平动气动力矩系数')
%% 功率消耗和效率
C_act=10.1e-9;   % nF
P=0.5*C_act*V_pp^2*f;
M_fwmav=80; % mg
% M_fwmav=87; % mg
eta=P/M_fwmav;  % 












