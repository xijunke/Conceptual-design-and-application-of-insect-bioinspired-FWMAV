%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lift and drag coefficients with Alpha——随攻角变化的升阻力系数
% lift_drag_coefficients_with_alpha
% clear all; clc;
% alpha=atan2(-w_y,w_z);    %Aerodynamic angle of attack  %这个要求事先求出绕各轴的翅运动角速率
beta=pi/4; f=0.25; %\Hz
delta=pi/4;    %相差delta:  提前模式:delta=pi/4; 延迟模式:delta=-pi/4; 对称模式:delta=0;
t=linspace(0,10,200);
alpha=(beta*sin(2*pi*f*t+delta))*180/pi;
% figure(1)
% plot(t,alpha,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 先使用函数fzero,求零点值，再使用函数find找到其位置索引
options=optimset('Display','off');
t1_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[3,4],options); 
t2_0=fzero(inline('45*sin(pi/4 + (pi*t)/2)','t'),[5,6],options);
% Result:   t1_0=3.5000;   t2_0=5.5000;
indx1=find(t>3.45 & t<=3.5);  %条件语句
indx2=find(t>5.45 & t<=5.5);
% Result:   indx1=70;   indx2=110;
alpha_2=alpha(1,indx1:indx2);      % 半个周期    %  alpha_2∈(0°,45°)
% figure(2)
% plot(t(1,indx1:indx2),alpha_2,'ro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 求出alpha_2的最大值和最小值及其位置索引
[alpha_2_max, locat_max]=max(alpha_2);   % 输出:  alpha_2_max =44.9716;  locat_max =22;
% [alpha_2_min locat_min]=min(alpha_2);  % 输出:  alpha_2_min =-2.3078;    locat_min =1;
m=length(alpha_2);                                                     %m =41
alpha_3=alpha_2(1,locat_max:m)+45;     % size(alpha_3)  % (1*20)                  %  alpha_3∈(45°,90°)
[alpha_new,index]=sort(alpha_3);           %按升序排列，返回两个数据
alpha_4=[0, alpha_2(1,2:locat_max-1), alpha_new];  % size(alpha_4)% (1*41)   %  alpha_4∈(0°,90°)
% figure(3)
% plot(t(1,indx1:indx2),alpha_4,'m-d')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 各种升阻力系数的对比分析
%% (1) 1999-Science-MH Dickinson——三维实验测试升阻力系数
% C_rot_theo=pi*(0.75-x_0nd(r));　  % 转动环量气动力系数
C_L_dickinson=0.225+1.58*sin((2.13*alpha_4-7.2)*pi/180);  
C_D_dickinson=1.92-1.55*cos((2.04*alpha_4-9.82)*pi/180); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (8) 波尔豪森采用前缘吸力比拟法建立的升力系数计算模型
% C_L_polhamus=K_p*sin(alpha)*(cos(alpha)).^2+K_v*(cos(alpha))*sin(alpha).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (9) 2014-JRSI-Nabawy
% C_La2d=2*pi;    % 单位是rad^(-1), 这里2*pi (rad^(-1))= 0.11 (deg^(-1));   
C_La2d=5.15662;  % 单位是rad^(-1), 针对典型昆虫雷诺数的平板, 建议使用 0.09 (deg^(-1)) % 采用前缘吸力比拟法建立的波尔豪森升力系数模型C_L_polhamus1完全重合
lambda=0.70; % 采用这个缩放因子与采用Dickinson拟合的系数计算而得的升重比和气动功率密度较接近,即P_asterisk=29.7503;  L =1.2761;
x=[3.004,0.8854,0.3289,0.356737];
% x =[3.6837,1.7620,2.0000,0.0055];
R_wing=x(1);
C_aver=x(2);
xr0=x(3);
C_maxyaxis=x(4);
C_periEdge=C_perimeter_to_Edge_correction(R_wing,C_aver,xr0,C_maxyaxis);
E=lambda*C_periEdge;   % E=C/2/R; 针对初始果蝇翅膀 E_C_peri=1.146071666; 
k=1.51;  % 针对果蝇
% for the span of a single  wing 
AR=3.40158;     % 注意这里的展弦比AR∈(3,5)是比较好的。% 针对果蝇
% AR =3.7643;
C_La_Ny=C_La2d/(E+k*C_La2d/(pi*AR));     % 由普朗特升力线理论获得的三维翅膀升力曲线斜率  % 输出: C_La_Ny =2.7506;
% C_L_nabawy=(0.5*C_La2d/(E+k*C_La2d/(pi*AR)))*sin(2*alpha_4);
C_L_nabawy=0.5*C_La_Ny*sin(2*alpha_4*pi/180);   
% C_D_nabawy=C_Lnabawy.*tan(alpha_4);   % C_Dnabawy=2*C_T*(sin(alpha))^2=C_La_Ny*(sin(alpha))^2;  
C_D_nabawy=C_La_Ny*(sin(alpha_4*pi/180)).^2;  
% 下面是采用前缘吸力比拟法建立的波尔豪森升力系数模型
% C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha)).*(cos(alpha)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha));
C_L_polhamus1=(0.5*C_La_Ny*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-k*C_La_Ny/(pi*AR)).*sin(alpha_4*pi/180)); 
% figure(11)
% plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_nabawy,'r-',alpha_4,C_D_nabawy,'b-',alpha_4,C_L_polhamus1,'g-')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (10) 2014-AST-Taha
a_0=2*pi;
%  AR being based on one wing: AR=R^2/S;
AR=2.8;  % 注意这里的展弦比AR＜3是比较好的。 
C_La_Ta=pi*AR/(1+sqrt((pi*AR/a_0)^2+1));     % 由扩展升力线获得的三维翅膀升力曲线斜率   % 输出:C_La_Ta =3.2334;
% C_L_taha=(0.5*pi*AR/(1+sqrt((pi*AR/a_0)^2+1))).*sin(2*alpha_4*pi/180);
C_L_taha=0.5*C_La_Ta*sin(2*alpha_4*pi/180);
C_D_taha=C_L_taha.*tan(alpha_4*pi/180);
% 下面是采用前缘吸力比拟法建立的波尔豪森升力系数模型
C_L_polhamus2=(0.5*C_La_Ta*sin(2*alpha_4*pi/180)).*(cos(alpha_4*pi/180)+(1-C_La_Ta/(pi*AR)).*sin(alpha_4*pi/180));
figure(12)
C_LD=plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_D_dickinson,'b*',alpha_4,C_L_taha,'r-',alpha_4,C_D_taha,'b',...
       alpha_4,C_L_nabawy,'k-',alpha_4,C_D_nabawy,'k-',alpha_4,C_L_polhamus1,'c-',alpha_4,C_L_polhamus2,'m-');
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha ) and \it C_D (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{D,dickinson}','\itC_{L,taha}','\itC_{D,taha}','\itC_{L,nabawy}',...
            '\itC_{D,nabawy}','\itC_{L,polhamus1}','\itC_{L,polhamus2}')
set(C_LD,'LineWidth',2)
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一部分——升力系数的对比分析
figure(13)    %  subplot(224)
C_L=plot(alpha_4,C_L_dickinson,'r*',alpha_4,C_L_nabawy,'k-',alpha_4,C_L_taha,'r:',...
                alpha_4,C_L_polhamus1,'c-',alpha_4,C_L_polhamus2,'m-');
xlabel('\it\alpha (deg.)')
ylabel('\it C_L (\alpha )')
title('Aerodynamic coefficients of lift \itvs. \alpha \rm for flapping wing')
legend('\itC_{L,dickinson}','\itC_{L,nabawy}','\itC_{L,polhamus1}','\itC_{L,taha}','\itC_{L,polhamus2}')
set(C_L,'LineWidth',2) 
grid on
%axis([xmin,xmax,ymin,ymax])
%% 第二部分——阻力系数的对比分析
figure(14)
C_D=plot(alpha_4,C_D_dickinson,'b*',alpha_4,C_D_nabawy,'k-',alpha_4,C_D_taha,'b:');
xlabel('\it\alpha (deg.)')
ylabel('\it C_D (\alpha )')
title('Aerodynamic coefficients of drag \itvs. \alpha \rm for flapping wing')
legend('\itC_{D,dickinson}','\itC_{D,nabawy}','\itC_{D,taha}')
set(C_D,'LineWidth',2)   
grid on
%% 坐标轴显示区间的设定
v_axis=axis;     %axis([xmin,xmax,ymin,ymax])
%Result: v_axis =  0   100     0     4
% v_axis(1)=0;          %指定x轴的最小值
% v_axis(2)=100;      %指定x轴的最大值
v_axis(3)=-0.5;          %指定y轴的最小值
v_axis(4)=3.5;           %指定y轴的最大值
axis(v_axis);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%