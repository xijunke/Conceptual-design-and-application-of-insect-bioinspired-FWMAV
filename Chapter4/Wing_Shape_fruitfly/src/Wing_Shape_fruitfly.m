% function wing_para_output=wing_shape_fruitfly_sixteen_good_2()
% 创建时间——2014年6月19日,0:16:04
% 修改时间——2014年6月14日,23:35:03
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];   % 变量
% x =[3.004,0.8854,0.3289,0]; 
% x =[3.004,0.8854,0.3289,0.356737];  % R_wing=3.004;  C_aver=0.8854;   % 原始果蝇翅膀
% x =[3.004,0.8854,0.3289,0.25];
x =[3.1267,1.1578,0.7896,0.1624]; 
% x=[2.5088,1.1424,0.7124,0.0458];
% x =[2,0.5,0.3289,0.0]; % R_wing=2;  C_aver=0.5;               % 最小翅膀——下限
%  x =[4,2,0.3289,0.0]; % R_wing=4;  C_aver=2;                 % 最大翅膀——上限
% x =[100,33.7,0.3289,0.0]; % R_wing=100;  C_aver=33.7;   % mm ——@wing 1
%%%%%%%%%%%%%%%%%%%%%%%%%
R_wing=x(1);
C_aver=x(2);
xr0=x(3);               % x-root offset  \mm
C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅形貌参数化——由Hedrick程序计算而得
%拟合得到 R_wingeff=3.3328-0.3289=3.004;  Hedrick程序计算而得: R_wingeff=3.007; Science数据是: R_wingeff=2.99; 
R_wingeff=3.004;    %有效翅膀长度(mm)  
% Hedrick程序计算而得: C_avereff=0.884;     % Science数据是C_avereff=0.9468mm;  
C_avereff=0.8854;  % 单位:mm---见后文—C_avereff=C_aver0 =0.8854;—由前后缘实际拟合曲线函数相差求积分均值而得
% 可以考虑已知展弦比，求解气动力最优的翅形貌
% AR=R_wingeff/C_avereff;     %aspect ratio: AR=R^2/A_w;  这里输出为: AR=3.40158;  % Science数据是: AR=3.1519;
% C_avereff=R_wingeff/AR;     %mean chord length: C_aver=A_w/R=R^2/AR/R=R/AR; % C_aver=0.884mm;
% A_w=R_wingeff^2/AR;        %Area of wing: 这里输出为: A_w=2.66mm^2;   %RJ Wood设计的翅膀: A_w=2.84345 mm^2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 第一种方式求解无量纲气动力分量——翅形貌因子——需要输入的参数数据
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  无量纲展向偏置距离
% yr=0;                          % y-root offset  \mm
% yr_nd=yr/C_avereff;   % y-root offset  无量纲弦向偏置距离;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 已知r2_nd,求解r1_nd, r3_nd, 和r3_nd3
r2_nd=0.5801;      %已知无量纲二阶面积矩r2_nd;    %由Hedrick程序计算而得non_dimensional_2ndMoment: 0.5801(更准) 
r1_nd=1.106*r2_nd^1.366;  %2013-ICRA-Deng XY   % 输出: r1_nd =0.5257
% r3_nd=0.9*r1_nd^0.581; %1984-JEB-Ellingdon;  % 输出: r3_nd =0.6194 % 由Hedrick程序计算而得 non_dimensional_3rdMoment: 0.6184(更准)     
% r3_nd3=r3_nd^3;                   % r3_nd3=r3_nd^3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff1=3.007;            % Hedrick程序计算而得
% R_wingeff1=R_wingeff;            % 
xr_nd1=xr/R_wingeff1;      % x-root offset  无量纲展向偏置距离
F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %无量纲气动力分量F_nd1, 这里输出: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 第二种方式求解无量纲气动力分量——wing_shape——翅前后缘拟合函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅根-偏离-坐标系原点的距离
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_maxylb =0.464385778290230;
% C_maxy25 =0.138924474377504;  % 针对程序wing_model_88_yaxis有: 第122行; C_maxy =0.1389; 
% C_maxyub =-0.186536829535222; 
% C_maxy=0;
C_maxy=C_maxyaxis;
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% 前后缘6阶多项式拟合函数
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
C_rx=yr_lead-yr_trail;              % 正确——量纲化实际弦长分布
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 采用前缘拟合函数求解——实际平均弦长=面积/R_wingeff
wing_aera=trapz(x,C_rx);             %输出: wing_aera =2.6597; % mm^2
C_aver0=wing_aera/R_wingeff;   % 输出量纲化平均弦长: C_avereff=C_aver0 =0.8854; % mm
%%%%%%%%%%%%%%%%%%%%%%%%%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% 采用多项式函数拟合获得——无量纲前后缘分布函数和弦长分布函数
% (a) 无量纲前缘分布函数
r_nd=(x-xr)/R_wingeff;
yr_leadnd0=yr_lead/C_avereff;
P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  无量纲后缘分布函数
yr_trailnd0=yr_trail/C_avereff;
P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  无量纲弦长分布函数
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % 多项式系数  % Cr_nd2=polyval(Coeff,r_nd1);
% cftool 
syms r_nd   % 无量纲弦长分布为6阶多项式——转换必须有这条指令
yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
Cr_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 注意——该段程序切记不得修改，前提只要保证输入正确的无量纲弦长分布即可。
%以下的公式应使用合理的无量纲的弦长分布公式C_nd
C_nd=Cr_nd;
R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); %二阶面积矩的回转半径的平方
R1nd1=double(int(r_nd*C_nd,r_nd,0,1));     %一阶面积矩的回转半径
S_nd=double(int(C_nd,r_nd,0,1));                %无量纲翅面积 % S_nd =1.0000;
disp(['二阶面积矩的回转半径的平方: r2_2nd=' num2str(R2nd2)  ' 量纲单位是mm^4'])
disp(['二阶面积矩的回转半径: r_2nd=' num2str(sqrt(R2nd2))  ' 量纲单位是mm^3'])
disp(['一阶面积矩的回转半径: r_1nd=' num2str(R1nd1)  ' 量纲单位是mm^3'])
disp(['无量纲翅面积Swing_nd=' num2str(S_nd)  ' 量纲单位是mm^2'])
% 假设xr_nd=0，则自变量r_nd的取值范围为:r_nd∈[0,1], 在取1时r_nd=R=11.95 /mm
fx1=(r_nd+xr_nd)^2*C_nd;    % 无量纲气动力F_nd的原始被积函数
% fx2=vpa(fx1,5)
fx3=expand(fx1);
F_ndTrans=double(int(fx3,r_nd,0,1));                    % Result: F_ndTrans =0.46391;
disp(['无量纲气动力F_ndTrans=' num2str(F_ndTrans)  ' 量纲单位是mm^4'])
F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;    %使用这句计算结果也正确; 输出:F_nd2 =0.46391;
disp(['无量纲气动力F_ndTrans=' num2str(F_nd2)  ' 量纲单位是mm^4'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输出:
% 二阶面积矩的回转半径的平方: r2_2nd=0.3358 量纲单位是mm^4
% 二阶面积矩的回转半径: r_2nd=0.57949 量纲单位是mm^3  ——————————需要该数据
% 一阶面积矩的回转半径: r_1nd=0.53032 量纲单位是mm^3
% 无量纲翅面积Swing_nd=1 量纲单位是mm^2
% 无量纲气动力F_nd=0.46391 量纲单位是mm^4    ——————————需要该数据
% 无量纲气动力F_nd=0.46391 量纲单位是mm^4    ——————————需要该数据
%  对比下面由 (1) 第一种方式求解获得的无量纲气动力分量
% F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %无量纲气动力分量F_nd1, 这里输出: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 果蝇比例翅膀——fruitgfly_scaled_wing
%%%%%%%%%%%%%%%%%%%%%%%%%
% PHI=74.5*pi/180;
% f3 =0.3378;                        % Hz
% xr=0.3289;                           % x-root offset  \mm                    
xr_nd1=xr0/R_wing;               % x-root offset  无量纲展向偏置距离
% xr0=14.5;                            % mm——参考Nano 17数据
% xr_nd1=xr0/R_wing;           % x-root offset  无量纲展向偏置距离; 这里xr_nd被更新 xr_nd=0.03222;@wing 2;
% 注意前面针对 F_ndTrans的计算进行了更新
F_ndTrans=R2nd2+2*xr_nd1*R1nd1+xr_nd1^2;    % 更新: F_ndTrans=0.5106;@wing 1; or F_ndTrans=0.3710;@wing 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 绘制拟合好的前后缘
% syms r_nd     % r_nd=(x-xr)/R_wingeff;  % r_nd∈(0,1)
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd+0.000835197;%@最大前缘点距翅根翅尖;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd+0.000893019;%@最大前缘点距翅根翅尖;
f_x_lead=C_aver*yr_leadnd;
f_x_trail =C_aver*yr_trailnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_x_lead1=inline(vectorize(f_x_lead),'r_nd');  %数值函数
f_x_trail1=inline(vectorize(f_x_trail),'r_nd');   %数值函数
r_nd1=linspace(0,1,200);     % x=linspace(R_proximal,R_distal,200);  % x=r_nd*R_wingeff+xr;
f_x_lead2=f_x_lead1(r_nd1); 
f_x_trail2=f_x_trail1(r_nd1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (一) 分辨平均弦长和翅膀扭转轴所在位置 
% C_avereff=0.8854;  % mm
% % (1) 最大前缘点对应的片条坐标和弦长         % 原始果蝇翅膀@C_maxy=0;
C_lead_ymax=max(f_x_lead2);                        % 输出: C_lead_ymax=0.3255; k_leadmax=644;% @最大前缘点到翅根翅尖连线的距离;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (2) 最小后缘点对应的片条坐标和弦长      % 原始果蝇翅膀@C_maxy=0;
C_trail_ymin=min(f_x_trail2);                         % 输出: C_trail_ymin = -0.1361;  k_trailmin =409;% @最小后缘点到翅根翅尖连线的距离;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_max_LtoT=C_lead_ymax-C_trail_ymin;      % C_max_LtoT =1.3018;% 原始果蝇翅膀@C_maxy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rx=r_nd1*R_wing+xr0;  % rx=linspace(0,1,200)*R_wingeff+xr; 
r_x_max=max(rx);          % r_x_max=100.3289;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cftool % 这种方式不提倡
% % (1) 拟合多项式函数——由f_x_lead2(r_x)拟合而得yr_leadnd2(r_x)
r_x=linspace(xr0,r_x_max,200);
yr_leadnd2=-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166;
% Linear model Poly6: f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 +p6*x + p7
% Coefficients (with 95% confidence bounds):
% p1=-2.307e-009(-2.307e-009, -2.307e-009);p2=9.027e-007(9.027e-007, 9.027e-007); p3=-0.0001408(-0.0001408, -0.0001408)
% p4=0.01104(0.01104, 0.01104); p5=-0.4529(-0.4529, -0.4529); p6=9.311(9.311, 9.311); p7=-67.75  (-67.75, -67.75)
% Goodness of fit: SSE: 1.817e-022; R-square: 1;  Adjusted R-square: 1; RMSE: 9.704e-013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
% % (2) 拟合多项式函数——由f_x_trail2(r_x)拟合而得yr_trailnd2(r_x)
yr_trailnd2= -9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034;
% Linear model Poly6: f(x)=p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x + p7
% Coefficients (with 95% confidence bounds):
% p1 =-9.314e-010  (-9.314e-010, -9.314e-010);p2 =4.891e-007  (4.891e-007, 4.891e-007);p3 =-9.514e-005  (-9.514e-005, -9.514e-005);
% p4 =0.008779  (0.008779, 0.008779); p5 =-0.3877  (-0.3877, -0.3877); p6 = 6.714  (6.714, 6.714); p7 = -38.66  (-38.66, -38.66)
% Goodness of fit: SSE: 2.433e-022; R-square: 1;  Adjusted R-square: 1;  RMSE: 1.123e-012;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)   % 绘制拟合好的前后缘 & 翅膀放大之后的压心位置
hold on;
%绘制拟合好的前后缘—前后缘的纵轴坐标当横坐标，横轴坐标当纵坐标——类似于扭转90°
plot(rx,f_x_lead2,'k-',rx,f_x_trail2,'k-','LineWidth',3);  hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制坐标系
%% 翅肩部——Shoulder-X轴——Z轴
quiver(-0.25,0,...                                                                % 起点
            4.25,0,1,'k-','LineWidth',2);   hold on;                    % 终点   1倍显示——X轴     
quiver(0,-1.5,...                                                                 % 起点
           0,2.5,1,'k-','LineWidth',2);   hold on;                         % 终点   1倍显示 ——Z轴
% 翅根部——root-x轴——x轴    
% quiver(rx(1,1), f_x_lead2(1,1),...                                                   % 起点
%            rx(1,length(rx)),0,1,'r-.','LineWidth',2);   hold on;    %终点 ——x'轴
quiver(rx(1,1), f_x_lead2(1,1),...                                                   % 起点
           rx(1,length(rx))-0.5,0,1,'r-.','LineWidth',2);   hold on;    %终点 ——x'轴
quiver(rx(1,1), f_x_lead2(1,1)-1,...                                                % 起点
           0,2*abs(f_x_lead2(1,1))+1.75,1,'r-.','LineWidth',2);   hold on;% 终点 ——z'轴
% pitch-轴的绘制       
% quiver(0, C_maxyaxis*C_max_LtoT,...                                  % 起点
%            rx(1,length(rx))+xr0,0,1,'r-.','LineWidth',2);   hold on;    %终点 ——x'轴
quiver(0, C_maxyaxis*C_max_LtoT,...                                  % 起点
           rx(1,length(rx))+xr0-0.5,0,1,'r-.','LineWidth',2);   hold on;    %终点 ——x'轴    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
% plot(r_x,yr_leadnd2,'-.r',r_x,yr_trailnd2,'-.r','LineWidth',3); hold on
xlabel('右侧翅膀的展向(y坐标)')
ylabel('右侧翅膀的弦向(x坐标)')
title('采用前后缘拟合函数绘制右侧翅膀形貌, 扭转轴的位置和压心随时间和空间的分布')
grid off
axis([-0.3,4.75,-1.25,1.05])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 绘制翅膀质心坐标
XC=1.920243385-xr; 
% C_maxy=0.138924474377504;
% YC=-0.149785466+C_maxy;    % YC=-0.288709941;
YC=-0.149785466;    % YC=-0.288709941;
% R_ratio=(R_wing+xr1)/(R_wingeff+xr1);    % R_wingeff=3.004;  % R_wing=100;   % mm  %  XXXX
R_ratio=R_wing/R_wingeff;    % R_wingeff=3.004;    % R_wing=100;   % mm   % R_ratio =33.2889;
C_ratio=C_aver/C_avereff;     % C_avereff=0.8854;    % C_aver=33.7;   % mm   % C_ratio =38.0619;
% delta=xr0-xr;
% XC_tran=R_ratio*XC+xr+delta;      % XC_tran =53.3030
XC_tran=R_ratio*XC+xr0;      % XC_tran =53.3030
YC_tran=C_ratio*YC;              % YC_tran =-5.7011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_leadnd2=inline('-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166','r_x');
yr_trailnd2= inline('-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034','r_x');
Area=quadl(yr_leadnd2,0.3289,100.3289)-quadl(yr_trailnd2,0.3289,100.3289);
integrand_xc=inline('r_x.*((-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166)-(-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034))','r_x'); %数值函数
integrand_zc=inline('(-2.307e-009*r_x.^6+7.065e-007*r_x.^5-8.38e-005*r_x.^4+0.004742*r_x.^3-0.1288*r_x.^2+1.698 *r_x-0.5166).^2-(-9.314e-010*r_x.^6+4.099e-007*r_x.^5-6.329e-005*r_x.^4+0.004316*r_x.^3 -0.1158*r_x.^2+0.02573*r_x+0.034).^2','r_x'); %数值函数
X_C=quadl(integrand_xc,xr0,r_x_max)/Area+xr0;                 % X_C =53.8010;
Z_C=quadl(integrand_zc,xr0,r_x_max)/(2*Area);                   % Z_C =-5.4617;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(X_C,Z_C,'rd',XC_tran,YC_tran,'dk','LineWidth',4)
plot(XC_tran,YC_tran,'dk','LineWidth',4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wingeff=R_wing;
C_avereff=C_aver;
%% 第一部分——平动环量气动力和力矩参数
% (1) 平动环量气动力参数——绕翅平面下的法向
% Rou=0.88*10^(3)*10^(-3);     % mineral oil——矿物油的密度
% % % Rou=1.225;                         % 单位是Kg/m^3  ——空气的密度
Rou=1.225*10^(-3);           % 单位是Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3——空气的密度
% Coeff_liftdragF_N=6.8201e-012——单位是:mg*mm^-3*mm^4= 10^(-9) kg*m
Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_ndTrans;  % mg*mm
% (2) 平动环量气动力矩参数——绕翅平面下的展向轴
% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]=mg.mm^2: [10^12*10^(-12)]= mg.mm^2
M_xaercoeff=(1/2)*Rou*C_avereff^2*R_wingeff^3*F_ndTrans;   % 旋转轴气动力矩系数  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) 平动环量气动力矩参数——I1y——绕翅平面坐标系下的弦向轴
fx4=(r_nd+xr_nd)^3*C_nd;                                 % 这里C_nd由前后缘多项式函数之差给出
fr_nd5=expand(fx4);
I1=double(int(fr_nd5,0,1));                                 % 无量纲,量纲化单位是mm^5; 
I1y=(1/2)*Rou*C_avereff*R_wingeff^4*I1;         % 单位是 mg.mm^2  % I1y=0.0162
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二部分—— 针对气动阻尼力矩——求解有效力臂的无量纲位置Y_nd
% % 翅根-偏离-坐标系原点的距离
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
% x=linspace(R_proximal,R_distal,200);
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 方案(1)——由前后缘函数的无量纲化yr_leadnd——yr_trailnd——求得
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd+0.000835197;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd+0.000893019;
%%%%%%%%%%%%%%%%%%%
C_maxy=0.138924474377504; % z_r, 即0.25C_max_LtoT-Pitch轴到翅根的z-向距离 % C_max_LtoT =1.3018;% 原始果蝇翅膀@最大前缘点到最小后缘点的距离@C_maxy=0;
C_maxynd=C_maxy/C_avereff;
yr_leadnd=yr_leadnd-C_maxynd;
yr_trailnd=yr_trailnd-C_maxynd;
%%%%%%%%%%%%%%%%%%%
% yr_nd1=expand((yr_leadnd^4+yr_trailnd^4)/4);     %被积函数2
% Z_rnd=double(int(yr_nd1,r_nd,0,1));                         % 输出: Y_nd=0.08802(old);  %  mm  % Z_rnd =0.1626;
% 方案(2)——由前缘函数的无量纲化yr_leadnd和无量纲化弦长分布C_nd求得
% y0=yr_trailnd;         % 此时输出Z_nd=0.08802; ——见方案(1)
y0=yr_leadnd-C_nd;   % 此时输出Z_nd=0.08802;
y1=yr_leadnd;
yr_nd2=(y1^4+y0^4)/4;   % 注意这里的扭转轴直接通过翅尖和翅根的连线，前缘函数的翅膀展向区间唯一
Z_rnd=double(int(yr_nd2,r_nd,0,1));                          % 输出: Z_nd=0.08802(old);   %  mm  % Z_rnd =0.1627;
disp(['有效力臂的无量纲位置Z_nd=' num2str(Z_rnd)  ' 量纲单位是mm'])  % 结果可由Maple求得
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 转动气动阻尼力矩参数——绕翅平面下的展向轴
% 下面为转动气动阻尼力矩系数% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]= mg.mm^2
M_xrdcoeff=(1/2)*Rou*C_avereff^4*R_wingeff*Z_rnd;  %M_xrdcoeff=9.9528e-005=0.0001;
% (2) 转动气动阻尼力矩参数——绕翅平面下的弦向轴
% 下面为转动气动阻尼力矩系数% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]= mg.mm^2
fx16=(r_nd+xr_nd)*C_nd^3;                             % 这里C_nd由前后缘多项式函数之差给出
fr_nd17=expand(fx16);
I8=double(int(fr_nd17,0,1));                                  % 无量纲,量纲化单位是mm^5;
I8z=I8;           % 单位是 mg.mm^2
X_rnd=I8z;
M_zrdcoeff=(1/6)*Rou*C_avereff^3*R_wingeff^2*X_rnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第三部分——转动环量气动力和力矩参数
fx6=(r_nd+xr_nd)*C_nd^2;                                  % 这里C_nd由前后缘多项式函数之差给出   % 无量纲气动力F_ndRot的原始被积函数
fr_nd7=expand(fx6);
% F_ndRot1=double(int(fr_nd7,r_nd,0,1))           % 这个积分也行哦
F_ndRot=double(int(fr_nd7,0,1));                        % 无量纲,量纲化单位是mm^4; F_ndRot=0.7895
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 转动环量气动力参数——F_yrotcoeff——绕翅平面下的法向   % 单位是mg*mm^-3*mm^4= kg.m 10^(-9)
F_yrotcoeff=(1/2)*Rou*C_avereff^2*R_wingeff^2*F_ndRot;  % 单位是mg*mm
% (2) 转动环量气动力矩参数——M_xRotcoeff——绕翅平面下的展向轴
% mg*mm^-3*mm^5=Kg.m^2:[10^(-12)]=[10^12*10^(-12)] mg.mm^2
M_xRotcoeff=(1/2)*Rou*C_avereff^3*R_wingeff^2*F_ndRot;   % 转动环量气动力矩系数——绕翅平面下的扭转轴
% (3) 转动环量气动力矩参数——I6y——绕翅平面下的弦向轴
fx8=(r_nd+xr_nd)^2*C_nd^2;                            % 这里C_nd由前后缘多项式函数之差给出
fr_nd9=expand(fx8);
I2=double(int(fr_nd9,0,1));                                 % 无量纲,量纲化单位是mm^5;
I2y=(1/2)*Rou*C_avereff^2*R_wingeff^3*I2;     % 单位是 mg.mm^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第四部分——虚质量效应的系数求解和虚质量转动惯量——虚质量气动力和力矩参数
% C_maxy=0.138924474377504;
% C_maxynd=C_maxy/C_avereff;
yr_hnd=C_nd/2-yr_leadnd+C_maxynd; %转动轴偏离弦向中点的偏移量坐标,在扭转轴向上偏移C_maxy之后, 虚质量作用力力臂将全变为正值
% 非零虚质量系数
% a=C_nd/2;
% lambda_z=pi*Rou*a^2;
% lambda_zw=-pi*Rou*a^2*yr_hnd;
% lambda_w=pi*Rou*a^2*yr_hnd^2+pi*Rou*a^4/8;
% dOmega0=-r*(domega_y-omega_x*omega_z);           %法向加速度
% Z0=-lambda_z*dOmega0-lambda_zw*domega_x;      %虚质量气动力
% Y0=0;                                                                            %虚质量气动力
% M0=-lambda_zw*dOmega0-lambda_w*domega_x;    %虚质量气动力矩
%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 虚质量气动力矩参数——虚质量转动惯量I_xzam——绕翅平面下的扭转轴
fx10=(r_nd+xr_nd)*C_nd^2*yr_hnd;
fr_nd11=expand(fx10);
I_xzamnd=int(fr_nd11,0,1);
I3=double(I_xzamnd);                                         % 无量纲,量纲化单位是mm^5;          % I1 =0.2666;  
I_xzam=pi*Rou*C_avereff^3*R_wingeff^2*I3/4; % 单位是 mg.mm^2;   % I_xzam =0.0007(old);  I_xzam =0.0021;
% (2) 虚质量气动力矩参数——虚质量转动惯量I_xxam——绕翅平面下的扭转轴
fx12=C_nd^2*(yr_hnd^2+C_nd^2/32);
fr_nd13=expand(fx12);
I_xxamnd=int(fr_nd13,0,1);
I4=double(I_xxamnd);                                         % 无量纲,量纲化单位是mm^5;          % I2 =0.1637;  
I_xxam=pi*Rou*C_avereff^4*R_wingeff*I4/4;      % 单位是 mg.mm^2;   % I_xxam=0.0002(old);  I_xxam =6.0288e-004;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) 虚质量气动力参数——I5z
I5z=pi*Rou*C_avereff^2*R_wingeff^2*F_ndRot/4;    % 单位是 mg.mm
I_xzam2=I_xzam+I5z*C_maxynd^2;                   % I_xzam2 =0.0023;
% (4) 虚质量气动力参数——I6z
fx14=C_nd^2*yr_hnd;
fr_nd15=expand(fx14);
I6=double(int(fr_nd15,0,1));                               % 无量纲,量纲化单位是mm^4; 
I6z=pi*Rou*C_avereff^3*R_wingeff*I6/4;           % 单位是 mg.mm  % I6z = 0.0011;
I_xxam2=I_xxam+I6z*C_maxynd^2;                   % I_xxam2 =6.3082e-004;
% (5) 虚质量气动力矩参数——虚质量转动惯量I7y——绕翅平面下的弦向轴
 I7y=pi*Rou*C_avereff^2*R_wingeff^3*I2/4;     % 单位是 mg.mm^2
%% (a) 弦向压心——展向轴气动力矩/法向气动力% 量纲化需要乘以*C_avereff  or c(r)@r=(R+x_rnd)*r_xcopnd_tr;...
c_zcopnd_tr=I1/F_ndTrans;                            % Y_rcpnd_transaver=0.1715; % *C_avereff  or C(r) @r=R*r_xcopnd_tr;...
c_zcopnd_rot=I2/F_ndRot;                             % Y_rcpnd_rotaver=0.1677;   % *C_avereff  or C(r) @r=R*r_xcopnd_rot;...
% 弦向虚质量力集成的平均压心——对应的是每个片条的中弦点
c_zcopnd_addaver=-0.3359; %c_zcopnd_add=M_xam./F_yadd1;%c_zcopnd_addaver=mean(c_zcopnd_add); 
% c_zcopnd_addtr=I3/F_ndRot;   % c_zcopnd_addtr=0.1587;    % *C_avereff or C(r) @r=R*r_xcopnd_addtr;...
% c_zcopnd_addrot=I4/I6;           % c_zcopnd_addrot=0.4816;  % *C_avereff or C(r) @r=R*r_xcopnd_addrot;...
%% (b) 展向压心——弦向轴气动力矩/法向气动力 % 量纲化需要乘以*R_wingeff 而不是  *(R_wingeff+xr)  
% format long
r_xcopnd_tr=I1/F_ndTrans;   % r_xcopnd_tr =0.788691874094779; % r_xcopnd_tr= 0.7887;          % *R_wingeff
r_xcopnd_rot=I2/F_ndRot;    % r_xcopnd_rot =0.712638285096407;  % r_xcopnd_rot =0.7126;        % *R_wingeff  
% 展向虚质量力集成的平均压心——对应一个特征展向位置
r_xcopnd_addaver=0.7934;  % r_xcopnd_add=M_zadd./F_yadd1;  % r_xcopnd_addaver=mean(r_xcopnd_add); 
% r_xcopnd_addtr=I2/F_ndRot;   % r_xcopnd_addtr=0.7126;     % *R_wingeff
% r_xcopnd_addrot=I3/I6;           % r_xcopnd_addrot=0.5837;   % *R_wingeff
% format short
%% (c) 展向压心——法向轴气动力矩/法向气动力
% r_ycopnd_tr1=I1/F_ndTrans;                          % r_ycopnd_tr= 0.7887;
% r_ycopnd_rot1=I2/F_ndRot;                           % r_ycopnd_rot =0.7126;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 用于Aero_F3_fruitfly_exp & Aero_M_fruitfly2——实验测试的数据
% % 20141122-修改后的输出
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%