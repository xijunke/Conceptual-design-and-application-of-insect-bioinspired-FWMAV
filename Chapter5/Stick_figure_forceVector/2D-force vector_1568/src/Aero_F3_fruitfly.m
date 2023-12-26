%% 气动力: Solution of the aerodynamic force
clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用含翅形貌参数化和气动力系数的函数
% wing_para_output=[R_wingeff,C_avereff,F_nd,Y_rcpnd,Y_rnd,I_xxam,I_xyam,I3,I3y,I4y];
wing_para=wing_shape_fruitfly();      %调用函数wing_shape;  
% size(wing_para)
r2_nd=0.5801; 
R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;        % 单位是 mm
C_avereff=wing_para(1,2);  % C_avereff =0.8854;          % 单位是 mm
F_nd=wing_para(1,3);         % F_nd = 0.4639;                % 量纲单位是mm^4
Y_rnd=wing_para(1,4);        % Y_rnd =0.1402;               % 量纲单位是mm
% I_xyam=wing_para(1,5);      % I_xyam =-0.0012;            % 单位是 mg.mm^2;  
% I_xxam=wing_para(1,6);      % I_xxam =2.5080e-004;    % 单位是 mg.mm^2;  
I3=wing_para(1,7);              % I3 =0.7485;                      % 无量纲，量纲化单位是 mm^4
I3y=wing_para(1,8);            % I3y = 0.0051;                    % 单位是 mg.mm
I4y=wing_para(1,9);            % I4y =-5.3461e-004;          % 单位是 mg.mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',dalpha',ddphi',ddpsi',ddalpha',C_L',C_D',C_N1',C_N2',C_N3',C_T',alpha2'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly();      %调用函数kenimatics_wing_and_AoA;  % size(wing_kenimatics)
size(wing_kenimatics)
t=wing_kenimatics(:,1);                % 单位是ms
phi=wing_kenimatics(:,2);        % 拍打角——单位是rad
psi=wing_kenimatics(:,3);            % 拍打角——单位是rad
alpha1=wing_kenimatics(:,4);      % alpha1=pi/2+psi.*sign(dphi);    %几何攻角——弧度制             %输出———全正几何攻角
alpha2=wing_kenimatics(:,5);      % alpha2=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi=wing_kenimatics(:,6);          % 单位是rad/s
dpsi=wing_kenimatics(:,7);          % 单位是rad/s——可能会出问题—alpha—XXXXXX
% dalpha1=wing_kenimatics(:,8);      % 单位是rad/s——%输出——有正有负
ddphi=wing_kenimatics(:,9);       % 单位是rad/s^2
ddpsi=wing_kenimatics(:,10);     % 单位是rad/s^2——可能会出问题—alpha—XXXXXX
% ddalpha1=wing_kenimatics(:,11);  % 单位是rad/s^2——%输出——有正有负
C_L=wing_kenimatics(:,12);          % 可能会出问题—alpha—XXXXXX
C_D=wing_kenimatics(:,13);         % 可能会出问题—alpha—XXXXXX
C_N1=wing_kenimatics(:,14);   % 可能会出问题—alpha—XXXXXX
C_N2=wing_kenimatics(:,15);       % 可能会出问题—alpha—XXXXXX
% C_N3=wing_kenimatics(:,16);   % 可能会出问题—alpha—XXXXXX
% C_T=wing_kenimatics(:,17);     % 可能会出问题—alpha—XXXXXX
alpha=alpha2;
% C_N=C_N2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀坐标系下的角速率和角加速率——————这组数据来自翅2DOF运动
% 翅膀坐标系下的角速度
omega_x=-dpsi; 
omega_y=dphi.*sin(psi); 
omega_z=dphi.*cos(psi);
omega_h=dphi;       % 铰链的角速度% omega_h=-sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
% 翅膀坐标系下的角加速度——用于虚质量力的计算
domega_x=-ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一种气动力——翅坐标系下：平动气动环量力——升阻力分量；作用点：片元弦长中点
%% 输入参数为：4个
% C_avereff;  R_wingeff;    F_nd;                ——————这组数据来自翅形貌参数化;% F_nd=0.50236;
% omega_h;  C_L(alpha);  C_D(alpha);       ——————这组数据来自翅2DOF运动铰链角速率, 攻角
%%
Rou=1.225;        %这里单位是Kg/m^3——可换算成——10^6/(10^3)^3=10^-3mg/mm^3    
g=9.821/10^6;   % 这里重力加速度:g=9.821N/kg=9.821/10^6=9.821e-006   N/mg  ——g的国际单位是m/s^2或N/kg
LD_aerocoeff=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);  % 单位是：kg/m^3*mm*m^3= 10^(-12) kg.m
lift_wingframe=sign(alpha).*omega_h.^2.*C_L*LD_aerocoeff/g;      % 单位是N / (N/mg) =mg
drag_wingframe=omega_h.^2.*C_D*LD_aerocoeff/g;
f=188.7;  T=1/f;          % Hz
Phi=max(phi);              % Phi*180/pi=73.4427
S=0.0266*100*10^(-6);      % m^2
U_ref=2*f*Phi*r2_nd*R_wingeff*10^(-3);   % m/s   % U_ref =0.8430
C_L=(lift_wingframe*g)/(Rou*S*U_ref^2);
figure(30)
plot(t/T,C_L,'r-');   
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对瞬时气动升阻力使用梯形积分函数trapz进行数值积分，除以周期, 求解平均升阻力
lift_aver=trapz(t,lift_wingframe)/T                                    % 
drag_aver=trapz(t,drag_wingframe)/T                             %这里的阻力也是去了绝对值的！
drag_instabs_aver=trapz(t,abs(drag_wingframe))/T;        %这里的阻力也是去了绝对值的！
lift2drag_aver=lift_aver/drag_instabs_aver                 %升阻比
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%输出结果:
% lift_aver =1.0913
% drag_aver =0.1573
% lift2drag_aver =1.1037
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)
F_LD=plot(t/T,lift_wingframe,'r-',t/T,drag_wingframe,'k-');   
xlabel('\itNormalized time')
ylabel('\itAerodynamic force (mg)')                    %单位换算到mg了哦
title('Aerodynamic lift and drag force \itvs. \itt \rm for flapping wing')
% legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\it|F_D|)')
legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\itF_D)')
set(F_LD,'LineWidth',2)
grid on
% axis([-inf,inf,-900,900])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：平动气动环量力(1)——法向力；作用点：片元弦长中点
F_ny_coeff=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);        % 单位是kg/m^3*m*m^3=kg.m
F_ny1=omega_h.^2.*C_N1*F_ny_coeff*10^6;  % 单位是uN
F_ny2=omega_h.^2.*C_N2*F_ny_coeff*10^6;  % 单位是uN————这个更靠谱
figure(14)
plot(t/T,F_ny1,'r-',t/T,F_ny2,'g-')
xlabel('\itNormalized time')
ylabel('\itF_{ny_1} and F_{ny_2} (uN)')
legend('F_{ny_1}','F_{ny_2}')
title('平动气动环量力(法向)随时间的变化规律')   % 平动气动环量力(法向)随时间的变化规律
grid on
% 平动气动环量力/平动法向力分量——法向气动力系数
% C_Nt=F_ny1./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %单位是uN
% figure(15)                                                       
% plot(t/T,C_Nt,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt')
% legend('C_Nt')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：平动气动环量力(2)——法向力和切向力；作用点：片元弦长中点
% % 无量纲气动力分量F_nd=R2nd2+2*xr_nd*R1nd1+xr_nd^2;   在翅根到翅肩的无量纲距离xr_nd=0时; F_nd=R2nd2
% R2_nd2=0.2350;      % 由无量纲气动力分量F_nd=R2nd2计算气动力的结果不够准确
% F_ny_total=-1/2*Rou*R_wingeff^3*C_avereff*R2_nd2*C_N3.*dphi.*abs(dphi)*(10^-3*10^-9*10^3);                    %法向，单位是mN
% F_Ttz=-1/2*Rou*R_wingeff^3*C_avereff*R2_nd2*C_T.*dphi.*abs(dphi).*psi.*abs(psi)*(10^-3*10^-9*10^3);  %弦向，单位是mN
% figure(16)
% F_LD2=plot(t/T,F_ny_total,'r-',t/T,F_Ttz,'k-');   
% xlabel('\itNormalized time')
% ylabel('\itAerodynamic force (mN)')
% title('Aerodynamic normal and tangential force \itvs. \itt \rm for flapping wing')
% legend('\itinstantaneous normal force (\itF_N)','\itinstantaneous tangential force (\itF_T)')
% set(F_LD2,'LineWidth',2)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二种气动力——翅坐标系下：虚质量气动力：法向力；作用点：片元弦长中点
%%% 拍打和被动扭转的协同作用？——虚质量力——参考SK Agrawal的文献推导过程，以便理解其物理意义
%% 输入参数为：10个
% I3y;   I4y;                                                         ——————这组数据来自翅形貌参数化
% domega_z; omega_x; omega_y;                     ——————这组数据来自翅2DOF运动翅角加速率, 角速率
% dphi; dalpha; ddphi; ddalpha; alpha;            ——————这组数据来自翅2DOF运动角速率, 角加速率, 攻角
%%  (1) 2010_JFM_Aeromechanics-虚质量气动力公式——% I3y和I4y单位是 mg.mm-换算到kg.m需要*10^(-9)
% F_yam=sign(phi).*(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x)*(10^(-9)*10^6);   %单位是uN
% F_yam1=(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x)*(10^(-9)*10^6); % 原文推出的公式
% 注意公式前(-)符——根据果蝇翅膀坐标系修改了角速度和角加速度
F_yam1=-(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x/4)*10^(-3);    %单位是uN
%% (2)下面是2001-JEB-Sanjay SP-虚质量气动力公式
% F_yam2=(I3y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I4y/4)*ddalpha/4)*10^(-3);% 原文推出的公式——不动
% 下面的alpha1为全正几何攻角，有正有负的dalpha1和ddalpha1
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))-(I4y/4)*ddalpha1/4)*10^(-3); 
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))+(I4y/4)*ddalpha1/4)*10^(-3); 
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))-(I4y/4)*ddpsi/4)*10^(-3); %这里的dpsi和ddpsi方向变化需要确定
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))+(I4y/4)*ddpsi/4)*10^(-3); %这里的dpsi和ddpsi方向变化需要确定
% 下面的alpha2有正有负, dalpha2和ddalpha2需要求导
% F_yam2=-(I3y*(ddphi.*sin(alpha2)+dphi.*dalpha2.*cos(alpha2))-(I4y/4)*ddalpha2/4)*10^(-3); 
i=length(alpha2);   % size(alpha2)=(2000*1)
dalpha2=[diff(alpha2);alpha2(i,1)];    
j=length(dalpha2);
ddalpha2=[diff(dalpha2);dalpha2(j,1)];  % size(ddalpha2)
F_yam2=-(I3y*(ddphi.*sin(abs(alpha2))+dphi.*dalpha2.*cos(abs(alpha2)))-(I4y/4)*ddalpha2/4)*10^(-3); 
figure(17)          % 该气动力有些小
hold on
plot(t/T,F_yam1,'r-',t/T,F_yam2,'b-')                    % 注意这里的domega_x=ddpsi与ddalpha的符号存在不同哦
xlabel('\itNormalized time')
ylabel('\itF_{am_y} (uN)')
legend('F_{am1_y}','F_{am2_y}')
title('虚质量气动力(法向)随时间的变化规律')          % 虚质量气动力(法向)随时间的变化规律
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-');     %画x-axis
grid on
%% 虚质量气动力/平动法向力分量——法向气动力系数
% C_Na=F_yam2./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %单位是uN
% figure(18)                                                       
% plot(t/T,C_Na,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Na')
% legend('C_Na')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第三种气动力—— 翅坐标系下： 旋转气动环量力：法向力；作用点：片元弦长中点
%% 输入参数为：3个
% C_avereff;   R_wingeff;   I3;                              ——————这组数据来自翅形貌参数化
% omega_h;  omega_x;                                       ——————这组数据来自翅2DOF运动铰链角速率, 翅角速率
%%
C_R=1.55;
% F_zrot=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*omega_x.*omega_h*I3*(10^(-12)*10^3);%单位是mN——未考虑方向
F_y_rotcoeff=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*I3*10^(-12);          % 单位是kg/m^3*mm^2*mm^2= kg.m 10^(-12)
% F_y_rot=sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6;    % 单位是uN
F_y_rot=omega_x.*omega_h*F_y_rotcoeff*10^6;    % 单位是uN
figure(19)                  % 该气动力有些小  
plot(t/T,F_y_rot,'b-')
xlabel('\itNormalized time')
ylabel('\itF_{y,rot} (uN)')
legend('\itF_{y,rot}')
title('旋转气动环量力(法向)随时间的变化规律')   % 旋转气动环量力(法向)随时间的变化规律
grid on
%% 旋转气动环量力/平动法向力分量——法向气动力系数               %  求解转动环量气动力系数意义不大价值
% figure(20)                                                                       
% % C_Nr=F_y_rot./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %单位是uN
% C_Nr=F_y_rot./(sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6);  %单位是mN
% plot(t/T,C_Nr,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nr')
% legend('C_Nr')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 各种机制的法向气动力系数以及总和
% C_Ns=C_Nt+C_Na+C_Nr;  %单位是mN
% figure(21)                                                                         %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Nt,'r-',t/T,C_Na,'g-',t/T,C_Nr,'b-',t/T,C_Ns,'k-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt & C_Na & C_Nr & C_Ns')
% legend('C_Nt','C_Na','C_Nr','C_Ns')
% grid on
%% 法向气动力系数总和分解到垂直于来流方向(C_L)以及平行于来流方向(C_D)
% C_Lv=sign(dphi).*C_Ns.*sin(psi);     % 垂直方向vertical气动力系数
% C_Dh=sign(dphi).*C_Ns.*cos(psi);    % 水平方向horizontal气动力系数
% figure(22)                                                                       %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Nt,'r-',t/T,C_Lv,'g-',t/T,C_Dh,'b-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt & C_Lv & C_Dh')
% legend('C_Nt','C_Lv','C_Dh')
% title('法向气动力系数总和分解到垂直于来流方向(C_L)以及平行于来流方向(C_D)') 
% grid on
% figure(23)  % 法向气动力系数总和分解到垂直于来流方向(C_L)                                                                       %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Lv,'g-')
% xlabel('\itNormalized time')
% ylabel('\itC_Lv')
% legend('C_Lv')
% title('法向气动力系数总和分解到垂直于来流方向(C_L)')   
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 总气动力—— 翅坐标系下：总(合成)法向气动力：平动气动环量力+旋转气动环量力+虚质量气动力；作用点：片元弦长中点
F_N=F_ny1+F_y_rot+F_yam1;    % 法向平动环量力F_ny由平动升阻力系数合成而得,%单位是mN
% F_N=F_ny_total+F_y_rot+F_yam;  % 法向平动环量力F_ny_total由平动法向气动力系数计算而得,%单位是mN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 将单个稳定的周期的总(合成)法向气动力数据输出到Excel表——mN
% 以备thr_dim_chord2程序的调用用于绘制球棍图
B=F_N;    % size(B)  % (2000*1)     %——mN
xlswrite('Forcenormal_oneT.xlsx',B,'sheet1','A1:A2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(24)
plot(t/T,F_N,'r-')
xlabel('\itNormalized time')
ylabel('\itF_N (mN)')
legend('F_N')
title('平动气动环量力+旋转气动环量力+虚质量气动力(法向)随时间的变化规律')
grid on
%% 翅坐标系下：片元法向合力分解到片元垂直方向和水平方向；作用点：片元弦长中点
F_vertical=abs(F_N.*cos(alpha));                          % 垂直方向vertical, %单位是mN
F_horizontal=sign(alpha).*F_N.*sin(alpha);   % 水平方向horizontal,%单位是mN
% F_vertical=abs(F_N.*cos(alpha));                          % 垂直方向vertical, %单位是mN
% F_horizontal=abs(sign(alpha).*F_N.*sin(alpha));                                   % 水平方向horizontal,%单位是mN
% F_v=(abs(F_N).*cos(alpha))*10^-3/g;     % 垂直方向vertical,%单位是mg
% F_h=(F_N.*sin(alpha))*10^-3/g;             % 水平方向horizontal,%单位是mg
F_vaver=trapz(t,F_vertical)/T                       % T——这里有问题哦
F_haver=trapz(t,F_horizontal)/T                       %这里的阻力也是取了绝对值的！
F_haverabs=trapz(t,abs(F_horizontal))/T;         %这里的阻力也是取了绝对值的！
F_v2haver=F_vaver/F_haverabs          %升阻比
%%%output: Units: mg
% F_vaver =10.8146
% F_haver =1.5524
% F_v2haver =0.9314
figure(25)
plot(t/T,F_vertical,'r-',t/T,F_horizontal,'b-')
xlabel('\itNormalized time')
ylabel('\itF_{vertical} & F_{horizontal} (mN)')
legend('F_{vertical}','F_{Phorizontal}')
title('法向合力分解到垂直方向和水平方向的分量力随时间的变化规律')   
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下： 由法向气动力系数求解: 垂直方向和水平方向气动力系数———这里的求解不正确
% C_N=F_N./((1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));   % 单位有问题哦XXXX
C_N=F_N./(-(1/2*Rou*omega_h.^2.*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));   % 单位有问题哦XXXX
C_v=F_vertical./(-sign(dphi).*(1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3)); 
C_h=F_horizontal./(-sign(dphi).*(1/2*Rou*sign(dphi).*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));
% % C_h=F_horizontal./((1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3)); 
% C_v=C_N.*cos(alpha); 
% C_h=C_N.*sin(alpha); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对瞬时气动力系数使用梯形积分函数trapz进行数值积分，除以周期, 求解平均垂直方向和水平方向气动力系数
C_vaver=trapz(t,C_v)/T                     
C_haver=trapz(t,C_h)/T                      %这里的阻力也是取了绝对值的！
C_habsaver=trapz(t,abs(C_h))/T;        %这里的阻力也是取了绝对值的！
C_v2haver=C_vaver/C_habsaver              %升阻比
% 输出
% C_vaver = 3.7772e+005
% C_haver =4.5096e+006
% C_v2haver =0.0540
% figure(26)
% plot(t/T,C_N,'g-')
% xlabel('\itNormalized time')
% ylabel('\itC_N')
% legend('C_N')
% title('法向合气动力系数随时间的变化规律')
% grid on
% figure(27)
% plot(t/T,abs(C_v),'r-',t/T,C_h,'b-')
% xlabel('\itNormalized time')
% ylabel('\itC_v & C_h')
% legend('C_v','C_h')
% title('垂直方向和水平方向气动力系数随时间的变化规律')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

