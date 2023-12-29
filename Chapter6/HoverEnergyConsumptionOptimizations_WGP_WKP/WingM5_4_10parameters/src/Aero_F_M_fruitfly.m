function F_M=Aero_F_M_fruitfly(x)
% function F_M=Aero_F_M_fruitfly()
% R_wing=3.004;
% C_aver=0.8854;
% xr0=0.3289;
% C_maxyaxis=0.25;
% R_wing=x(1);
% C_aver=x(2);
% xr0=x(3);
% C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=x(1);    T=1/f;          % Hz
phi_m=x(2);
epsilon=x(3);
psi_m=x(4);
zeta=x(5);
psi_0=x(6);
% %%%%%%%%%%%%%%%%%%
R_wing=x(7);
C_aver=x(8);
xr0=x(9);
C_maxyaxis=x(10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis);% 注意调用的函数是wing_shape_fruitfly_sixteen_good2
% (1) 平动环量气动力和力矩参数
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46391;  无量纲, 量纲化单位为mm^4
Coeff_liftdragF_N=wing_para(1,2);  % Coeff_liftdragF_N=0.00682;  %单位是mg*mm  %平动环量法向力
M_xaercoeff=wing_para(1,3);              % M_xaercoeff=0.006038;   %单位是: mg.mm^2   % 平动环量气动力矩参数——绕翅平面下的展向轴
% C_aver1=M_xaercoeff/Coeff_liftdragF_N;  % C_aver1=0.8854;
I1z=wing_para(1,4);                             % I1y=0.016158   % 单位是 mg.mm^2            % 平动环量气动力矩参数——绕翅平面下的弦向轴
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  无量纲, 量纲化单位为mm
M_xrdcoeff=wing_para(1,6);                % M_xrdcoeff=0.0001839; % 单位是mg.mm^2 %转动气动阻尼力矩参数—绕翅平面下的展向轴
% (2) 转动环量气动力和力矩参数
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74847;  无量纲, 量纲化单位为mm^4
F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.003243;  % 单位是 mg.mm   % 转动环量法向力
M_xRotcoeff=wing_para(1,9);             % M_xRotcoeff=0.002871;   % 单位是 mg.mm^2 % 转动环量气动力矩系数——绕翅平面下的展向轴
% C_aver2=M_xRotcoeff/F_yrotcoeff    % C_aver2=0.8854;
I2z=wing_para(1,10);                          % I2y=0.006943;        % 单位是 mg.mm^2            % 转动环量气动力矩参数——绕翅平面下的弦向轴
% (3) 虚质量气动力和力矩参数
I_xzam=wing_para(1,11);                    % I_xzam =0.001424  % 单位是 mg.mm^2  % 虚质量气动力矩参数——绕翅平面下的展向轴
I_xxam=wing_para(1,12);                    % I_xxam =0.000338  % 单位是 mg.mm^2  % 虚质量气动力矩参数——绕翅平面下的展向轴
I5y=wing_para(1,13);                          % I5z=0.0050926   % 单位是 mg.mm       % 虚质量气动力参数—法向力
I6y=wing_para(1,14);                          % I6z=0.00077164    % 单位是 mg.mm       % 虚质量气动力参数—法向力
I7z=wing_para(1,15);                          % I7y=0.0109056;      % 单位是 mg.mm^2   % 虚质量气动力矩参数——绕翅平面下的弦向轴
M_zrdcoeff=wing_para(1,16);             % M_zrdcoeff=0.001169; % 单位是 mg.mm^2 % 转动气动阻尼力矩参数—绕翅平面下的弦向轴
% C_max_LtoT=wing_para(1,17);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=x(1);  T=1/f;   % Hz
% phi_m=x(2);
% epsilon=x(3);    % K=x(3);
% phi_0=x(4);       % eta_m=x(4);
% psi_m=x(5);      % C_eta=x(5);
% zeta=x(6);         % Phi_eta=x(6);
% psi_0=x(7);       % eta_0=x(7);
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,epsilon,psi_m,zeta,psi_0);  %  6个参数的运动学
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0); %  7个参数的运动学
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim();
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
% f=188.7;  % f=234;
% T=1/f;
t=wing_kenimatics(:,1);                % 单位是ms
phi=wing_kenimatics(:,2);        % 拍打角——单位是rad
psi=wing_kenimatics(:,3);            % 拍打角——单位是rad
alpha2=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi=wing_kenimatics(:,5);          % 单位是rad/s
dpsi=wing_kenimatics(:,6);          % 单位是rad/s
ddphi=wing_kenimatics(:,7);       % 单位是rad/s^2
ddpsi=wing_kenimatics(:,8);     % 单位是rad/s^2
C_N=wing_kenimatics(:,9);   
% C_L=wing_kenimatics(:,10);          
% C_D=wing_kenimatics(:,11);   
% C_T=wing_kenimatics(:,12);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀坐标系下的角速率和角加速率——————这组数据来自翅2DOF运动
% 翅膀坐标系下的角速度
omega_x=dpsi;                     % 展向
omega_y=dphi.*sin(psi);       % 法向(初始向左)
omega_z=dphi.*cos(psi);      % 弦向朝上
omega_h=dphi;       % 铰链的角速度% omega_h=-sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
% 翅膀坐标系下的角加速度——用于虚质量力的计算
domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 气动攻角和当地流场速度的计算——取流产速度相对于刚体好了速度，所以整体加负号；
v_y_nonr=-omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;    % v_z=-r*dphi*sin(psi)
% alpha2=atan2(-v_y_nonr,v_z_nonr);   % 正确——注意与下文的alpha=atan2(omega_z,-omega_y)*180/pi; 不同
% % 由于alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %这里atan2给出象限正负值，尤其是alpha>pi/2时
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % 当地来流速度V_nonr=omega_h=dphi;   % 单位是 rad/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 翅坐标系下：法向气动力
% 第一种气动力——翅坐标系下：平动气动环量力——升阻力分量；作用点：片元弦长中点
F_ytran=-sign(alpha2).*V_nonr.^2.*abs(C_N)*Coeff_liftdragF_N*10^(-3);   % 单位是rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
% 第二种气动力—— 翅坐标系下： 旋转气动环量力：法向力；作用点：片元弦长中点
% C_R=1.55;
x_rd0=C_maxyaxis;
C_R=pi*(0.75-x_rd0);
F_yrot=C_R*omega_x.*V_nonr*F_yrotcoeff*10^(-3);      % 单位是rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
% 第三种气动力——翅坐标系下：虚质量气动力：法向力；作用点：片元弦长中点
F_yadd1=(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);  % 修改了原文推出的公式的正负号,且不含/4  % ——第二方案
% 第四种翅膀自身惯性力分量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_wing=2.4*10^-9;                    % kg
% x_com=1.920243385*10^-3;      %  m                           % 质心的展向坐标
% % z_com=-0.149785466+0.636=0.486215*10^-3;      % 到扭转轴的弦向距离
% z_com=0.149785466*10^-3;       %  m                          % 到扭转轴的弦向距离
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis);        %调用函数inertia_moment;
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];
x_com=I_moment(1,1)*10^(-3);    % 质心的展向坐标 % m;
z_com=I_moment(1,2)*10^(-3);    % 到扭转轴的弦向距离 % m;
m_wing=I_moment(1,3)*10^(-3);  % 翅膀形貌学参数变化变化后的质量 % kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % 翅膀自身惯性力—法向
% 总气动力—— 翅坐标系下：总(合成)法向气动力：平动气动环量力+旋转气动环量力+虚质量气动力；作用点：片元弦长中点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
plot(t/T,F_ytran,'r-',t/T,F_yrot,'b-.',t/T,F_inert_y,'m-',t/T,F_yadd1,'c-.','LineWidth',2.5)
legend('平动力','转动力','惯性力','虚质量力')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_N=F_ytran+F_yadd1+F_yrot+F_inert_y;    % 单位是rad^2*s^-2*kg*m=N=10^6uN
% F_N=F_ytran+F_yadd1+F_yrot;    % 单位是rad^2*s^-2*kg*m=N=10^6uN
B=[F_ytran,F_yadd1,F_yrot,F_inert_y];    % size(B)  % (1000*4)     %——mN
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal\optimal_wing_para_motion\Forcenormal_3T.xlsx',B,'sheet1','A1:D1000');
xlswrite('Forcenormal_3T.xlsx',B,'sheet1','A1:D1000'); % 保存输出数据
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：片元法向合力分解到片元垂直方向和水平方向；作用点：片元弦长中点
F_vertical=-sign(alpha2).*F_N.*cos(alpha2);       % 垂直方向vertical, %单位是uN  % 
F_verticalaver=trapz(t,F_vertical)/(3*T);              % F_verticalaver=12.3074uN; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_horizontal=F_N.*sin(abs(alpha2));                   % 水平方向horizontal,%单位是uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%虫体坐标系下：阻力(drag)=F_horiz_y 和侧边力(side force)F_horiz_X
% F_Z=F_vertical;
F_horiz_Y=cos(phi).*F_horizontal;                       % 水平方向_Y_axis推力
% F_horiz_X=sin(phi).*F_horizontal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1=length(F_horiz_Y);
for i=1:N1
    if F_horiz_Y(i)==0
        F_horiz_Y(i)=0.0000001;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%
F_horiz_Y_aver=trapz(t,abs(F_horiz_Y))/(3*T);
Rario_F_vertical_horiz_aver=F_verticalaver/F_horiz_Y_aver    % 升推比
F_vertical_to_horiz_rms=sqrt(trapz(t,(abs(F_vertical)./(abs(F_horiz_Y))-Rario_F_vertical_horiz_aver).^2)/(3*T)) % 升推比rms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 法向合力分解到垂直方向和水平方向的分量力随时间的变化规律与实验测试结果的对比
figure(10)
% plot(t/T,F_vertical,'r-',t/T,-sign(alpha2).*F_horiz_Y,'b-.','LineWidth',2.5)
plot(t/T,F_vertical,'r-',t/T,F_horiz_Y,'b-.','LineWidth',2.5)
xlabel('Normalized time')
ylabel('Force  ( \muN )')
legend('\itF_{\itvertical,Z} ','\itF_{\ithorizontal,Y}')
title('法向合力分解到垂直方向和水平方向的分量力随时间的变化规律')
grid on
t=t';
% axis([t(1,1)/T,t(1,length(t))/T,min(F_horiz_Y)-0.5,max(F_vertical)+0.5])
axis([t(1,1)/T,t(1,1)/T+1,min(F_horiz_Y)-0.5,max(F_vertical)+0.5])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 气动功率的计算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一模块——展向扭转轴——气动力矩分量%%%%
%%%%含七部分, 分别是：%%%%%%%%%%%%%%单位: (mN.mm)或(uN.m)
%%%%第一部分——平动环量产生的气动力矩%%%%%%
%%%%第二部分——转动环量产生的气动力矩%%%%%%
%%%%第三部分——转动气动阻尼力矩%%%%%%%%%
%%%%第四部分——转动虚质量力矩%%%%%%%%%%
%%%%第五部分——扭转铰链的弹性回复力矩%%%%%%
%%%%第六部分——翅膀重力矩%%%%%%%%%%%%
%%%%第七部分——翅膀自身惯性力矩%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一部分——平动环量产生的气动力矩—扭转轴—单位: (mN.mm)或(uN.m)
%% (3) 方案3——针对旋转轴气动力矩——变压心位置
k_xaero=1; 
N=length(t);
M_xtrans=zeros(N,1);
Y_rcpnd_trans=zeros(N,1);
for i=1:1:N
    Y_rcpnd_trans(i,1)=COP_Ycpnd2_TransCirc(alpha2(i,1),R_wing,xr0,C_maxyaxis); %调用函数COP_Ycpnd2_TransCirc求解净压心的无量纲位置Y_rcpnd; % 正负交替
    %  Y_rcpnd_trans=abs(Y_rcpnd_trans(i,1));  % 正
    % 下面的单位是 (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
    %  M_xaercoeff=0.0060;   %单位是: mg.mm^2
    % 平动环量产生的旋转轴气动力矩一开始是逆时针的
     % M_xtrans(i,1)=-k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);      
     M_xtrans(i,1)=-k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);  
end
%% 第二部分——转动环量产生的气动力矩——扭转轴
k_xRot=1;  
% C_R=1.55;    % 该系数可以修改
% x_rd0=C_maxyaxis;
% C_R=pi*(0.75-x_rd0);
M_xRotcoeff=k_xRot*C_R*M_xRotcoeff;
N=length(t);
M_xrotcirc=zeros(N,1);
Y_rcpnd_rot=zeros(N,1);
for i=1:1:N
    % 转动环量绕扭转轴气动力矩——调用函数COP_Ycpnd2_RotCirc求解净压心的无量纲位置Y_rcpnd_rot
    % 转动环量产生的—压心分布符合Dickinson函数 or 压心在中弦点 or 压心在c(r)/4处—扭转轴力矩
    Y_rcpnd_rot(i,1)=COP_Ycpnd2_RotCirc(alpha2(i,1),R_wing,xr0,C_maxyaxis); % 压心分布符合Dickinson函数: 
   % Y_rcpnd_rot=abs(Y_rcpnd_rot(i,1));
   % 下面的单位是 (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
   %  M_xrotcirc(i,1)=k_xRot*C_R*sign(alpha2(i,1)).*omega_x(i,1).*abs(V_nonr(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3); % 转动环量绕扭转轴气动力矩  
   M_xrotcirc(i,1)=k_xRot*C_R*omega_x(i,1).*abs(omega_h(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3);    % 转动环量绕扭转轴气动力矩
end
%% 第三部分——转动气动阻尼力矩—扭转轴—单位: (mN.mm)或(uN.m)
C_RD=1;  
M_xrd=-C_RD*omega_x.*abs(omega_x)*M_xrdcoeff*10^(-3);   % 阻尼力矩一开始是顺时针的 % 注意这个方向取决于omega_x
%% 第四部分——转动虚质量力矩+平动虚质量力矩—扭转轴—单位: (mN.mm)或(uN.m)
k_am=1;   % 虚质量力矩系数; 2.35是近似上限； psi_max =32.9386
% M_xam=-k_am*(-I_xzam*(domega_z+omega_x.*omega_y)+I_xxam*domega_x)*10^(-3);   % 初始逆时针(-)
M_xam=-k_am*(-I_xzam*(domega_z+omega_x.*omega_y)+I_xxam*domega_x)*10^(-3);   % 初始逆时针(-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(101)
plot(t/T,M_xtrans,'r--',t/T,M_xrotcirc,'b-',t/T,M_xrd,'g--',t/T,M_xam,'c-','LineWidth',2)
legend('平动力矩','转动力矩','转动阻尼力矩','虚质量力矩')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 绕扭转轴的功率计算
% 平动环量力矩功率  % uW
P_xtrans=-M_xtrans.*omega_x*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% 转动阻尼力矩功率  % uW
P_xrd=-M_xrd.*omega_x*10^-3;  
% 转动环量力矩功率  % uW  
P_xrotcirc=-M_xrotcirc.*abs(omega_x)*10^-3; 
% 虚质量力矩功率和翅膀自身惯性力矩功率  % uW
P_xam=-M_xam.*omega_x*10^-3;
P_aerox=P_xtrans+P_xrd+P_xrotcirc+P_xam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(102)
plot(t/T,P_xtrans,'r--',t/T,P_xrotcirc,'b-',t/T,P_xrd,'g--',t/T,P_xam,'c-',t/T,P_aerox,'k-','LineWidth',2)
legend('平动功率','转动功率','转动阻尼功率','虚质量功率','扭转轴总功率')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二模块——弦向转动轴——气动力矩分量%%%%%%
%%%%含七部分, 分别是：%%%%%%%%%%%%%%单位: (mN.mm)或(uN.m)
%%%%第一部分——平动环量产生的气动力矩%%%%%%%
%%%%第二部分——转动环量气动力矩%%%%%%%%%%%
%%%%第三部分——转动虚质量力矩%%%%%%%%%%%%
%%%%第四部分——转动阻尼力矩%%%%%%%%%%%%%
%%%%第五部分——拍打轴转动铰链的弹性回复力矩%%%%
%%%%第六部分——拍打轴驱动力矩%%%%%%%%%%%%
%%%%第七部分——翅膀自身惯性力矩%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一部分 平动环量法向气动力产生的力矩——绕翅平面下的弦向轴
k_ztrans=1;  
M_ztrans=k_ztrans*sign(alpha2).*(I1z.*C_N.*omega_h.^2)*10^(-3);   % I1z=0.0162; % 单位是 mg.mm^2=10^-3uN*mm=10^-6 mN*mm
%% 第二部分 转动环量法向气动力产生的力矩——绕翅平面下的弦向轴
C_R=1.55;    % 该系数有问题
k_zrot=1;
M_zrot=-k_zrot*I2z*C_R*omega_x.*omega_h*10^(-3);  % I2z=0.0069; 单位是: mg.mm^2*(rad*s^-1)^2=10^-6mN.mm=10^-3uN.mm
%% 第三部分 虚质量法向气动力产生的力矩——绕翅平面下的弦向轴 % I5y和I6y单位是 mg.mm=*10^(-9) kg.m  %单位是10^6uN
% k_za=0.35; 
k_za=1;  % 该系数有问题么？
% M_zadd=-k_za*(-I7z*(domega_z+omega_x.*omega_y)+I_xzam*domega_x)*10^(-3); % 原文推出的公式——更为合理些哦   % 初始逆时针(-)
M_zadd=k_za*(-I7z*(domega_z+omega_x.*omega_y)+I_xzam*domega_x)*10^(-3); 
%% 第四部分 转动阻尼力矩——绕翅平面下的弦向轴
% C_RD2=0.05*C_RD;
C_RD2=1*C_RD;  % 该系数有问题么？
M_zrd=-C_RD2*omega_x.*abs(omega_x)*M_zrdcoeff*10^(-3);   % 阻尼力矩一开始是顺时针的 % 注意这个方向取决于omega_x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103)
plot(t/T,M_ztrans,'r--',t/T,M_zrot,'b-',t/T,M_zrd,'g--',t/T,M_zadd,'c-','LineWidth',2)
legend('平动力矩','转动力矩','转动阻尼力矩','虚质量力矩')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 拍打轴的功率
% % 平动环量力矩功率  % uW
% P_ztrans=M_ztrans.*omega_z*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % 转动环量力矩功率  % uW  
% P_zrotcirc=M_zrot.*omega_z*10^-3; 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % 虚质量力矩功率和翅膀自身惯性力矩功率  % uW
% P_zam=M_zadd.*omega_z*10^-3;
% % P_inert_z=M_inert_z.*omega_z*10^-3;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % 转动阻尼力矩功率  % uW
% P_zrd=M_zrd.*omega_z*10^-3;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 平动环量力矩功率  % uW
P_ztrans=cos(psi).*M_ztrans.*omega_z*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
%%%%%%%%%%%%%%%%%%%%%%%%%
% 转动环量力矩功率  % uW  
P_zrotcirc=cos(psi).*M_zrot.*omega_z*10^-3; 
%%%%%%%%%%%%%%%%%%%%%%%%%
% 虚质量力矩功率和翅膀自身惯性力矩功率  % uW
P_zam=cos(psi).*M_zadd.*omega_z*10^-3;
% P_inert_z=M_inert_z.*omega_z*10^-3;
%%%%%%%%%%%%%%%%%%%%%%%%%
% 转动阻尼力矩功率  % uW
P_zrd=cos(psi).*M_zrd.*omega_z*10^-3;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_aeroz=P_ztrans+P_zrd+P_zrotcirc+P_zam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(104)
plot(t/T,P_ztrans,'r--',t/T,P_zrotcirc,'b-',t/T,P_zrd,'g--',t/T,P_zam,'c-',t/T,P_aeroz,'k-','LineWidth',2)
legend('平动功率','转动功率','转动阻尼功率','虚质量功率','拍打轴总功率')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀两自由度运动的功率计算——平动功率和扭转功率——扭转轴功率和拍打轴功率
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aerox=P_xtrans+P_xrd+P_xrotcirc+P_xam;
% P_totalx=P_xtrans+P_xrd+P_xrotcirc+P_xam+P_inert_x;
P_psix_total=P_aerox;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aeroz=P_ztrans+P_zrd+P_zrotcirc+P_zam;
% P_totalz=P_ztrans+P_zrd+P_zrotcirc+P_zam+P_inert_z;
P_phiz_total=P_aeroz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis);        %调用函数inertia_moment;
% uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];
Ix_inertia=I_moment(1,4)*10^(-3);     % 由平移轴定理获得  % g.mm^2=10^(-9)kg.m^2——for 10^-3uW
Iz_inertia=I_moment(1,5)*10^(-3);     % 由平移轴定理获得  % g.mm^2=10^(-9)kg.m^2——for 10^-3uW
% 下面——rad*s^-1*10^(-9)kg.m^2*rad*s^-2=10^(-9)(N=kg*m.s^-2)*(m.s^-1)=10^-9W=10^-6mW=10^-3uW
P_Ix_inertia=omega_x.*Ix_inertia.*domega_x;      % uW
P_Iz_inertia=omega_z.*cos(psi).*Iz_inertia.*domega_z;   % uW     % cos(psi).*
P_totalx=P_psix_total+P_Ix_inertia;
P_totalz=P_phiz_total+P_Iz_inertia;
% P_total=[P_totalx,P_totalz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aero_aver=trapz(t,(P_psix_total+P_phiz_total))/(3*T)
% P_inertia_aver=trapz(t,(P_Ix_inertia+P_Iz_inertia))/(3*T)
% P_aero_aver_to_inertia_aver=P_aero_aver/P_inertia_aver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_total=Aero_M_fruitfly2_exp(x);      % 调用函数——正功率——平动功率和扭转功率——扭转轴功率和拍打轴功率
% % size(P_total)              % (1000*2)
% P_totalx=P_total(:,1);    % (1000*1)
% P_totalz=P_total(:,2);    % (1000*1)
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total.xlsx',P_total,'sheet1','A1:B2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(P_psix_total);
for i=1:N
    if P_psix_total(i)<=0
        P_psix_total(i)=0;
    end
end
P_psix_total_posi=P_psix_total;% 正的扭转气动总功率
%%%%%%%%%%%%%%%%%%%
N=length(P_Ix_inertia);
for i=1:N
    if P_Ix_inertia(i)<=0
        P_Ix_inertia(i)=0;
    end
end
P_Ix_inertia_posi=P_Ix_inertia;% 正的扭转惯性功率
%%%%%%%%%%%%%%%%%%%
N=length(P_totalx);
for i=1:N
    if P_totalx(i)<=0
        P_totalx(i)=0;
    end
end
P_totalx_posi=P_totalx; % 正的扭转总功率=气动总功率+惯性功率
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:N  % N=length(P_totalx);
    if P_phiz_total(j)<=0
        P_phiz_total(j)=0;
    end
end
P_phiz_total_posi=P_phiz_total;% 正的拍打气动总功率
%%%%%%%%%%%%%%%%%%%
for j=1:N  % N=length(P_totalx);
    if P_Iz_inertia(j)<=0
        P_Iz_inertia(j)=0;
    end
end
P_Iz_inertia_posi=P_Iz_inertia;% 正的拍打惯性功率
%%%%%%%%%%%%%%%%%%%
for j=1:N  % N=length(P_totalx);
    if P_totalz(j)<=0
        P_totalz(j)=0;
    end
end
P_totalz_posi=P_totalz;% 正的拍打总功率=气动总功率+惯性功率
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_aero_aver_posi=trapz(t,(P_psix_total_posi+P_phiz_total_posi))/(3*T);
P_inertia_aver_posi=trapz(t,(P_Ix_inertia_posi+P_Iz_inertia_posi))/(3*T);
P_aero_aver_to_inertia_aver=P_aero_aver_posi/P_inertia_aver_posi       % 全正气动功率与全正惯性功率的比值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 针对最优翅膀形貌参数的翅扭转功率、拍打功率和总功率输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
plot(t/T,P_totalx_posi,'r--',t/T,P_totalz_posi,'b-','LineWidth',2)
hold on
P_total=P_totalx_posi+P_totalz_posi;
plot(t/T,P_total,'k-','LineWidth',2.5)
xlabel('Normalized time')
ylabel('Power output ( \muW )')
legend('\itP_{\itx,\psi}','\itP_{\itZ,\phi}','\itP_{\ittotal}')
title('针对最优翅膀形貌参数的翅扭转功率、拍打功率和总功率输出')   
grid on
% axis([t(1,1)/T,t(1,length(t))/T,0,max(P_total)+0.5])
axis([t(1,1)/T,t(1,1)/T+1,0,max(P_total)+0.5])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_total_posi=[P_totalx_posi,P_totalz_posi]; 
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total_positive.xlsx',P_total_posi,'sheet1','A1:B1000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=188.7; 
% f=x(1);
% T=1/f;  %翅拍频率 (Hz)和周期  % w =1185.6; 
% t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1——XXXXXX
% t_00= 1/(4*f);    % t_00= -1/(4*f);
% t=linspace(t_00,t_00+3*T,1000); 
P_totalx_aver=trapz(t,P_totalx_posi)/(3*T);  % 时均正功率——平动功率和扭转功率——扭转轴功率和拍打轴功率
P_totalz_aver=trapz(t,P_totalz_posi)/(3*T);  % 时均正功率——平动功率和扭转功率——扭转轴功率和拍打轴功率
P_total_aver=P_totalx_aver+P_totalz_aver;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_M=[F_verticalaver,P_total_aver];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%