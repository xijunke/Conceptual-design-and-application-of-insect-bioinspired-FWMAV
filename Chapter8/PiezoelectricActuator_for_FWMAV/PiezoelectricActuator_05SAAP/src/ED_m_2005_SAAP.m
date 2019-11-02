% 计算压电驱动器输出的尖端位移和堵死力
% 针对不同的尺寸外形和直流驱动电压值的情况来进行具体的计算
clear all;clc;
%% 碳纤维的厚度
t_cf=40;  % um
% t_cf=60;  % um
%% (1) 第一种尺寸%%%%%%%%%%%%
l_pzt=8;
l_ext=4;
w0=3.5;
w1=1;
%% (2) 第二种尺寸%%%%%%%%%%%%
% l_pzt=10;
% l_ext=10;
% w0=2.5;
% w1=1.5;
%%%%%%%%%%%%%%
t_cf=t_cf*(10^-6);
l_pzt=l_pzt*10^-3;
l_ext=l_ext*10^-3;
w0=w0*10^-3;
w1=w1*10^-3;
%%%%%%%%%%%%%%
%% pzt——压电陶瓷的参数
t_pzt=127e-6;    % um
Rou_pzt=7800;   %密度kg/m
U=250; % V
% U=300; % V
E3=U/t_pzt;            %外加驱动的电场强度v
E=6.2e10;   v=0.31;%弹性和剪切模量Gpa &泊松比 其中G12=24e9
%E=E/(1-v^2);
d31=-320e-12; d32=d31; %压电常数pm/v
d=[d31;d32;0];
%% cf——碳纤维的力学性能参数
Rou_cf=1500;   %密度kg/m
E1=300e9;  E2=7e9; G12=5e9;   %弹性和剪切模量Gpa
v12=0.35; v21=E2*v12/E1;%泊松比
%% gf——玻璃碳纤维的力学性能参数
t_gf=120e-6;     % um
Rou_gf=1600;   %密度kg/m
l_gf=l_ext;
%%%%%%%%%%%%%%%%%%%%%%%%%
wnom=(w0+w1)/2;
l_r=l_ext/l_pzt;
w_r=w0/wnom;
m_actuator=(2*Rou_pzt*t_pzt*l_pzt*wnom+Rou_cf*t_cf*l_pzt*wnom+Rou_cf*t_cf*l_ext*w1+2*Rou_gf*t_gf*l_gf*w1)
d_r=(Rou_cf*t_cf+2*Rou_gf*t_gf)/(2*Rou_pzt*t_pzt+Rou_cf*t_cf)
%%%%%%%%%%%%%%%%%%%%%%%%%
%% 刚度矩阵的元素   各同异性材料PZT
Q112=E/(1-v^2);  Q222=Q112;
Q122=E*v/(1-v^2);
Q662=24E9; G=Q662;
%Q662=E/(2*(1+v)); G=Q662;
%刚度矩阵的元素   各向异性材料CF
Q11=E1/(1-v12*v21);
Q12=v12*E2/(1-v12*v21);
Q22=E2/(1-v12*v21);
Q66=G12;
% PZT layer刚度矩阵
Q0=[Q112 Q122 0
      Q122 Q222 0
      0      0    Q662];
% CF layer刚度矩阵
Q1=[Q11 Q12 0
      Q12 Q22 0
      0      0    Q66];
% 调整刚度矩阵
Qmid=Q1;
Qfirst=Q0;
Qthird=Q0;
%综合模型的各种系数
g_delta=(1+2*l_r);
g_c=8*(1-w_r)^3;
g_d=6*(w_r-1)*(3+4*l_r-2*w_r-4*l_r*w_r);
g_e=3*(-2-2*l_r+w_r+2*l_r*w_r)^2*log((2-w_r)/w_r);
G_l_ext=(g_d+g_e)/g_c;
G_F=g_delta/G_l_ext;
G_U=g_delta*G_F;
 
z0=-t_cf/2-t_pzt; 
z1=-t_cf/2; 
z2=t_cf/2; 
z3=t_cf/2+t_pzt;  

A=Qfirst*(z1-z0)+Qmid*(z2-z1)-Qthird*(z3-z2);
B=(Qfirst*(z1^2-z0^2)+Qmid*(z2^2-z1^2)+Qthird*(z3^2-z2^2))*(1/2);
D=(Qfirst*(z1^3-z0^3)+Qmid*(z2^3-z1^3)+Qthird*(z3^3-z2^3))*(1/3);
C=pinv([A B;B D]);

Fp=E3*(z1-z0)*Qfirst*d+E3*(z3-z2)*Qthird*d;  
Mp=E3*1/2*(z1^2-z0^2)*Qfirst*d-E3*1/2*(z3^2-z2^2)*Qthird*d;
Np=[Fp;Mp];
P=C(4,1)*Np(1)+C(4,2)*Np(2)+C(4,4)*Np(4)+C(4,5)*Np(5);
% ED_m=3/8*P*P/C(4,4)/(2*Rou_pzt*t_pzt+Rou_cf*t_cf)*G_U  % 能量密度
ED_m=3*P.^2*l_pzt*wnom*G_U./(8*C(4,4)*m_actuator)  % 能量密度
 
 td=P*l_pzt^2/2*g_delta;
 bf=3*P*wnom/(2*C(4,4)*l_pzt)*G_F;
 t_cf=t_cf*1e6
 td=td*1e6
 bf=bf*1e3
 %% 第一种尺寸%%%%%%%%%%%%%%%
% (1) output:  U=250; % V t_cf=40;  % um
% m_actuator =3.8518e-005;  % mg
% d_r =0.2175;
% ED_m =1.2078;   % 实际测试时的选择，使用250v驱动，能量密度最优。
% t_cf = 40;
% td = 399.8900;
% bf =232.6649;
% (2) output:  U=250; % V t_cf=60;  % um
% m_actuator =3.9178e-005;
% d_r =0.2289;
% ED_m =1.2017;
% t_cf =60.0000;
% td =361.4155;
% bf =260.5220;
% (3) output:  U=300; % V t_cf=40;  % um
% m_actuator =3.8518e-005;
% d_r =0.2175;
% ED_m =1.7392; % 能量密度最大
% t_cf =40;
% td =479.8680;
% bf =279.1979;
% (4) output:  U=300; % V t_cf=60;  % um
% m_actuator =3.9178e-005;
% d_r =0.2289;
% ED_m =1.7304;
% t_cf =60.0000;
% td =433.6986;
% bf = 312.6264;  

%% 第二种尺寸%%%%%%%%%%%%%%%
% (1) output:  U=250; % V  t_cf=40;  % um
% m_actuator =4.7484e-005;
% d_r =0.2175;
% ED_m =1.0873;
% t_cf =40;
% td =937.2421;
% bf =110.1760;
% (2) output:  U=250; % V  t_cf=60;  % um
% m_actuator =4.8534e-005;
% d_r =0.2289;
% ED_m =1.0766;
% t_cf = 60.0000;
% td =847.0677;
% bf =123.3674;
% (3) output:  U=300; % V t_cf=40;  % um
% m_actuator =4.7484e-005;  % mg
% d_r =0.2175;
% ED_m =1.5658;
% t_cf =40;
% td =1.1247e+003;
% bf =132.2112;
% (4) output: U=300; % V t_cf=60;  % um
% m_actuator =4.8534e-005;
% d_r =0.2289;
% ED_m =1.5503;
% t_cf = 60.0000;
% td =1.0165e+003;
% bf =148.0409;

