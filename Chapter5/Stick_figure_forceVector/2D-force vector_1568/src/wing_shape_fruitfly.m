function wing_para_output=wing_shape_fruitfly()
% 创建时间——2014年6月19日,0:16:04
% 修改时间——2014年6月14日,23:35:03
% clear all;clc;
%% 翅形貌参数化——由Hedrick程序计算而得
%拟合得到 R_wingeff=3.3328-0.3289=3.004;  Hedrick程序计算而得: R_wingeff=3.007; Science数据是: R_wingeff=2.99; 
R_wingeff=3.004;    %有效翅膀长度(mm)  
% Hedrick程序计算而得: C_avereff=0.884;     % Science数据是C_avereff=0.9468mm;  
C_avereff=0.8854;  %见后文—C_avereff=C_aver =0.8854;—由前后缘实际拟合曲线函数相差求积分均值而得
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
xr_nd1=xr/R_wingeff1;      % x-root offset  无量纲展向偏置距离
F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %无量纲气动力分量F_nd1, 这里输出: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 第二种方式求解无量纲气动力分量——wing_shape——翅前后缘拟合函数
% 求取前后缘的交点坐标 
% yr_lead= -2.607e-015*x.^11.78+0.8139-0.806122+0.73;    % 前缘的拟合函数——乘方函数
% yr_lead= -2.607e-015*x.^11.78+0.737778;                          % 由上面一条指令获得
% yr_trail=-0.0017*x.^3+0.1073*x.^2-1.3182*x+0.7783;        % 后缘的拟合函数——三阶多项式
% p(x)=-0.0017*x.^3+0.1073*x.^2-1.3182*x+0.7783+2.607e-015*x.^11.78-0.737778;  %采用Maple求根
% fsolve(p2(x) = 0, x, x = 10 .. 18);     %采用Maple求根   % x=16.13067567;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅根-偏离-坐标系原点的距离
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
x=linspace(R_proximal,R_distal,200);
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413;  
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282;
% figure(1)  % 翅形貌——采用前缘拟合函数
% plot(x,yr_lead,'r-',x,yr_trail,'b-')
% xlabel('展向r (mm)')
% ylabel('前缘和后缘弦长 (mm)')
% legend('前缘拟合曲线','后缘拟合曲线')
% title('R_wingeff=3.004mm: 翅几何形貌的前缘和后缘拟合曲线')
C_rx=yr_lead-yr_trail;              % 正确——量纲化实际弦长分布
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 采用前缘拟合函数求解——实际平均弦长=面积/R_wingeff
wing_aera=trapz(x,C_rx);             %输出: wing_aera =2.6597; % mm^2
C_aver=wing_aera/R_wingeff;   % 输出量纲化平均弦长: C_avereff=C_aver =0.8854;
% % 采用第二种积分方法求解，貌似不对————XXXX
% yr_lead1=inline('-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413','x'); 
% yr_trail1=inline('-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282','x'); 
% wing_aera1=abs(quadl(yr_lead1,R_proximal,R_distal))+abs(quadl(yr_trail1,R_proximal,R_distal)); %输出: wing_aera1 =3.1119; % mm^2
% C_aver1=wing_aera1/R_wingeff;    % 输出: C_aver1 =1.0359;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 求解——无量纲弦长分布多项式函数
r_nd1=linspace(0,1,200);    % r=r_nd*R_wingeff;
Cr_nd1=C_rx/C_aver;          % 无量纲弦长分布  
% cftool
% 针对r_nd1和Cr_nd1采用函数cftool调用拟合工具箱进行无量纲弦长分布函数拟合，获得如下K阶多项式
% 输出结果
% Linear model Poly6:
%      f(x) = p1*x^6 + p2*x^5 + p3*x^4 + p4*x^3 + p5*x^2 + 
%                     p6*x + p7
% Coefficients (with 95% confidence bounds):
%        p1 =      -40.83  (-40.83, -40.83)
%        p2 =       87.21  (87.21, 87.21)
%        p3 =      -59.43  (-59.43, -59.43)
%        p4 =       11.86  (11.86, 11.86)
%        p5 =      -3.754  (-3.754, -3.754)
%        p6 =       4.938  (4.938, 4.938)
%        p7 = -5.782e-005  (-5.782e-005, -5.782e-005)
% Goodness of fit:
%   SSE: 4.697e-026
%   R-square: 1
%   Adjusted R-square: 1
%   RMSE: 1.56e-014
% 获得无量纲弦长分布——6阶多项式
% Cr_nd2=-40.83*r_nd1.^6 + 87.21*r_nd1.^5-59.43*r_nd1.^4+11.86*r_nd1.^3-3.754*r_nd1.^2+4.938*r_nd1-5.782e-005;
% figure(2)
% subplot(211)
% plot(x,C_rx,'r-')
% xlabel('展向r (mm)')
% ylabel('量纲化实际弦长(mm) ')
% legend('量纲化实际弦长')
% title('R_wingeff=3.004mm: 翅几何形貌的量纲化实际弦长分布')
% subplot(212)
% plot(r_nd,Cr_nd2,'b-')
% xlabel('展向\itr_{\rmnd} (I)')
% ylabel('无量纲弦长(I)')
% legend('无量纲弦长')
% title('R_wingeff_nd=1: 翅几何形貌的无量纲弦长分布')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 无量纲气动力分量(nondimention_aerodynamic_component)的求解——F_nd
Coeff=polyfit(r_nd1,Cr_nd1,6);  % 多项式系数  % Cr_nd2=polyval(Coeff,r_nd1);
syms r_nd
Cr_nd2=vpa(poly2sym(Coeff,r_nd),6);   % 无量纲弦长分布为6阶多项式——这里Cr_nd3=Cr_nd2;只是自变量由r_nd1变成了r_nd
% Cr_nd2 =-40.827*r_nd^6+87.2061*r_nd^5-59.4281*r_nd^4+11.8648*r_nd^3-3.75417*r_nd^2+4.938*r_nd-0.0000578229;
C_nd=Cr_nd2;
% 注意——该段程序切记不得修改，前提只要保证输入正确的无量纲弦长分布即可。
%以下的公式应使用合理的无量纲的弦长分布公式C_nd
R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1)); %二阶面积矩的回转半径的平方
R1nd1=double(int(r_nd*C_nd,r_nd,0,1));     %一阶面积矩的回转半径
S_nd=double(int(C_nd,r_nd,0,1));                %无量纲翅面积
disp(['二阶面积矩的回转半径的平方: r2_2nd=' num2str(R2nd2)  ' 量纲单位是mm^4'])
disp(['二阶面积矩的回转半径: r_2nd=' num2str(sqrt(R2nd2))  ' 量纲单位是mm^3'])
disp(['一阶面积矩的回转半径: r_1nd=' num2str(R1nd1)  ' 量纲单位是mm^3'])
disp(['无量纲翅面积Swing_nd=' num2str(S_nd)  ' 量纲单位是mm^2'])
% C_nd=vpa(C_nd,5)  %函数vpa将符号表达式中的数值(常转换为两个整数的比值,即分数)，转换为十进制小数表示。
% 假设xr_nd=0，则自变量r_nd的取值范围为:r_nd∈[0,1], 在取1时r_nd=R=11.95 /mm
fx2=(r_nd+xr_nd)^2*C_nd;    % 无量纲气动力F_nd的原始被积函数
% fx3=vpa(fx2,5)
fx4=expand(fx2);
F_nd=double(int(fx4,r_nd,0,1));                    % Result: F_nd =0.46392;
disp(['无量纲气动力F_nd=' num2str(F_nd)  ' 量纲单位是mm^4'])
F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;    %使用这句计算结果也正确; 输出:F_nd2 =0.46392;
disp(['无量纲气动力F_nd=' num2str(F_nd2)  ' 量纲单位是mm^4'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输出:
% 二阶面积矩的回转半径的平方: r2_2nd=0.3358 量纲单位是mm^4
% 二阶面积矩的回转半径: r_2nd=0.57949 量纲单位是mm^3  ——————————需要该数据
% 一阶面积矩的回转半径: r_1nd=0.53033 量纲单位是mm^3
% 无量纲翅面积Swing_nd=1 量纲单位是mm^2
% 无量纲气动力F_nd=0.46392 量纲单位是mm^4    ——————————需要该数据
% 无量纲气动力F_nd=0.46392 量纲单位是mm^4    ——————————需要该数据
%  对比下面由 (1) 第一种方式求解获得的无量纲气动力分量
% F_nd1=r2_nd^2+2*xr_nd1*r1_nd+xr_nd1^2;   %无量纲气动力分量F_nd1, 这里输出: F_nd1 =0.4635
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 针对气动阻尼力矩——求解有效力臂的无量纲位置Y_nd
% % 翅根-偏离-坐标系原点的距离
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
% x=linspace(R_proximal,R_distal,200);
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413;  
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282;
syms r_nd
yr=0;                              % 扭转轴通过翅根与翅尖的连线—\mm
yr_nd=yr/C_avereff;       % y-root offset  无量纲弦向偏置距离 yr_nd = 0.0214;
yr_lead1=-0.08249*r_nd^6+0.9167*r_nd^5-4.04*r_nd^4+8.872*r_nd^3-10.06*r_nd^2+5.674*r_nd-0.413;
yr_leadnd=yr_lead1/C_avereff;
% % 方案(1)——由前后缘函数的无量纲化yr_leadnd——yr_trailnd——求得
% yr_trail1=-0.0333*r_nd^6+0.504*r_nd^5-2.795*r_nd^4+7.258*r_nd^3-8.769*r_nd^2+3.739*r_nd+0.1282;
% yr_trailnd=yr_trail1/C_avereff;
% yr_nd1=expand((yr_leadnd^4+yr_trailnd^4)/4);     %被积函数2
% Y_rnd=double(int(yr_nd1,r_nd,0,1));                         % 输出: Y_nd=0.15341;  %  mm
% 方案(2)——由前缘函数的无量纲化yr_leadnd和无量纲化弦长分布C_nd求得
y0=yr_nd+yr_leadnd-C_nd;
y1=yr_nd+yr_leadnd;
yr_nd2=(y1^4+y0^4)/4;   % 注意这里的扭转轴直接通过翅尖和翅根的连线，前缘函数的翅膀展向区间唯一
Y_rnd=double(int(yr_nd2,r_nd,0,1));                           % 输出: Y_rnd2=0.1402;
disp(['有效力臂的无量纲位置Y_nd=' num2str(Y_rnd)  ' 量纲单位是mm'])  % 结果可由Maple求得
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 虚质量效应的系数求解和虚质量转动惯量——虚质量气动力和力矩
% yr_lead1=-0.08249*r_nd^6+0.9167*r_nd^5-4.04*r_nd^4+8.872*r_nd^3-10.06*r_nd^2+5.674*r_nd-0.413;
% yr_leadnd=yr_lead1/C_avereff;
% yr=0;                              % 扭转轴通过翅根与翅尖的连线—\mm
% yr_nd=yr/C_avereff;       % y-root offset  无量纲弦向偏置距离 yr_nd = 0.0214;
yr_hnd=C_nd/2-yr_leadnd-yr_nd;    %转动轴偏离弦向中点的偏移量坐标
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
% 虚质量气动力矩参数——虚质量转动惯量
fx5=(r_nd+xr_nd)*C_nd^2*yr_hnd;
fr_nd1=expand(fx5);
fx6=C_nd^2*(yr_hnd^2+C_nd^2/32);
fr_nd2=expand(fx6);
I_xzamnd=int(fr_nd1,0,1);
I_xxamnd=int(fr_nd2,0,1);
I1=double(I_xzamnd);                                          % 无量纲,量纲化单位是mm^5;          % I1 =0.2666;  
I2=double(I_xxamnd);                                          % 无量纲,量纲化单位是mm^5;          % I2 =0.1637;  
Rou=1.225*10^(-3);           %单位是Kg/m^3=10^6/(10^3)^3=10^(-3)mg/mm^3
I_xzam=pi*Rou*C_avereff^3*R_wingeff^2*I1/4;  % 单位是 mg.mm^2;   % I_xzam =2.3951;
I_xxam=pi*Rou*C_avereff^4*R_wingeff*I2/4;      % 单位是 mg.mm^2;   % I_xxam = 0.3279;
%%%%%%%%%%%%%%%%%%%%%%%%
% 虚质量气动力参数
fx7=(r_nd+xr_nd)*C_nd^2;                                  % 这里C_nd由前后缘多项式函数之差给出
fr_nd3=expand(fx7);
fx8=C_nd^2*yr_hnd;
fr_nd4=expand(fx8);
I3=double(int(fr_nd3,0,1));                                  % 无量纲,量纲化单位是mm^4; 在转动气动力中也用到了I3
I4=double(int(fr_nd4,0,1));                                  % 无量纲,量纲化单位是mm^4; 
I3z=pi*Rou*C_avereff^2*R_wingeff^2*I3/4;       % 单位是 mg.mm
I4z=pi*Rou*C_avereff^3*R_wingeff*I4/4;            % 单位是 mg.mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para_output=zeros(1,9);
wing_para_output=[R_wingeff,C_avereff,F_nd,Y_rnd,I_xzam,I_xxam,I3,I3z,I4z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
