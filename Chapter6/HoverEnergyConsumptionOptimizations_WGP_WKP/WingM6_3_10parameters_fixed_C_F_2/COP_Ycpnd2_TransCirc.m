function Y_rcpnd_Trans=COP_Ycpnd2_TransCirc(alpha,R_wing,xr0,C_maxyaxis)  % 输入变量不含有C_aver;
% function Y_rcpnd_Trans=COP_Ycpnd2_TransCirc(alpha,C_maxyaxis)
% 修改时间——2014年12月21日,23:36
% 修改时间——2015年01月20日,11:16
% 修改时间——2015年05月20日,11:16
% 修改时间——2015年05月28日,17:27
% 修改时间——2015年05月29日,17:27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[3.004,0.8854,0.3289,0.255]; 
% % x=[3.004,0.8854,0.3289,0.35]; % %Ratio_leadmax=0.4644/1.3018=0.356737; 
% R_wing=x(1);
% % C_aver=x(2);
% xr0=x(3);
% C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 针对旋转轴气动力矩——求解净压心的无量纲位置Y_cpnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R=16.0148-0.88=15.1348;
% 当 r_nd∈(0.88/15.1348,13.4427/15.1348)
% clear all; clc;
R_wingeff=3.004;          %有效翅膀长度(mm) 
% C_avereff=0.8854;     % mm 
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  无量纲展向偏置距离   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅根-偏离-坐标系原点的距离
R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
x=linspace(R_proximal,R_distal,200);
x_mod_Root=0.636;                                            % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C_maxylb=0.464385778290230;
% C_maxy25=0.138924474377504;  % 针对程序wing_model_88_yaxis有: 第122行; C_maxy =0.1389; 
% C_maxyub=-0.186536829535222; 
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 未优化的果蝇翅膀形貌参数——翅根和翅尖连线的扭转轴定出的前缘
% 原始果蝇翅膀 R_wing=3.004;  C_aver=0.8854;   
% C_lead_ymax=0.4644;   % C_trail_ymin =-0.8374;  
C_max_LtoT= 1.3018;    % @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; 
% x_start=[3.004,0.8854,0.3289,0.356737]; % 初始值. %针对翅根和翅尖连线的扭转轴定出的前缘
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 未优化的果蝇翅膀形貌参数——扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时
% C_lead_ymax =0.3255;  % C_trail_ymin =-0.9764;  
% C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;
% x_start=[3.004,0.8854,0.3289,0.25];     % 初始值.  %针对扭转轴位于最大前缘点和最小后缘点坐标间距0.25倍时
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_pitchaxis=Ratio_leadmax-C_maxyaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 针对翅根和翅尖连线的扭转轴定出的前缘     
% C_maxy=C_lead_ymax-C_maxy*C_max_LtoT;  %转动轴偏离弦向中点的偏移量坐标,在扭转轴向上偏移 -C_maxy之后 % XXX
yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root;  % 有量纲
yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root;  % 有量纲
% C_rx=yr_lead-yr_trail;      % 正确——量纲化实际弦长分布
% 采用前缘拟合函数求解——实际平均弦长=面积/R_wingeff
% wing_aera=trapz(x,C_rx);             %输出: wing_aera =2.6597; % mm^2
% C_aver=wing_aera/R_wingeff    % 输出量纲化平均弦长: C_avereff=C_aver =0.8854; % mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) 无量纲前缘分布函数
r_nd=(x-xr)/R_wingeff;
% yr_leadnd0=yr_lead/C_avereff;
yr_leadnd0=yr_lead/C_max_LtoT;
% P_coeff_lead=polyfit(r_nd,yr_leadnd0,6);
% (b)  无量纲后缘分布函数
% yr_trailnd0=yr_trail/C_avereff;
yr_trailnd0=yr_trail/C_max_LtoT;
% P_coeff_trail=polyfit(r_nd,yr_trailnd0,6);
% (c)  无量纲弦长分布函数
Cr_nd=yr_leadnd0-yr_trailnd0;
P_coeff_Cr=polyfit(r_nd,Cr_nd,6);  % 多项式系数  % Cr_nd2=polyval(Coeff,r_nd1);
syms r_nd   % 无量纲弦长分布为6阶多项式——转换必须有这条指令
% yr_leadnd=vpa(poly2sym(P_coeff_lead,r_nd),6); 
% yr_trailnd=vpa(poly2sym(P_coeff_trail,r_nd),6);
C_nd=vpa(poly2sym(P_coeff_Cr,r_nd),6);  
% % 方案(1)——由前后缘函数的无量纲化yr_leadnd——yr_trailnd——求得
% 下面是翅前缘函数——针对扭转轴的位置不同需要进行分段函数处理吗?
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% 无量纲弦长分布为6阶多项式
% C_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 弦向压心点的分布—决定扭转轴力臂的方向改变
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();        %调用函数kenimatics_wing_and_AoA
% % size(wing_kenimatics)                  %  (1000*12)
% alpha=wing_kenimatics(:,4);
% alpha=pi/2.4;                                   % 单位是rad
% alpha=59.265*pi/180;                      % alpha=62.775*pi/180;
% alpha=(90-45)*pi/180;
d_cprnd=0.82*abs(alpha)/pi+0.05;     % 旋转轴气动力矩的弦向压心位置; %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% yr_cpnd=yr_nd+yr_leadnd-C_nd*d_cprnd; 
%%第一方案——不够合理哦
%  if d_cprnd>0.25          % 在轴之后——问题在于d_cprnd是向量，而不是一个数，例如假设(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=C_nd*d_cprnd-yr_leadnd;     %  z3=-(d_cp-C_025);    % 正号
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  elseif d_cprnd<=0.25; % 在轴之前——问题在于d_cprnd是向量，而不是一个数，例如假设(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=yr_leadnd-C_nd*d_cprnd;    %  z3=C_025-d_cp;          % 正号
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  end
%%第二方案——可以考虑选取相应的展向气动压心位置所对应的片条进行分析——弦向压心 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 多项式值比值——不可行
% Ratio_axis=vpa(yr_leadnd./C_nd, 5);  % vpa(f,d); d是指有效数字的位数
% f_1=vpa(yr_leadnd,5)
% f_11=factor(f_1)
% f_2=vpa(C_nd,5)
% f_22=factor(f_2)
% Ratio_axis=vpa(f_11./f_22, 5)
% simple(Ratio_axis)
% Ratio_axis =
% -(1.0*(- 68.464*r_nd^6 + 208.3*r_nd^5 - 245.23*r_nd^4 + 137.47*r_nd^3 - 36.859*r_nd^2 + 4.7923*r_nd + 0.0008352))...
%    /(40.826*r_nd^6 - 87.204*r_nd^5 + 59.427*r_nd^4 - 11.864*r_nd^3 + 3.7541*r_nd^2 - 4.9379*r_nd + 0.000057821)  
%%% ——不可行
%  if d_cprnd>yr_leadnd/C_nd            % 在轴之后—— 问题在于d_cprnd是向量，而不是一个数，例如0.25
%      yr_cpnd=C_nd*d_cprnd-yr_leadnd;     %  z3=-(d_cp-C_025);     % 正号   
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  elseif d_cprnd<=yr_leadnd/C_nd    % 在轴之前
%      yr_cpnd=yr_leadnd-C_nd*d_cprnd;    %  z3=C_025-d_cp;          % 正号
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 选取相应的展向气动压心位置所对应的片条进行分析——弦向压心 
xr_nd_vari=xr0/R_wing;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx1=(r_nd+xr_nd_vari)^2*C_nd;    % 无量纲气动力F_nd的原始被积函数
% fx2=vpa(fx1,5)
fx3=expand(fx1);
F_ndTrans=double(int(fx3,r_nd,0,1));   % Result: F_ndTrans =0.46391;
%%%%%%%%%%%%%%%%%%%%%%%%
fx4=(r_nd+xr_nd_vari)^3*C_nd;  % 这里C_nd由前后缘多项式函数之差给出
fr_nd5=expand(fx4);
I1=double(int(fr_nd5,0,1));       % 无量纲,量纲化单位是mm^5;   I1 =0.3659;
r_xcopnd_tr=I1/F_ndTrans;       % r_xcopnd_tr= 0.7887;       %    *R_wingeff
%%%%%%%%%%%%%%%%%%%%%%%%%
% % R_tr=(R_wingeff+xr)*r_xcopnd_tr;
% % r_nd_tr=(R_tr-xr)/R_wingeff;                % r_nd_tr =0.7656;
% r_nd_tr =0.788691874094779;  % ≈ 0.7887; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ratio_axis=vpa(yr_leadnd./C_nd, 5);        % vpa(f,d); d是指有效数字的位数
% Ratio=inline(vectorize(Ratio_axis),'r_nd');
% Ratio_axistr=Ratio(r_nd_tr)                         % Ratio_axistr =0.2901; ————————一个非常重要的参数
% % % 下面的解释程序详细些
% % yr_leadnd1=vpa(yr_leadnd, 5);          % vpa(f,d); d是指有效数字的位数
% % yr_leadnd=inline(vectorize(yr_leadnd1),'r_nd');
% % yr_leadnd_tr=yr_leadnd(r_nd_tr)       % yr_leadnd_tr =0.3398; ————————一个非常重要的参数
% % C_nd1=vpa(C_nd, 5);                         % vpa(f,d); d是指有效数字的位数
% % C_nd=inline(vectorize(C_nd1),'r_nd');
% % C_nd_tr=C_nd(r_nd_tr)                        % C_nd_tr =1.1714; ————————一个非常重要的参数
% % Ratio_axistr=yr_leadnd_tr/C_nd_tr         % Ratio_axistr =0.2901; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 确定平动气动力展向压心处对应片条的无量纲前缘和无量纲弦长 @ r_nd_tr =0.7656;
% yr_leadnd_tr1=inline(vectorize(yr_leadnd),'r_nd');
% yr_leadnd_tr=yr_leadnd_tr1(r_xcopnd_tr);   % yr_leadnd_tr =0.5027;   % yr_leadnd_tr =0.3359;
C_nd_tr1=inline(vectorize(C_nd),'r_nd');
C_nd_tr=C_nd_tr1(r_xcopnd_tr);      % C_nd_tr =1.2034 = C_copnd_tr =1.0655/0.8854=1.2034;  
% alpha_crit_abs=(Ratio_axis-0.05)*pi/0.82*180/pi  % 当d_cprnd=Ratio_axis时, 确定临界绝对值角度alpha_crit_abs=81.2374;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if d_cprnd>Ratio_axis  %    % 在轴之后——问题在于d_cprnd是向量，而不是一个数，例如0.25
%      yr_cpnd=-(C_nd_tr*d_cprnd-yr_leadnd_tr); %z3=-(d_cp-C_025);   % 负号  %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% elseif d_cprnd<=Ratio_axis  % 在轴之前
%      yr_cpnd=yr_leadnd_tr-C_nd_tr*d_cprnd;    %z3=C_025-d_cp;       % 正号  %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if d_cprnd>C_maxyaxis  %    % 在轴之后——问题在于d_cprnd是向量，而不是一个数，例如0.25
%      yr_cpnd=-(d_cprnd-C_maxyaxis); 
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% elseif d_cprnd<=C_maxyaxis  % 在轴之前
%      yr_cpnd=C_maxyaxis-d_cprnd; 
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% end
 yr_cpnd=C_nd_tr*(C_maxyaxis-d_cprnd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yr_cpnd=yr_leadnd_tr-C_nd_tr*d_cprnd-delta_pitchaxis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx1=(r_nd+xr_nd)^2*C_nd;                  % 无量纲气动力F_ndTrans的原始被积函数
fx3=expand(fx1);
F_ndTrans=double(int(fx3,r_nd,0,1));           % Result: F_ndTrans=0.46392
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 无量纲气动力分量(nondimention_aerodynamic_component)的求解——F_ndTrans
% 注意——该段程序切记不得修改，前提只要保证输入正确的无量纲弦长分布即可。
%以下的公式应使用合理的无量纲的弦长分布公式C_nd
% R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1));   % 二阶面积矩的回转半径的平方
% R1nd1=double(int(r_nd*C_nd,r_nd,0,1));        % 一阶面积矩的回转半径
% % S_nd=double(int(C_nd,r_nd,0,1));               % 无量纲翅面积
% F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;      % 使用这句计算结果也正确; 输出:F_nd2 =0.5024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_Trans=int(yr_cpnd*fx3,r_nd,0,1)/F_ndTrans;
% disp(['净压心的无量纲位置Y_cpnd_Trans(alpha)=' num2str(Y_rcpnd_Trans)  ' 量纲单位可能是mm'])
% Y_rcpnd_Trans=abs(double(int(yr_cpnd*fx3,r_nd,0,1))/F_ndTrans);
Y_rcpnd_Trans=double(int(yr_cpnd*fx3,r_nd,0,1))/F_ndTrans;


