function I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis)   %调用函数inertia_moment;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moment_inertia_scaling_law
R_wingeff=3.004;     % mm    % R_ratio =33.2889;
C_avereff=0.8854;    % mm    % C_ratio =38.0619;
R_ratio=R_wing/R_wingeff;
C_ratio=C_aver/C_avereff; 
%% 翅膀形貌学参数变化是的质心位置坐标
% 初始模型翅膀的质心坐标——下面的数据是针对扭转轴为翅根和翅尖的连线时翅膀的质心位置坐标
% 翅膀质心位置坐标: mm
% XC=1.920243385;      % X=1.920243385;          % 展    向
% YC=0.000000068;      % Y=0.000000068;          % 厚度方向——不作考虑
% ZC=-0.149785466;    %  Z=-0.149785466         % 弦 向  % 到扭转轴的相对距离: 用于实际计算
xr1=0.3289;  % 针对右侧翅膀——Z_0轴偏离翅根距离(在翅根左侧)
XC=1.920243385-xr1; % 1.5914
% C_maxy=0.138924474377504; % 针对程序wing_model_88_yaxis有: 第122行; C_maxy =0.1389; 
% YC=-0.149785466-C_maxy;    % YC=-0.288709941;
YC=-0.149785466;
% 翅膀形貌学参数变化时——质心坐标
XC_tran=R_ratio*XC+xr0;           % XC_tran =53.3030@R_wing=100; & xr1=0.3289; % mm;
YC_tran=C_ratio*YC;                    % YC_tran =-5.7011@C_aver=33.7; % mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀形貌学参数变化时的质心处的——惯性矩
% 初始模型翅膀的——质心惯性矩
% Mass_Properties
% M=0.000002251;         % 质量 % unit: g
% 惯性矩（质心）         % g.mm^2
% Ixc=0.000000217;        % g.mm^2
% Iyc=0.000001362;    % g.mm^2
% Izc=0.000001145;        % g.mm^2 
% % 惯性积（质心）     % g.mm^2
% Iyzc=-0.000000000;     
% Ixzc=0.000000053; 
% Ixyc=-0.000000000; 
Ixc_small=0.000000215;     % g.mm^2
Izc_small=0.000001129;      % g.mm^2
% 翅膀形貌学参数变化时的——质心惯性矩
Ixc_big=R_ratio*C_ratio^3*Ixc_small;    % Ixc_big =0.3946@R_wing=100; % mm;   % g.mm^2
Izc_big=R_ratio^3*C_ratio*Izc_small;    % Izc_big =1.5852@C_aver=33.7; % mm;   % g.mm^2
%% 翅膀形貌学参数变化变化前后——质量
M_small=0.000002237;  % g
M_big=R_ratio*C_ratio*M_small;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀形貌学参数变化时的相对于拍打轴和扭转轴的——惯性矩——平移轴定理
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xr_initial=0.3289;                        % x-root offset  \mm                    
% x_com_zaxis=XC-xr_initial+xr0;   % ————————————————————被更新
x_com_zaxis=XC_tran;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZC_tran=YC_tran;
% x_0∈[0,0.5];  % 扭转轴距前缘的无量纲距离
% C_maxylb=0.464385778290230;    % C_maxylb =0.4644;
% C_maxy25=0.138924474377504;   % C_maxy25 =0.1389; 
% C_maxyub=-0.186536829535222; % C_maxyub =-0.1865
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %针对翅根和翅尖连线的扭转轴定出的前缘
Ratio_leadmax=0.356737;  % 无量纲
% C_lead_ymax=max(f_x_lead);   % 输出: C_lead_ymax=0.4644; k_leadmax=644;
% C_lead_ymax=0.464385778290230;               % mm——翅膀形貌学参数变化时最大前缘点z坐标
wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis);
C_max_LtoT=wing_para(1,17);
z_maxlead=Ratio_leadmax*C_max_LtoT;            % mm——翅膀形貌学参数变化时最大前缘点z坐标
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_maxlead_xaxis=C_maxyaxis*C_max_LtoT;        % mm; 扭转轴距前缘的量纲距离
delta=z_maxlead-z_maxlead_xaxis;  % 翅根和翅尖连线轴到新扭转轴之间的量纲距离
z_com_zaxis=abs(abs(ZC_tran)+delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 惯性矩 (WCS)   g.mm^2
% Ix=Ixc+M*ZC^2         % 由平移轴定理获得
% Iy=Iyc+M*XC^2         % 由平移轴定理获得
% Iz=Izc+M*XC^2         % 由平移轴定理获得
Ix_inertia=Ixc_big+M_big*z_com_zaxis^2;     % 由平移轴定理获得       % g.mm^2
Iz_inertia=Izc_big+M_big*x_com_zaxis^2;     % 由平移轴定理获得       % g.mm^2
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];  % XXX
I_moment=[x_com_zaxis,z_com_zaxis,M_big,Ix_inertia,Iz_inertia];  % 这个输出更为合理些
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



