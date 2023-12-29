function obj_function=morphology_para_opt_obj_function1(x)
% 功率调用函数和气动力调用函数一起调用同一个函数Aero_F_M_fruitfly——程序较简洁
% morphology_para_opt_obj_function
% morphology parametric optimization objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 变量
% % f∈[0,inf];  phi_m∈[0,pi/2];  K∈[0,1]; % eta_m∈[0,pi]; C_eta∈[0,inf];  Phi_eta∈[-pi,pi]; eta_0∈[eta_m-pi,pi-eta_m];
% eta_m=x(8);
% eta_0min=eta_m-pi;
% eta_0max=pi-eta_m;
% LB = [2, 0.5, 0, 0,0,0,0,0,0,-pi,eta_0min];            % Lower bound       % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
% UB = [4, 2, 2, 0.5,inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound      % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第一部分——气动功率
% 原始果蝇翅膀 R_wing=3.004;  C_aver=0.8854;   
% C_lead_ymax =0.4644;  C_trail_ymin =-0.8374;  C_max_LtoT = 1.3018; @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT;
% %Ratio_leadmax=0.4644/1.3018=0.356737; %针对翅根和翅尖连线的扭转轴定出的前缘/最大前缘点和最新后缘点坐标间距=0.3567
% x=[3.004,0.8854,0.3289,0.356737]; % 初始值. % 未优化的果蝇翅膀形貌参数—翅根和翅尖连线的扭转轴定出的前缘/最大前缘点和最新后缘点坐标间距=0.3567
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%针对扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时
% x=[3.004,0.8854,0.3289,0.25]; % 初始值. % 未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时
% x=[3.004,0.8854,0.3289,0.5]; % 初始值. % 未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.5倍时
% x=[100,33.7,0.3289,0.356737];  %验证被放大后的翅膀的转动惯量——非常棒
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 变量单位: mm & 扭转轴距前缘的量纲距离[0,0.6509]——注意这里可能不用——XXX
% x=[3.004,0.8854,0.3289,0.1389];     % 未优化的果蝇翅膀形貌参数
% C_maxy25_nd=(0.464385778290230-0.138924474377504)/1.3018=0.3255/1.3018;  %C_maxy25_nd=0.25;
% x=[3.004,0.8854,0.3289,0.3255]; %翅膀被放大之前,扭转轴在翅膀最大前缘点到最小后缘点的距离1/4位置@0.25*C_max_LtoT, 见上一句指令的解释
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %针对翅根和翅尖连线的扭转轴定出的前缘
% x=[3.004,0.8854,0.3289,0.4644];     %翅膀被放大之前,扭转轴在翅根翅尖连线上@0.356737*C_max_LtoT, 见上一句指令的解释
% x=[3.004,0.8854,0.3289,0];          % 翅膀被放大之前,扭转轴在翅膀前缘最大值点处@0*C_max_LtoT
% x=[3.004,0.8854,0.3289,0.6509]; % 翅膀被放大之前,扭转轴在翅膀最大前缘点到最小后缘点的距离一半位置@0.5*C_max_LtoT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2014-Science-实验测试而得运动学数据；翅膀形貌学参数为：
% % x=[3.004,0.8854,0.3289,0.25]; %初始值. %未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时——结果
% % P_asterisk=31.8058; L=1.3940; % 未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时——结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%@ options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400);
% x =[3.6470,1.9176,0.6487,0.4258]; 
%@options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.3086,0.9363,0.4644]; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.2895,0.9238,0.4644]; %ConstraintFunction=@aspectratio_constraint; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fmincon,fminconOptions});
% x =[2.0000,1.1161,1.0711,-0.1865];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % 不合理
% x=[4.0000,1.4283,0.9992,0.5000];  % 不合理
% x =[4,2,0.8188,0.5];                       % L = 1.0009;
% x=[2.2899,1.0561,0.3788,0.0001];  % 20150204-Hybrid GA-Fmincon计算结果
% x =[3.1133,0.5938,0.3297,0.2450]; % 20150204-Hybrid GA-Fminsearch计算结果
% x=[2.5088,1.1424,0.7124,0.0458];   % 20150209-Hybrid GA-Fminsearch计算结果
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[108.5082,1.7828,3.0150,0.7854,1.0431,0.6109,0.0749]; %——20140301-hybrid_GA_fminsearch_WingM4
% output: P_asterisk =27.4317;   L =1.0000;   delta =2.7839e-005;
% x_start= [120.4862, 1.5400,2.8107,-0.7854,0.9935,0.0114,0.0381]; %——20140304-hybrid_GA_fminsearch_WingM4
% output: P_asterisk =19.4628;   L =1.0000; delta =-1.2695e-009;
% x =[3.1935,1.6054,2.1531,0.5000,1.2853,0.8770,0.6763,0.1370,1.5444,-2.0757,0.2672]; %20150315-hybrid_GA_fminsearch_WingM4_3-11变量-结果
% output: P_asterisk =13.5366;    L =0.9762;  delta =-0.0238;
% x=[3.1940,1.7834,1.9737,0.4248,2.1375,0.3594,0.7267,-0.4506,0.8049,-2.9541,0.1517];
% output: P_asterisk =3.7620;     L =0.9998;       delta =-1.9874e-004;
% x=[108.5082,1.7828,3.0150,0.7854,1.0431,0.6109,0.0749 ,2.5088,1.1424,0.7124,0.0458]; %——20140301-hybrid_GA_fminsearch_WingM4
% output: P_asterisk =6.0113;   L =0.6944   delta =-0.3056; % 采用这组数据组合得到的结果是升重比小于1
% x=[120.4862,1.5400,2.8107,-0.7854,0.9935,0.0114,0.0381,2.5088,1.1424,0.7124,0.0458];%——20140304-hybrid_GA_fminsearch_WingM4
% output: P_asterisk =3.9471;   L =0.5966 ;  delta =-0.4034; % 采用这组数据组合得到的结果是升重比小于1
% x=[2.5265,0.5467,1.2822,0.6052,1.2692,-1.0154,0.0470,2.0970,1.3003,1.9931,0.0008];
% output: P_asterisk =4.8268;      L =1.0001;      delta =6.9477e-005;
% x =[1.0944,0.7427,0.8889,0.5062,1.9926,3.5444,1.6864,2.3574,0.4257];
% output: P_asterisk =3.5631;  L =1.0007;  delta =7.1246e-004;
% x=[120.4862,1.5400,2.8107,0.9935,0.0114,2.5088,1.1424,0.7124,0.0458];
% x=[126.0625,1.3884,2.3453,0.8991,0.2622,2.9992,1.1188,0.8973,0.2395];
% output: P_asterisk =2.0616;    L =0.9999;   delta =-1.1719e-004;
% x=[129.8096,1.5350,2.3679,0.8461,-0.0488,3.1267,1.1578,0.7896,0.1624];
% output: P_asterisk =3.9132;     L =0.9999;   delta =-7.2005e-005;
% x=[188.6980,1.2630,2.6270,0.4302,1.0689,0.0852,-0.0563,2.5604,0.9805,0.8083,0.0946]; % 11个变量
% output: P_asterisk =5.7957;   L =1.0000;  delta =4.4414e-006;
% x =[189.6461,1.4858,2.7456,0.5867,1.3172,-0.0779,0.0181,2.5223,1.0679,0.7204,0.2377];
% output: P_asterisk = 9.3661 L =0.9999;  delta =-1.1967e-004;
% x=[189.0216,1.4714,2.9520,0.0080,1.3754,-0.3428,0.0222,2.4904,1.1268,0.7390,0.2835];%20150428_WingM4_4_3_11variable_group-结果
% output: P_asterisk =10.1381;     L =1.0000;         delta =2.6063e-005;
% x =[74.7130,1.4522,2.5332,0.1185,0.7491,0.0379,0.0081,3.5848,1.4845,1.9894,0.3193];%20150519-WingM4_5_11variable_group-结果
% output: P_asterisk =88.6723;     L =2.2113;  delta =1.2113;
% x=[35.4018,1.5769,-0.5250,-0.3650,0.8572,-3.1414,0.0167,3.9583,1.8190,1.9995,-0.0486]; % 11个变量
% x=[35.4018,1.5769,-0.5250,0.8572,-3.1414,0.0167,3.9583,1.8190,1.9995,-0.0486]; % 10个变量
% x =[48.9196,1.4864,2.5901,0.9381,-0.1470,-0.0008,2.7621,2.0000,1.9383,-0.0080];% 10个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[45.4147,1.5696,1.8489,1.0385,-0.9685,-0.0463,3.5193,1.8689,1.3038,-0.0008];% 10个变量
% P_asterisk =6.9607;     L =1.0002;   delta =1.6715e-004;   penaltyfun1 =0.3343; 
% penaltyfun2 =4.5714;   AR =2.5807;   Re =126.7144;     penaltyfun5 =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[45.4271,1.5679,1.8515,1.0364,-0.9666,-0.0470,3.5258,1.8657,1.3081,0.0122];% 10个变量
% P_asterisk =7.0267;   L =1.0000; delta =2.6792e-005;   penaltyfun1 =0.0536;
% penaltyfun2 =0;    AR =2.5909;   Re =126.6284;    penaltyfun5 =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[38.7818,1.5259,2.6496,0.9744,-0.0440,0.0015,3.6837,1.7620,2.0000,0.0055];
% psi_0_deg=0.0015*180/pi; 
% P_asterisk =9.3866;   L =1.0000;    delta =-3.3726e-005;   penaltyfun1 =0.0675;   
% penaltyfun2 =0;      AR =3.2257;   Re =160.1730(103.8108);        penaltyfun5 =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[127.25,1.5259,2.6496,0.9744,-0.0440,0.0015,3.004,0.8854,0.3289,0.356737];
% Rario_F_vertical_horiz_aver =1.3266; F_vertical_to_horiz_rms =3.2174; 
% P_aero_aver_to_inertia_aver =11.7265;  P_asterisk =23.3323; L =1.0051;
% x =[188.7,1.443,2.2757,1.1398,-0.1799,-0.0203,3.004,0.8854,0.3289,0.356737];
% Rario_F_vertical_horiz_aver =1.1781;  F_vertical_to_horiz_rms =3.0437;
% P_aero_aver_to_inertia_aver =10.2407;  P_asterisk =59.8796;  L =1.8416;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % f=185.5;        % P_asterisk =15.5043; L =1.0053;
% f=188.7;            % P_asterisk =16.0184;  L =1.0223;
% phi_m=pi/3;  % P_asterisk =15.5043; L =1.0053;
% K=0.0001;
% eta_m=1.0157;  % eta_m=45*pi/180;  
% C_eta=2.375; 
% Phi_eta=-pi/2;  
% x=[f,phi_m,K,eta_m,C_eta,Phi_eta,3.004,0.8854,0.3289,0.356737];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[-29.6664,0.0000,1.0278,1.5708,25.9878,-1.2863,3.6594,1.7727,1.9991,0.0002]; % AR =3.1921; % fval =2.0574e+003 +1.7778e-007i;
% x =[-29.6676,0.1662,0.3033,1.5707,25.9866,0.0444,3.6592, 1.7728,1.9995,0.0002]; % AR =3.1920; % fval =2.0014e+003;
% x =[49.225,1.0472,0.99,0.7854,25.9878,-1.2863,3.6594,1.7727,1.9991,0.0002]; % P_asterisk =11.0877; L =1.0000; Re =139.7495;
% x=[-36.4926,-0.8963,0.0000,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003]; % fval =1.1492e+003; AR =2.9015;
% x=[-36.4926,-0.8963,0.00001,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003]; % P_asterisk =8.0650;L =1.0000;AR =2.9015; Re =100.0018;
% x=[36.4926,0.8963,0.00001,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003]; % P_asterisk =11.8987; L = 0.7799; Re =100.0018;
% x=[41.3222,0.8963,0.001,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003];% P_asterisk =17.2758; L =1.0000; Re =113.2364;
% % x =[36.9382,0.9944,0.3065,1.1108,26.1599,0.4563,3.8949,1.9631,1.7982,0.0487];%P_asterisk =11.2214; L =0.9999; delta =-8.8309e-005; AR =2.9; Re =110.9498;
% x =[36.94,0.9944,0.3065,1.1108,26.1599,0.4563,3.8949,1.9631,1.7982,0.0487]; % P_asterisk =11.2231; L =1.0000; delta =9.1452e-006; AR =2.9; Re =110.9552;
% % x =[36.4001,1.437,0.3706,0.8834,26.4481, -0.7086,3.9809,1.964,1.7827,0.0000];% P_asterisk =15.5950;L =1.0001;delta =5.1088e-005;AR =2.9346;Re =160.0269;
% x =[36.4,1.4367,0.3706,0.8834,26.4481, -0.7086,3.9809,1.964,1.7827,0.0000]; % P_asterisk =15.5910; L =1.0000; delta =-4.1946e-005; AR =2.9346; Re =159.9930;
% x =[-31.0614,1.5655,0.0339,1.0516,9.2061,1.9760,3.9880,1.9642,1.9748,0.0145]; % P_asterisk =5.0899; L =1.0001; delta =5.0325e-005; AR =3.0358; Re =-153.9247;
%x=[49.0006,1.4463,0.0759,0.5151,0.3459,-1.1952,3.8347,1.9398,1.9360,0.0104];% P_asterisk =41.1636; L =1.0000; delta =3.8690e-005; AR =2.9749; Re =214.4092;
%x=[50.5894,1.5061,0.3189,0.5494,0.8681,-1.3308,3.4759,1.7892,1.9606,0.0113];% P_asterisk =34.4173;L =1.0000;delta =-5.2861e-006; AR =3.0385; Re =200.3041;
%x=[50.6778,1.1762,0.5131,0.8191,2.5339,-1.4619,3.6447,1.8919,1.9353,0.0000];%P_asterisk =12.7395;L =1.0000delta =-6.0800e-008; AR =2.9494; Re =170.0776;
% x=[50.6778,1.1762,0.5131,0.8191,2.5339,-1.4619,3.6447,1.8919,1.9353,0.0000];%P_asterisk=12.7395;L=1.0000delta =-6.0800e-008; AR =2.9494; Re =170.0776;下一行为精确值
% x=[5.0677789e+001  1.1762423e+000  5.1307852e-001  8.1914557e-001  2.5339154e+000 -1.4619157e+000  3.6447180e+000  1.8919010e+000  1.9353045e+000  1.2250774e-007];
% x=[50.2509,1.1885,0.5673,0.8808,2.6712,-1.5708,3.6218,1.8656,1.9472,-0.0000];% P_asterisk=10.6327;L=1.0000;delta =-3.1183e-005;AR =2.9851;Re =167.7027;  下一行为精确值
% x =[5.0250945e+001  1.1885097e+000  5.6727703e-001  8.8076414e-001  2.6712153e+000 -1.5707963e+000  3.6218464e+000  1.8655888e+000  1.9471743e+000 -7.0895123e-007];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[188.0,1.125,0.01,1.2201,5,-1.5427,3.004,0.8854,0.3289,0.356737]; % P_asterisk =9.5682; L =1.0045;
% x=[187.9159,1.1926,0.0154,1.2716,5.0,-1.5167,2.8477,0.9105,0.3713,0.3225];% 下一行为精确值
% x=[1.8791586e+002,1.1926443e+0,1.5415423e-002,1.2716489e+0,5.0000000e+0,-1.5166711e+0,2.8477116e+0,9.1047551e-001,3.7127239e-001,3.2250352e-001];
% x =[48.4556,1.1912,0.6248,0.9703,2.9299,-1.4950, 3.6133,1.9328,1.9920,-0.0000];P_asterisk =8.7136; L =1.0000 delta =-4.7214e-008; AR =2.9000;Re =169.0186; % 下一行为精确值
% x=[4.8455586e+001,1.1912465e+0,6.2476075e-001,9.7025351e-001,2.9298698e+0,-1.4949918e+0,3.6132652e+0,1.9328419e+0,1.9919775e+0,-2.6197345e-006];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_start=[46.545,1.30898,0.0154,1.2716,2.5,-1.5167,3.9809,1.964,1.7827,0.0000];  % P_asterisk =4.7922; L =1.0020;
% x =[46.6104,1.3076,0.1891,1.2715,2.5103,-1.5169,3.9813,1.96,1.7845,0.0001]; % P_asterisk =4.7867;L =1.0000; delta = -2.9524e-008; AR =2.9417; Re =186.1461;% 下一行为精确值
x =[4.6610447e+001,1.3075639e+0,1.8911207e-001,1.2714657e+0,2.5103035e+0, -1.5168925e+0,3.9812634e+0,1.9600006e+0,1.7844515e+0,7.9074669e-005];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[3.004,0.8854,0.3289,0.356737];
% % 下面的数据是实测原始翅膀形貌和翅膀运动计算而得
% % 升阻比、功率比、功率密度以及升重比、展弦比和雷诺数。
% Rario_F_vertical_horiz_aver =1.3142;    % 升推比
% F_vertical_to_horiz_rms =0.8448;          % 升推比rms
% P_aero_aver_to_inertia_aver =13.5511; % 全正气动功率与全正惯性功率的比值
% P_asterisk =22.3727;  % 功率密度
% L =1.3319;                  % 升重比
% AR =3.7643;               % 展弦比
% Re =155.8311;           % 雷诺数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
Force_Power=Aero_F_M_fruitfly(x); % 调用函数: 求垂直方向的平均力和翅运动平均正功率
% Force_Power=Aero_F_M_fruitfly();      % 调用函数: 求垂直方向的平均力和翅运动平均正功率
F_vertical_aver=2*Force_Power(1);     % 两个翅膀——垂直方向的平均力
P_total_aver=2*Force_Power(2);         % 两个翅膀——分别两个自由度翅运动产生的正功率
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_insect=1.8;         % M_fly=1.8mg;  % 2014-Science-MH Dickinson     % P_asterisk =31.8058; % 量纲单位=uW/mg;
% m_insect=0.72;    % 0.72 mg  % 2007-JFM-Wang ZJ                           % P_asterisk =80.6538; % 量纲单位=uW/mg;
% m_insect=0.96;    % 0.96 mg  % 2005-JEB-MH Dickinson                    % P_asterisk =60.4903; % 量纲单位=uW/mg;
% m_insect=1.05;    % 1.05mg; % 1997-JEB-Fritz-Olaf lehmann              % P_asterisk =55.3054; % 量纲单位=uW/mg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x是变量variable=[R_wingeff,C_avereff,xr,C_maxy]; 
% P_total_aver=Aero_F_M_fruitfly(x);  %调用函数: 正功率——平动功率和扭转功率——扭转轴功率和拍打轴功率
P_asterisk=P_total_aver/m_insect       % 功率密度——量纲单位=uW/mg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二部分——penaltyfun1——升重比约束——惩罚函数
% % x=[3.004,0.8854,0.3289,0.4644];  % 未优化的果蝇翅膀形貌参数——翅根和翅尖连线的扭转轴定出的前缘
% x=[3.004,0.8854,0.3289,0.356737]; % 初始值. % 未优化的果蝇翅膀形貌参数—翅根和翅尖连线的扭转轴定出的前缘/最大前缘点和最新后缘点坐标间距=0.3567
% % x=[3.004,0.8854,0.3289,0.1389];  % 未优化的果蝇翅膀形貌参数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2014-Science-实验测试而得运动学数据；翅膀形貌学参数为：
% % x=[3.004,0.8854,0.3289,0.25]; %初始值. %未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时——结果
% % P_asterisk=31.8058; L=1.3940; % 未优化的果蝇翅膀形貌参数—扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时——结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @ options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400);
% x =[3.6470,1.9176,0.6487,0.4258]; 
%@options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.3086,0.9363,0.4644]; 
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fminsearch,fminsearchOptions});
% x =[4.0000,1.2895,0.9238,0.4644]; %ConstraintFunction=@aspectratio_constraint;  
% @options=gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'HybridFcn',{@fmincon,fminconOptions});
% x =[2.0000,1.1161,1.0711,-0.1865];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % 不合理
% x=[4.0000,1.4283,0.9992,0.5000];  % L =-2.6131;  % 不合理
% x =[4,2,0.8188,0.5];                       % L = 1.0009;
% x=[2.2899,1.0561,0.3788,0.0001];  % 20150204-Hybrid GA-Fmincon计算结果
% x =[3.1133,0.5938,0.3297,0.2450]; % 20150204-Hybrid GA-Fminsearch计算结果
% x=[2.5088,1.1424,0.7124,0.0458];   % 20150209-Hybrid GA-Fminsearch计算结果
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;   % N*kg^-1
m_insect=1.8; 
M_insectweight=m_insect*g;                  % M_weight=10.3±1.27uN@1997-JEB-Fritz-Olaf lehmann
% L会改变,因为扭转轴的位置为改变转动虚质量力系数I6z的大小;
% 此外还会引起阻尼力矩系数(Z_rnd,M_xrdcoeff)和虚质量力矩系数(I_xzam,I_xxam)的变化;
% x是变量variable=[R_wingeff,C_avereff,xr,C_maxy]; % F_verticalaver=12.3074uN; 
% F_vertical_aver=Aero_F_M_fruitfly(x); %调用函数:
% F_verticalaver=trapz(t,F_vertical)/(3*T);
L=F_vertical_aver/M_insectweight     % 升重比——无量纲: L =1.3940; 
delta=L-1                                        % delta =-1.0332e-004;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 构建惩罚函数——penaltyfun1
r=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (L-1)<2.9*10^(-15)&&(L-1)>0
    penaltyfun1=0;
% elseif L<0 
%     penaltyfun1=r*abs((L-1));
elseif (L-1)>=2.9*10^(-15)||(L-1)<=0
    penaltyfun1=r*abs((L-1));
else
    disp('垂直方向力有错');
end
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 该方案不可行的哦
% if ((L-1)>(2.9*10^(-15)) || (L-1)<0)
%     penaltyfun1=r*(L-1);    % penaltyfun1 =-1.6561e+004; 
% else
%    penaltyfun1=0;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% penaltyfun1=r*heaviside(1-L);        % heaviside(x) has the value 0 for x < 0;         1 for x > 0;       and 0.5 for x = 0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第三部分——penaltyfun2——优化参数变量的上下界约束——惩罚函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x=[4.0000,0.5000,1.3364,0.5000];  % 不合理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcn_con=sum_constraint(x);
% 构建惩罚函数——penaltyfun2
s=2000;
penaltyfun2=s*fcn_con;  % penaltyfun2 =663.6000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 第四部分——penaltyfun3——展弦比约束或者Rossby数约束——惩罚函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% penaltyfun3=aspectratio_constraint(x);  % 展弦比约束(即Rossby值约束)
% penaltyfun4=Re_constraint(x);                % Re值约束
% obj_function=P_asterisk+penaltyfun1+penaltyfun2+penaltyfun3+penaltyfun4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% obj_function=P_asterisk+penaltyfun1+penaltyfun2;
penaltyfun5=AR_Re_constraint(x);       % 展弦比约束(即Rossby值约束) 和 Re数约束
obj_function=P_asterisk+penaltyfun1+penaltyfun2+penaltyfun5;
% toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%