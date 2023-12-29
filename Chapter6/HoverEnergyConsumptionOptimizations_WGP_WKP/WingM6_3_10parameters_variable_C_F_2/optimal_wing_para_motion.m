% optimal_wing_para_motion
%% 第一部分——变量的约束边界
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 翅膀长度和展弦比——R_wingeff,C_avereff和AR
% R_wingeff=3.004;      % 有效翅膀长度(mm)  
% C_avereff=0.8854;     % mm
% (2) 可以考虑已知展弦比, 求解气动力最优的翅形貌
% AR=R_wingeff/C_avereff;     %aspect ratio: AR=R^2/A_w;  这里输出为: AR=3.40158;  % Science数据是: AR=3.1519;
% C_avereff=R_wingeff/AR;     %mean chord length: C_aver=A_w/R=R^2/AR/R=R/AR; % C_aver=0.884mm;
% A_w=R_wingeff^2/AR;        %Area of wing: 这里输出为: A_w=2.66mm^2;   %RJ Wood设计的翅膀: A_w=2.84345 mm^2 
% (3) 翅根部偏离距离——xr
% xr=0.3289;                     % x-root offset  \mm
% xr_nd=xr/R_wingeff;      % x-root offset  无量纲展向偏置距离
% (4) 扭转轴的位置——C_maxy
% C_max_LtoT=C_lead_ymax-C_trail_ymin;       % C_max_LtoT =1.3018;
% x_0lb=0*C_max_LtoT;                                    % x_0lb=0;
% x_025=0.25*C_max_LtoT;                              % x_025=0.3255;
% x_0ub=0.5*C_max_LtoT;                                % x_0ub=0.6509;
% C_lead_ymax=max(f_x_lead2); % 输出: C_lead_ymax=0.4644; k_leadmax=644;%原始果蝇翅膀@最大前缘点到翅根翅尖连线的距离;@C_maxy=0;
% % C_maxylb =0.464385778290230;
% C_maxy25 =0.138924474377504;  % 针对程序wing_model_88_yaxis有: 第122行; C_maxy =0.1389; 
% C_maxyub =-0.186536829535222; 
% C_maxy=C_maxylb;
% C_maxy=C_maxy25;
% C_maxy=C_maxyub;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable=[R_wingeff,C_avereff,xr,C_maxy];        % 变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5) 变量的约束边界——下面的排序依次是R_wingeff, C_avereff, xr和x_0(扭转到前缘最大值点的距离)
% nvars = 4;         % Number of variables 
% R_wingeff∈[2,4]*10^(-3);
% C_avereff∈[0.5,2]*10^(-3);
% xr∈[0,2]*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0∈[0,0.5];  % 扭转轴距前缘的无量纲距离——% 原始果蝇翅膀 R_wing=3.004;  C_aver=0.8854; 
% C_maxylb =0.464385778290230;                   % C_maxylb=0.4644;
% C_maxy25 =0.138924474377504;                  % C_maxy25=0.1389; 
% C_maxyub =-0.186536829535222;                % C_maxyub=-0.1865
%%%%%%%%%%%%%%%%%%%%%%%%%
% % variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 变量
% % f∈[0,inf];  phi_m∈[0,pi/2];  K∈[0,1]; % eta_m∈[0,pi]; C_eta∈[0,inf];  Phi_eta∈[-pi,pi]; eta_0∈[eta_m-pi,pi-eta_m];
% eta_m=x(8);
% eta_0min=eta_m-pi;
% eta_0max=pi-eta_m;
% LB = [2, 0.5, 0, 0,0,0,0,0,0,-pi,eta_0min];            % Lower bound       % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
% UB = [4, 2, 2, 0.5,inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound      % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二部分——目标函数——非线性不等式约束-Hybrid-GA优化
tic
% matlabpool open 12
% matlabpool open local2 12
ObjectiveFunction=@morphology_para_opt_obj_function1;  % fitness and constraint functions——调用目标函数
nvars=10;   % Number of variables 
% % variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 变量
% % f∈[0,inf];  phi_m∈[0,pi/2];  K∈[0,1]; % eta_m∈[0,pi]; C_eta∈[0,inf];  Phi_eta∈[-pi,pi]; eta_0∈[eta_m-pi,pi-eta_m];
% eta_m=x(8);
% eta_0min=eta_m-pi;
% eta_0max=pi-eta_m;
% LB = [2, 0.5, 0, 0,0,0,0,0,0,-pi,eta_0min];            % Lower bound       % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
% UB = [4, 2, 2, 0.5,inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound      % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第三部分——nonlinear inequality constraints——调用非线性不等式约束
% X_AspR=[R_wingeff,C_avereff,xr];               % AR=(R_wingeff+xr)/C_avereff; % AR∈[2.5,3.5]
% ConstraintFunction=@aspectratio_constraint;  % ——展弦比对变量的线性约束
% NonLConstraint=@NonL_constraint;                    % ——Re雷诺数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第四部分——混合遗传算法(Hybrid Function)-GA+fminsearch
% 原始果蝇翅膀 R_wing=3.004;  C_aver=0.8854;   
% C_lead_ymax =0.4644;  C_trail_ymin =-0.8374;  C_max_LtoT = 1.3018; @C_maxy=0;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %针对翅根和翅尖连线的扭转轴定出的前缘
% x_start=[3.004,0.8854,0.3289,0.356737]; % 初始值. % 未优化的果蝇翅膀形貌参数——翅根和翅尖连线的扭转轴定出的前缘
% C_lead_ymax =0.3255;  C_trail_ymin =-0.9764;  C_max_LtoT =1.3018;
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.3255/1.3018=0.25;%针对扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时
% x_start=[3.004,0.8854,0.3289,0.25];     % 初始值. % 未优化的果蝇翅膀形貌参数——扭转轴位于最大前缘点和最新后缘点坐标间距0.25倍时
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 纯粹遗传算法最小化优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = gaoptimset('InitialPopulation',x_start);  % (1)——可行——第四次计算(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'Display','iter');  % (2)——可行
% options = gaoptimset(options,'InitialPopulation',x_start,'TimeLimit',14400);  % (2)——可行——第四次计算(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 混合遗传算法(Hybrid Function)-GA+fminsearch
% % 自动由fminsearch跳转至fmincon搜索——可行——结果不好——第四次计算(3)
% options = gaoptimset('InitialPopulation',x_start,'TimeLimit',10800,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr});
% % fminsearchOptions = optimset('algorithm','Nelder-Mead simplex direct search');
% % options=gaoptimset(options,'HybridFcn',{@fminsearch,x_start, fminsearchOptions}); 
% fminsearchOptions=optimset('Display','iter'); % fminsearch——无约束多变量函数最小化
% options=gaoptimset(options,'HybridFcn',{@fminsearch,fminsearchOptions}); % fminsearch——无约束多变量函数最小化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 采用fminsearch局部搜索——不可行
% options=gaoptimset('PopulationSize',200,'InitialPopulation',[],'TimeLimit',10800,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr});
% options=gaoptimset('PopulationSize',100,'InitialPopulation',x_start,'TimeLimit',10800,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr}); %可行
% fminsearchoptions = optimset('Display','iter','TolCon',1e-10);% fminsearch——无约束多变量函数最小化 %可行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options=gaoptimset('PopulationSize',100,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'UseParallel','always');%并行计算——可行
% fminsearchoptions = optimset('Display','iter','TolCon',1e-10,'TolFun',1e-6,'UseParallel','always');%fminsearch—无约束多变量函数最小化 %并行计算——可行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 方案1——不设置初始值
% x=[45.4147,1.5696,1.8489,1.0385,-0.9685,-0.0463,3.5193,1.8689,1.3038,-0.0008]; % 初始种群1
% P_asterisk =6.9607;     L =1.0002;   delta =1.6715e-004;   penaltyfun1 =0.3343; 
% penaltyfun2 =4.5714;   AR =2.5807;   Re =126.7144;     penaltyfun5 =0;
% x=[45.4271,1.5679,1.8515,1.0364,-0.9666,-0.0470,3.5258,1.8657,1.3081,0.0122];  % 初始种群2
% P_asterisk =7.0267;   L =1.0000; delta =2.6792e-005;   penaltyfun1 =0.0536;
% penaltyfun2 =0;    AR =2.5909;   Re =126.6284;    penaltyfun5 =0;
% final_pop=x;  % options = gaoptimset('InitialPop', final_pop);
% options=gaoptimset('InitialPop', final_pop,'Display','iter','PopulationSize',100,'Generations',300,... 
% options=gaoptimset('Display','iter','PopulationSize',100,'Generations',500,...   % 回复继续搜索，将代数=200修改为代数=5+10
%     'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%     'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always'); %并行计算——可行
% options=gaoptimset('Display','iter','PopulationSize',100,'Generations',200,...   % 回复继续搜索，将代数=200修改为代数=5+10
%     'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%     'UseParallel','always'); %并行计算——可行
% 方案2——设置初始值
% x_start=[2.5088,1.1424,0.7124,0.0458,108.5082,1.7828,3.0150,0.7854,1.0431,0.6109,0.0749];
% x_start=[3.1935,1.6054,2.1531,0.5000,1.2853,0.8770,0.6763,0.1370,1.5444,-2.0757,0.2672]; %20150315-hybrid_GA_fminsearch_WingM4_3-11变量-结果
% options=gaoptimset('PopulationSize',150,'InitialPopulation',x_start,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'UseParallel','always'); %并行计算——可行
% x_start=[-29.6664,0.0000,1.0278,1.5708,25.9878,-1.2863,3.6594,1.7727,1.9991,0.0002]; % 遗传350代+2000代局部搜索的结果
% x_start=[-29.6676,0.1662,0.3033,1.5707,25.9866,0.0444,3.6592,1.7728,1.9995,0.0002];  % 继续遗传50代+2000代局部搜索的结果
% x_start=[-36.4926,-0.8963,0.0000,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003]; % 继续遗传100代+3000代局部搜索的结果
% x_start=[-36.8222,-0.8560,-0.0005,0.8286,27.7493,-0.3343,4.0041,1.9705,2.0000,0.0019]; % 继续遗传20代+2000代局部搜索的结果
% x_start=[36.4926,0.8963,0.0001,0.8794,26.3600,-0.1500,3.9772,1.9745,1.7520,0.0003]; % 由上上一行指令修改而来+继续遗传20
% x_start =[12.658,0.0675,1.1021,0.9084,16.4536,-0.0232,2.0000,1.9935,1.0759,0.0214];
% x_start=[5.7434,1.309,1.0386,1.3407,3.4500,-1.0554,3.2549, 1.9932,1.8395,0.2887];
% x_start=[188.0,1.155,0.01,1.2201,5,-1.5427,3.004,0.8854,0.3289,0.356737]; % P_asterisk =9.8620; L =1.0050;  % 继续遗传50代+2000代局部搜索的结果
x_start=[60.255,1.30898,0.01,1.2201,2.5,-1.5427,3.9809,1.964,1.7827,0.0000];  % P_asterisk=6.8426; L =1.0013;
options=gaoptimset('Display','iter','PopulationSize',100,'InitialPopulation',x_start,'Generations',50,...   % 回复继续搜索，将代数=200修改为代数=5+10
    'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
    'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always'); %并行计算——该方案种群过大
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_start=[33.4675,0.9944,0.3065,1.1108,26.1599,-pi/4,3.8949,1.9631,1.7982,0.0487]; % 接力继续遗传10代+2000代局部搜索的结果
% options=gaoptimset('Display','iter','PopulationSize',100,'InitialPopulation',x_start,'Generations',10,...   % 回复继续遗传5代
%     'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%     'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always'); %并行计算——可行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options=gaoptimset('Display','iter','PopulationSize',100,'InitialPopulation',x_start,'Generations',20,...   % 回复继续遗传50代+继续遗传100+继续遗传20
%     'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%     'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always'); %并行计算——可行
% options=gaoptimset('Display','iter','PopulationSize',100,'InitialPopulation',x_start,'Generations',4,...   % 回复继续遗传50代+继续遗传100+继续遗传20
%     'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
%     'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always'); %并行计算——可行
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fminsearchoptions =
% optimset('Display','iter','MaxFunEvals',14000,'TolCon',1e-10,'TolFun',1e-15,'UseParallel','always'); % 可行feasible
% fminsearchoptions=optimset('Display','iter-detailed','MaxFunEvals',2500,'MaxIter',2500,'TolX',1e-15,'TolFun',1e-15,'UseParallel','always');% 初始搜索1
% fminsearchoptions=optimset('Display','iter-detailed', 'MaxFunEvals',3500,'MaxIter',3500,'TolX',1e-15,'TolFun',1e-15,'UseParallel','always');% 初始搜索2
% fminsearchoptions=optimset('Display','iter-detailed','TolX',1e-15,'TolFun',1e-15,'UseParallel','always'); % 初始搜索3
% fminsearchoptions=optimset('Display','iter-detailed','MaxFunEvals',3000,'MaxIter',3000,'TolX',1e-15,'TolFun',1e-15,'UseParallel','always');% 初始搜索4-迭代次数过大
% fminsearchoptions=optimset('Display','iter-detailed','MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-15,'TolFun',1e-15,'UseParallel','always');% 初始搜索5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fminsearchoptions=optimset('Display','iter-detailed', 'MaxFunEvals',2000,'MaxIter',2000,'TolX',1e-15,'TolFun',1e-15,'UseParallel','always');% 初始搜索6                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fminsearch——Find minimum of unconstrained multivariable function using derivative-free method
options =gaoptimset(options,'Hybridfcn',{@fminsearch,fminsearchoptions}); % 寻找无约束多变量最小值——不可取
% [x,fval,exitflag,output,population,scores]=ga(ObjectiveFunction,nvars,[],[],[],[],[],[],NonLConstraint,options)% ——不可取
[x,fval,exitflag,output,population,scores]=ga(ObjectiveFunction,nvars,[],[],[],[],[],[],[],options)
% matlabpool close
matlabpool close force local2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wing=x(7);   % mm
C_aver=x(8);    % mm
xr0=x(9);
R_wing1=R_wing*10^-3;  %   R_wing=3.004*10^-3;            % m
C_aver1=C_aver*10^-3;   %  C_aver=0.8854*10^-3;            % m
xr01=xr0*10^-3;
f=x(1);               % f=188.7;                  % Hz  
phi_aver=x(2);   % phi_aver = 1.1488;  % rad   % phi_aver*180/pi  % ans =65.8189;
PHI=2*phi_aver;
nu=1.48*10^-5;      % 空气运动粘度 m^2/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%
AR=(R_wing+xr0)/C_aver
Re=(2*PHI*(R_wing1+ xr01)*f*C_aver1)/(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optimal_variablex=x;
object_value=fval;
save('result_x.txt', 'optimal_variablex', '-ASCII')
save('result_fval.txt', 'object_value', '-ASCII')
output_1=[population,scores];
% xlswrite('D:\cal_result\WingM6_3_10variable_group_Wang_wingbeatM_variable_C_F\population_scores.xlsx',output_1,'sheet1','A1:K100');
xlswrite('population_scores.xlsx',output_1,'sheet1','A1:K100');
elapsedTime =toc/3600;  
disp(['Caculation finished and the elapsedTime= ' num2str(elapsedTime ,'%5.4f') 'hours'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% 第三部分——目标函数——含边界、线性等式和非线性不等式约束-Hybrid-GA优化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 采用fminunc局部搜索——不可行
% options = gaoptimset('InitialPopulation',x_start,'TimeLimit',14400,'PlotFcns', {@gaplotbestf,@gaplotstopping});
% fminuncOptions=optimset('Display','iter', 'LargeScale','off'); % fminunc——无约束多变量函数最小化
% options=gaoptimset(options,'HybridFcn',{@fminunc, fminuncOptions});  % fminunc——无约束多变量函数最小化
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 采用fmincon局部搜索——可行——结果很好——第四次计算(4)
tic
matlabpool open 8
ObjectiveFunction=@morphology_para_opt_obj_function2;  % fitness and constraint functions——调用目标函数
nvars=11;   % Number of variables 
LB =[150, 0, -pi, -pi/4, 0, -pi,-pi/4, 2, 0.5, 0, 0];              % Lower bound   % 11个变量
UB=[260, pi/2, pi, pi/4, pi/2, pi, pi/4, 4, 2, 2, 0.5];          % Upper bound   % 11个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第三部分——nonlinear inequality constraints——调用非线性不等式约束
% X_AspR=[R_wingeff,C_avereff,xr];               % AR=(R_wingeff+xr)/C_avereff; % AR∈[2.5,3.5]
% ConstraintFunction=@aspectratio_constraint;  % ——展弦比对变量的线性约束
NonLConstraint=@NonL_constraint;                    % ——Re雷诺数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  A=[0, 0, 0, 0, 0, 0, 0, 1, -3.5, 1, 0; ...  % for R-3.5*C_aver+xr0<=0;  % XXX
%       0, 0, 0, 0, 0, 0, 0, -1, 2.5, -1, 0];    % for -R+2.5*C_aver-xr0<=0;  % XXX
 A=[0, 0, 0, 0, 0, 0, 0, 1, -3.5, 0, 0; ...  % for R-3.5*C_aver<=0;
      0, 0, 0, 0, 0, 0, 0, -1, 2.5, 0, 0];    % for -R+2.5*C_aver<=0;
 b=[0; 0];
x_start=[188.7,1.5400,2.8107,0, 0.9935,0.0114, 0, 2.5088,1.1424,0.7124,0.0458];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = gaoptimset('Display','iter','PopulationSize',100,'InitialPopulation',x_start,'TimeLimit',10800,...
    'PlotFcns',{@gaplotbestf,@gaplotmaxconstr,@gaplotbestindiv,@gaplotstopping},'UseParallel','always');
% fminconOptions=optimset('Display','iter','TolCon',1e-10); % fmincon——有约束非线性多变量最小化——算法不选择
fminconOptions=optimset('Display','iter','TolCon',1e-15,'TolFun',1e-15); % fmincon——有约束非线性多变量最小化——算法不选择
% fminconOptions=optimset('Display','iter','Algorithm','active-set'); % fmincon——有约束非线性多变量最小化——算法选择
options=gaoptimset(options,'HybridFcn',{@fmincon,fminconOptions}); % fmincon——有约束非线性多变量最小化
[x,fval,exitflag,output,population,scores]=ga(ObjectiveFunction,nvars,A,b,[],[],LB,UB,NonLConstraint,options)
matlabpool close
toc
%%
% output_1=[population,scores];
% xlswrite('D:\kxj\WingM5_2_10variable_group\population_scores.xlsx',output_1,'sheet1','A1:L100');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next we run the GA solver.
% [x,fval]=ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[],options)
% output_1=[population,scores];
% xlswrite('D:\KXJ\hybrid_GA_fminsearch_WingM4_4_2\population_scores.xlsx',output_1,'sheet1','A1:L100');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%