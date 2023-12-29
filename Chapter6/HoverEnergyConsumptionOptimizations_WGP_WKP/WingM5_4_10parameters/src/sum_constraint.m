function y=sum_constraint(x)
% 建立约束之外的惩罚函数
% variable=[f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 7个变量
% variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];   % 11个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 变量的约束边界——下面的排序依次是R_wingeff, C_avereff, xr和x_0(扭转到前缘最大值点的距离)
% nvars = 4;         % Number of variables 
% R_wingeff∈[2,4]*10^(-3);
% C_avereff∈[0.5,2]*10^(-3);
% xr∈[0,2]*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0∈[0,0.5];  % 距前缘的无量纲距离
% C_maxylb =0.464385778290230;    % C_maxylb =0.4644;
% C_maxy25 =0.138924474377504;   % C_maxy25 =0.1389; 
% C_maxyub =-0.186536829535222; % C_maxyub =-0.1865
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二方案————variable=[f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 7个变量
% % variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];   % 11个变量
% % f∈[0,inf];  phi_m∈[0,pi/2];  K∈[0,1]; % eta_m∈[0,pi]; C_eta∈[0,inf];  Phi_eta∈[-pi,pi]; eta_0∈[eta_m-pi,pi-eta_m];
% eta_m=x(4);
% eta_0min=eta_m-pi;
% eta_0max=pi-eta_m;
% % LB = [2, 0.5, 0, 0,0,0,0,0,0,-pi,eta_0min];            % Lower bound       % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]   % 11个变量
% % UB = [4, 2, 2, 0.5,inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound      % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]   % 11个变量
% LB = [0,0,0,0,0,-pi,eta_0min];                 % Lower bound   % 7个变量
% % UB = [inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound   % 7个变量
% UB = [400,pi/2,1,pi/2,200,pi,eta_0max];   % Upper bound   % 7个变量——注意3个数据修改上限,不同于2007-JFM-Wang ZJ
% % X_AspR=[R_wingeff,C_avereff,xr];     % AR=(R_wingeff+xr)/C_avereff; % AR∈[2,8]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第四方案——翅膀参数(4个变量)—— 人为设计谐波运动(harmonic_motion)(7个变量)——(f,phi_m,epsilon,phi_0,psi_m,zeta,psi_0)
% (1) 翅膀形貌学参数+翅膀运动参数
% LB = [2, 0.5, 0, 0,0,0,-pi,-pi/4,0,-pi,-pi/4];                   % Lower bound   % 11个变量
% UB=[4, 2, 2, 0.5,400,pi/2,pi,pi/4,pi/2,pi,pi/4];             % Upper bound   % 11个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 翅膀运动参数+翅膀形貌学参数——(f,phi_m,epsilon,psi_m,zeta)+(R_wing,C_aver,xr0,C_maxyaxis)
% LB = [0,0,-pi,0,-pi,2, 0.5, 0, 0];              % Lower bound   % 9个变量
% UB=[300,pi/2,pi,pi/3,pi,4, 2, 2, 0.5];          % Upper bound   % 9个变量
% LB =[150, 0, -pi, -pi/4, 0, -pi,-pi/4, 2, 0.5, 0, 0];              % Lower bound   % 11个变量
% UB=[260, pi/2, pi, pi/4, pi/2, pi, pi/4, 4, 2, 2, 0.5];          % Upper bound   % 11个变量
% (f,phi_m,epsilon,psi_m,zeta,psi_0,R_wing,C_aver,xr0,C_maxyaxis) % 10个变量
LB =[0,     0,    -pi, 0,    -pi,  -pi/4, 2, 0.5, 0, 0];              % Lower bound   % 10个变量
UB=[300, pi/2, pi, pi/2, pi,   pi/4,  4, 2,   2, 0.35];          % Upper bound   % 10个变量
%x=[35.4018,1.5769,-0.5250,-0.3650,0.8572,-3.1414,0.0167,3.9583,1.8190,1.9995,-0.0486]; % 11个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
con_min=LB;
con_max=UB;
zeta=zeros(N,1);
% y=zeros(N,1);
for i=1:N
    if x(i) < con_min(i)
        zeta(i)=abs(con_min(i)-x(i))/(con_max(i)-con_min(i));
    elseif x(i) > con_max(i)  % x=[4.0000,0.5000,1.3364,0.5000];  penaltyfun2=663.6000=2000*sum(zeta);则y=sum(zeta)=0.3318;
        zeta(i)=abs(x(i)-con_max(i))/(con_max(i)-con_min(i));
    else % x(i)>=con_min(i) && x(i)<=con_max(i);  % 无需这句表达式
        zeta(i)=0;  % disp('变量未超出边界约束');  
    end
end
 y=sum(zeta);
end