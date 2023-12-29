function penaltyfun5=AR_Re_constraint(x)
% variable=[f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0,R_wingeff,C_avereff,xr,C_maxy];    % 输入11个变量
% 2009-JEB-David Lentink: AR∈[1,5], 最好限制在AR∈[2.5,3.5]
% 2004-JEB-Dickinson: Re∈[100,1000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x =[2.8632;0.8950;0.1133;0.3662];
R_wing=x(7);   % mm
C_aver=x(8);    % mm
xr0=x(9);
% C_maxy=x(10);
%%%%%%%%%%%%%%%%%%
% R_wing=x(1);   % mm
% C_aver=x(2);    % mm
% xr0=x(3);
% % C_maxy=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_wing1=R_wing*10^-3;  %   R_wing=3.004*10^-3;            % m
C_aver1=C_aver*10^-3;   %  C_aver=0.8854*10^-3;            % m
xr01=xr0*10^-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % LB = [2, 0.5, 0, 0];      % Lower bound       % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
% UB = [4, 2, 2, 0.5];      % Upper bound      % 单位: mm & 扭转轴距前缘的无量纲距离[0,0.5]
% LB =[150, 0, -pi, -pi/4, 0, -pi,-pi/4, 2, 0.5, 0, 0];              % Lower bound   % 11个变量
% UB=[260, pi/2, pi, pi/4, pi/2, pi, pi/4, 4, 2, 2, 0.5];          % Upper bound   % 11个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=x(1);               % f=188.7;                  % Hz  
phi_aver=x(2);   % phi_aver = 1.1488;  % rad   % phi_aver*180/pi  % ans =65.8189;
PHI=2*phi_aver;
nu=1.48*10^-5;      % 空气运动粘度 m^2/s
%  x=[3.004,0.8854,0.3289,0.25]; % AR=3.76429;
% AR=(R_wing1+ xr01)/C_aver1     % AR =2.2102;   % AR =3.3928(XXX);
% % Re1=(2*PHI*R_wing1^2*f)/(nu*AR);   % Re1 =135.3017;
% Re=(2*PHI*(R_wing1+ xr01)*f*C_aver1)/(nu)      % Re =210.9388;  % Re∈[100,1000];  % 2004-JEB-Dickinson: Re∈[100,1000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (R_wing-3.5*C_aver)>0 || (-R_wing+2.5*C_aver)>0 || (100*nu-2*PHI*f*R_wing1*C_aver1)>0 || (2*PHI*f*R_wing1*C_aver1-1000*nu)>0
% if (R_wing+xr0-3.5*C_aver)>0 || (-R_wing-xr0+2.5*C_aver)>0 || (100*nu-2*PHI*f*R_wing1*C_aver1)>0 || (2*PHI*f*R_wing1*C_aver1-1000*nu)>0
 if (R_wing+xr0-5*C_aver)>0 || (-R_wing-xr0+2.9*C_aver)>0 || (100*nu-2*PHI*f*(R_wing1+xr01)*C_aver1)>0 || (2*PHI*f*(R_wing1+xr01)*C_aver1-3000*nu)>0 
    penaltyfun5=2000;  % 构建惩罚函数——penaltyfun5  ——% disp('展弦比约束失效(即Rossby值)或者Re值不合理');
else
    penaltyfun5=0;
end