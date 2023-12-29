function y=sum_constraint(x)
% ����Լ��֮��ĳͷ�����
% variable=[f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 7������
% variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];   % 11������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������Լ���߽硪�����������������R_wingeff, C_avereff, xr��x_0(Ťת��ǰԵ���ֵ��ľ���)
% nvars = 4;         % Number of variables 
% R_wingeff��[2,4]*10^(-3);
% C_avereff��[0.5,2]*10^(-3);
% xr��[0,2]*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% x_0��[0,0.5];  % ��ǰԵ�������پ���
% C_maxylb =0.464385778290230;    % C_maxylb =0.4644;
% C_maxy25 =0.138924474377504;   % C_maxy25 =0.1389; 
% C_maxyub =-0.186536829535222; % C_maxyub =-0.1865
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ�������������variable=[f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];     % 7������
% % variable=[R_wingeff,C_avereff,xr,C_maxy,f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0];   % 11������
% % f��[0,inf];  phi_m��[0,pi/2];  K��[0,1]; % eta_m��[0,pi]; C_eta��[0,inf];  Phi_eta��[-pi,pi]; eta_0��[eta_m-pi,pi-eta_m];
% eta_m=x(4);
% eta_0min=eta_m-pi;
% eta_0max=pi-eta_m;
% % LB = [2, 0.5, 0, 0,0,0,0,0,0,-pi,eta_0min];            % Lower bound       % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]   % 11������
% % UB = [4, 2, 2, 0.5,inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound      % ��λ: mm & Ťת���ǰԵ�������پ���[0,0.5]   % 11������
% LB = [0,0,0,0,0,-pi,eta_0min];                 % Lower bound   % 7������
% % UB = [inf,pi/2,1,pi,inf,pi,eta_0max];   % Upper bound   % 7������
% UB = [400,pi/2,1,pi/2,200,pi,eta_0max];   % Upper bound   % 7����������ע��3�������޸�����,��ͬ��2007-JFM-Wang ZJ
% % X_AspR=[R_wingeff,C_avereff,xr];     % AR=(R_wingeff+xr)/C_avereff; % AR��[2,8]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ķ�������������(4������)���� ��Ϊ���г���˶�(harmonic_motion)(7������)����(f,phi_m,epsilon,phi_0,psi_m,zeta,psi_0)
% (1) �����òѧ����+����˶�����
% LB = [2, 0.5, 0, 0,0,0,-pi,-pi/4,0,-pi,-pi/4];                   % Lower bound   % 11������
% UB=[4, 2, 2, 0.5,400,pi/2,pi,pi/4,pi/2,pi,pi/4];             % Upper bound   % 11������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ����˶�����+�����òѧ��������(f,phi_m,epsilon,psi_m,zeta)+(R_wing,C_aver,xr0,C_maxyaxis)
% LB = [0,0,-pi,0,-pi,2, 0.5, 0, 0];              % Lower bound   % 9������
% UB=[300,pi/2,pi,pi/3,pi,4, 2, 2, 0.5];          % Upper bound   % 9������
% LB =[150, 0, -pi, -pi/4, 0, -pi,-pi/4, 2, 0.5, 0, 0];              % Lower bound   % 11������
% UB=[260, pi/2, pi, pi/4, pi/2, pi, pi/4, 4, 2, 2, 0.5];          % Upper bound   % 11������
% (f,phi_m,epsilon,psi_m,zeta,psi_0,R_wing,C_aver,xr0,C_maxyaxis) % 10������
LB =[0,     0,    -pi, 0,    -pi,  -pi/4, 2, 0.5, 0, 0];              % Lower bound   % 10������
UB=[300, pi/2, pi, pi/2, pi,   pi/4,  4, 2,   2, 0.35];          % Upper bound   % 10������
%x=[35.4018,1.5769,-0.5250,-0.3650,0.8572,-3.1414,0.0167,3.9583,1.8190,1.9995,-0.0486]; % 11������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
con_min=LB;
con_max=UB;
zeta=zeros(N,1);
% y=zeros(N,1);
for i=1:N
    if x(i) < con_min(i)
        zeta(i)=abs(con_min(i)-x(i))/(con_max(i)-con_min(i));
    elseif x(i) > con_max(i)  % x=[4.0000,0.5000,1.3364,0.5000];  penaltyfun2=663.6000=2000*sum(zeta);��y=sum(zeta)=0.3318;
        zeta(i)=abs(x(i)-con_max(i))/(con_max(i)-con_min(i));
    else % x(i)>=con_min(i) && x(i)<=con_max(i);  % ���������ʽ
        zeta(i)=0;  % disp('����δ�����߽�Լ��');  
    end
end
 y=sum(zeta);
end