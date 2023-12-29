function I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis)   %���ú���inertia_moment;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% moment_inertia_scaling_law
R_wingeff=3.004;     % mm    % R_ratio =33.2889;
C_avereff=0.8854;    % mm    % C_ratio =38.0619;
R_ratio=R_wing/R_wingeff;
C_ratio=C_aver/C_avereff; 
%% �����òѧ�����仯�ǵ�����λ������
% ��ʼģ�ͳ����������ꡪ����������������Ťת��Ϊ����ͳ�������ʱ��������λ������
% �������λ������: mm
% XC=1.920243385;      % X=1.920243385;          % չ    ��
% YC=0.000000068;      % Y=0.000000068;          % ��ȷ��򡪡���������
% ZC=-0.149785466;    %  Z=-0.149785466         % �� ��  % ��Ťת�����Ծ���: ����ʵ�ʼ���
xr1=0.3289;  % ����Ҳ��򡪡�Z_0��ƫ��������(�ڳ�����)
XC=1.920243385-xr1; % 1.5914
% C_maxy=0.138924474377504; % ��Գ���wing_model_88_yaxis��: ��122��; C_maxy =0.1389; 
% YC=-0.149785466-C_maxy;    % YC=-0.288709941;
YC=-0.149785466;
% �����òѧ�����仯ʱ������������
XC_tran=R_ratio*XC+xr0;           % XC_tran =53.3030@R_wing=100; & xr1=0.3289; % mm;
YC_tran=C_ratio*YC;                    % YC_tran =-5.7011@C_aver=33.7; % mm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����òѧ�����仯ʱ�����Ĵ��ġ������Ծ�
% ��ʼģ�ͳ��ġ������Ĺ��Ծ�
% Mass_Properties
% M=0.000002251;         % ���� % unit: g
% ���Ծأ����ģ�         % g.mm^2
% Ixc=0.000000217;        % g.mm^2
% Iyc=0.000001362;    % g.mm^2
% Izc=0.000001145;        % g.mm^2 
% % ���Ի������ģ�     % g.mm^2
% Iyzc=-0.000000000;     
% Ixzc=0.000000053; 
% Ixyc=-0.000000000; 
Ixc_small=0.000000215;     % g.mm^2
Izc_small=0.000001129;      % g.mm^2
% �����òѧ�����仯ʱ�ġ������Ĺ��Ծ�
Ixc_big=R_ratio*C_ratio^3*Ixc_small;    % Ixc_big =0.3946@R_wing=100; % mm;   % g.mm^2
Izc_big=R_ratio^3*C_ratio*Izc_small;    % Izc_big =1.5852@C_aver=33.7; % mm;   % g.mm^2
%% �����òѧ�����仯�仯ǰ�󡪡�����
M_small=0.000002237;  % g
M_big=R_ratio*C_ratio*M_small;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����òѧ�����仯ʱ��������Ĵ����Ťת��ġ������Ծء���ƽ���ᶨ��
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xr_initial=0.3289;                        % x-root offset  \mm                    
% x_com_zaxis=XC-xr_initial+xr0;   % ����������������������������������������������
x_com_zaxis=XC_tran;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZC_tran=YC_tran;
% x_0��[0,0.5];  % Ťת���ǰԵ�������پ���
% C_maxylb=0.464385778290230;    % C_maxylb =0.4644;
% C_maxy25=0.138924474377504;   % C_maxy25 =0.1389; 
% C_maxyub=-0.186536829535222; % C_maxyub =-0.1865
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ratio_leadmax=C_lead_ymax/C_max_LtoT; %Ratio_leadmax=0.4644/1.3018=0.356737; %��Գ���ͳ�����ߵ�Ťת�ᶨ����ǰԵ
Ratio_leadmax=0.356737;  % ������
% C_lead_ymax=max(f_x_lead);   % ���: C_lead_ymax=0.4644; k_leadmax=644;
% C_lead_ymax=0.464385778290230;               % mm���������òѧ�����仯ʱ���ǰԵ��z����
wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis);
C_max_LtoT=wing_para(1,17);
z_maxlead=Ratio_leadmax*C_max_LtoT;            % mm���������òѧ�����仯ʱ���ǰԵ��z����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_maxlead_xaxis=C_maxyaxis*C_max_LtoT;        % mm; Ťת���ǰԵ�����پ���
delta=z_maxlead-z_maxlead_xaxis;  % ����ͳ�������ᵽ��Ťת��֮������پ���
z_com_zaxis=abs(abs(ZC_tran)+delta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���Ծ� (WCS)   g.mm^2
% Ix=Ixc+M*ZC^2         % ��ƽ���ᶨ����
% Iy=Iyc+M*XC^2         % ��ƽ���ᶨ����
% Iz=Izc+M*XC^2         % ��ƽ���ᶨ����
Ix_inertia=Ixc_big+M_big*z_com_zaxis^2;     % ��ƽ���ᶨ����       % g.mm^2
Iz_inertia=Izc_big+M_big*x_com_zaxis^2;     % ��ƽ���ᶨ����       % g.mm^2
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];  % XXX
I_moment=[x_com_zaxis,z_com_zaxis,M_big,Ix_inertia,Iz_inertia];  % ��������Ϊ����Щ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



