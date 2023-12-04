% function obj_function=objfunction_aeroF_stdev(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ԭʼ�ļ�������Aero_F3_fruitfly_exp
% objfunction_aeroF_stdev
% (1)�漰���׼��̬����Ԥ����õĳ�������ϵ����������ʵ����Ի�е��Ӭģ�ͻ�õ�����֮��ı�׼����
% (2) ���Կ��ǳ��������Ͳ����ǳ����������ֱ�����Ż�
% (3) �漰����ĸ�������ϵ���ĳͷ����Լ��;����F_N=F_ytran+F_yrot+F_yadd1+F_inert_y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
% x=[1,1,0.356737,1,1];
% x=[1,1,1,1,0];  % obj_function =40.1316;
% k_C_N=x(1);
% k_C_R=x(2);
% x0_nd=x(3);  % Ťת���λ�� x0_nd=0.356737;
% k_add=x(4);
% k_inert=x(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;
% % x=[1,1.281,1,1];
x=[1,1.8,0.35,1]; % ��x=[1,1.8,0.15,1]; Ҫ�õ�Ŷ
% x=[0.85,0.75,0.05,1]; 
% x =[1,0.5,0.571,0.01];  % obj_function =36.6581;
% x =[0.5,0.5,0.2000];  %fval =18.7559;
% x =[0.5, 0.5,1.0889]; % fval = 12.6197;
% x =[1,0.5,0.5,0.01];
% x =[1,0.5,0.5,0.5,0];
% x =[1.0000,1.0000,0.7500,0.5643,1.5000];  % fval =35.5047;  % Elapsed time is 583.787292 seconds. %F_Nû�п��ǹ�������k_inert��������1.5����ʵ����������
% x =[1.0000,1.6003,0.7500,0.5574,0.0100]; % fval =35.5099;
% x= [1.0000,0.5000,0.3634, 0.5000];  % fval =36.8167;
% x =[1.0000,0.5000,0.3358,0.5000];       % fval=36.5981; 
% x =[1.0000,1.6003,0.7500,1]; 
k_C_N=x(1);
k_C_R=x(2);
% x0_nd=x(3);  % Ťת���λ�� 
x0_nd=0.356737; % Ťת���λ�� 
k_add=x(3);
k_inert=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������: Solution of the aerodynamic force
% �޸�ʱ�䡪��2014��12��20��,11:54����ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ú�����ò��������������ϵ���ĺ���
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
wing_para=wing_shape_fruitfly_sixteen_good();   %���ú���wing_shape_fruitfly;  % size(wing_para)
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46392;  ������, ���ٻ���λΪmm^4
Coeff_liftdragF_N=wing_para(1,2);      % Coeff_liftdragF_N=0.00682;  %��λ��mg*mm
% M_xaercoeff=wing_para(1,3);          % M_xaercoeff=0.006038;   %��λ��: mg.mm^2
% I1z=wing_para(1,4);                         % I1y=0.016158   % ��λ�� mg.mm^2
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  ������, ���ٻ���λΪmm
% M_xrdcoeff=wing_para(1,6);           % M_xrdcoeff=0.0001839; % ��λ��mg.mm^2
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74851;  ������, ���ٻ���λΪmm^4
F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.0032433;  % ��λ�� mg.mm
% M_xRotcoeff=wing_para(1,9);         % M_xRotcoeff=0.002871;   % ��λ�� mg.mm^2
% I2z=wing_para(1,10);                      % I2y=0.006943;        % ��λ�� mg.mm^2
% I_xzam=wing_para(1,11);                % I_xzam=0.001424  % ��λ�� mg.mm^2
% I_xxam=wing_para(1,12);                % I_xxam=0.000338  % ��λ�� mg.mm^2
I5y=wing_para(1,13);                          % I5z=0.0050926   % ��λ�� mg.mm
I6y=wing_para(1,14);                          % I6z=0.00077164  % ��λ�� mg.mm
% I7z=wing_para(1,15);                      % I7y=0.0109056;        % ��λ�� mg.mm^2
% M_zrdcoeff=wing_para(1,16);         % M_zrdcoeff=0.0011169; % ��λ�� mg.mm % ת�������������ز������Ƴ�ƽ���µ�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % ��λ�� mg.mm        % ������I3yӦ�ø�ΪI5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % ��λ�� mg.mm         % ������I4yӦ�ø�ΪI6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % ��λ�� mg.mm^2    % ������I5zӦ�ø�ΪI7z
% I1z=wing_para(1,4);                         % I1y=0.0162;        % ��λ�� mg.mm^2    % ������I7zӦ�ø�ΪI1z
% I2z=wing_para(1,10);                       % I2y=0.0069;        % ��λ�� mg.mm^2    % ������I6zӦ�ø�ΪI2z 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(); %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,11)
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                % ��λ��ms
phi=wing_kenimatics(:,2);        % �Ĵ�ǡ�����λ��rad
psi=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
alpha=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi=wing_kenimatics(:,5);          % ��λ��rad/s
dpsi=wing_kenimatics(:,6);          % ��λ��rad/s
ddphi=wing_kenimatics(:,7);       % ��λ��rad/s^2
ddpsi=wing_kenimatics(:,8);     % ��λ��rad/s^2
C_L=wing_kenimatics(:,9);          
C_D=wing_kenimatics(:,10);     
C_N1=wing_kenimatics(:,11);   
C_T=wing_kenimatics(:,12);   
C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
% �������ϵ�µĽ��ٶ�
omega_x=dpsi;                     % չ��
omega_y=dphi.*sin(psi);       % ����(��ʼ����)
omega_z=dphi.*cos(psi);      % ������
omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% �������ϵ�µĽǼ��ٶȡ����������������ļ���
domega_x=ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% �������Ǻ͵��������ٶȵļ��㡪��ȡ�����ٶ�����ڸ�������ٶȣ���������Ӹ��ţ�
v_y_nonr=-omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;    % v_z=-r*dphi*sin(psi)
% alpha2=atan2(-v_y_nonr,v_z_nonr);   % ��ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ
% % ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % ���������ٶ�V_nonr=omega_h=dphi;   % ��λ�� rad/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xr=0.3289;  % x-root offset  \mm      % R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;   % ��λ�� mm
% R_ref=(xr+R_wingeff)*10^(-3);  
% R_ref_non=3.3293*10^(-3);
R_ref_non=1;   % r2_nd=0.5801;   
v_y=R_ref_non*omega_z;
v_z=-R_ref_non*omega_y;
V_ref=sqrt(v_y.^2+v_z.^2);
f=188.7;  T=1/f;          % Hz
V_ref_aver=trapz(t,V_ref)/(3*T);    % ���ο��ٶ�: V_ref_aver=867.0901rad/s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ������������������ϵ�£�ƽ�������������������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��4��
% C_avereff;  R_wingeff;    F_nd;                �����������������������Գ���ò������;% F_nd=0.50236;
% omega_h;  C_L(alpha);  C_D(alpha);       �����������������������Գ�2DOF�˶�����������, ����
%%  Coeff_liftdragF_N   %��λ  mg.mm=*10^(-9)kg.m
% g=9.821*10^(-6);   % �����������ٶ�:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
g=9.821;   % �����������ٶ�:g=9.821N/kg=9.821*10^(-6)=9.821e-006   N/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
% lift_inst=V_nonr.^2.*C_L*Coeff_liftdragF_N*10^(-3)/g;      % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
% drag_inst=V_nonr.^2.*C_D*Coeff_liftdragF_N*10^(-3)/g;  % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
lift_inst=V_nonr.^2.*C_L*Coeff_liftdragF_N*10^(-3);      % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
drag_inst=V_nonr.^2.*C_D*Coeff_liftdragF_N*10^(-3);  % ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��˲ʱ����������ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ��������
lift_aver=trapz(t,lift_inst)/(3*T);                            % lift_aver =1.0913 
drag_aver=trapz(t,drag_inst)/(3*T);                      % drag_aver =0.9887
drag_instabs_aver=trapz(t,abs(drag_inst))/(3*T);        
lift2drag_aver=lift_aver/drag_instabs_aver;         % ������������ֵ: lift2drag_aver =1.1037������������ϵ���ı�ֵ��û��ʵ�����壿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(10)
% F_LD=plot(t/T,lift_inst,'r-',t/T,drag_inst,'b:','LineWidth',2);   
% xlabel('\itNormalized time')
% ylabel('\itAerodynamic force (mg)')                    %��λ���㵽mg��Ŷ
% title('Aerodynamic lift and drag force \itvs. \itt \rm for flapping wing')
% % legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\it|F_D|)')
% legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\itF_D)')
% set(F_LD,'LineWidth',2)
% grid on
% axis([0.9,3.0,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƽ������������(1)���������������õ㣺ƬԪ�ҳ��е�
% C_N=C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %��������ϵ���ϳɡ���2010-JFM-RJ Wood
% Coeff_liftdragF_N=6.8201e-012������λ��:kg/m^3*mm*mm^3= 10^(-12) kg.m
% F_ytran=-sign(alpha).*V_nonr.^2.*C_N*Coeff_liftdragF_N/g;  %��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg������������ʾ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(11)
% plot(t/T,F_ytran,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} (mg)')
% legend('F_{y,tran}')
% title('ƽ������������(����)��ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ��λ��rad^2*s^-2*kg*m / (m*s^-2)=mg
F_ytran=-k_C_N*sign(alpha).*V_nonr.^2.*abs(C_N)*Coeff_liftdragF_N*10^(-3);   % ��λ��rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
F_ytran_aver=trapz(t,abs(F_ytran))/(3*T);   % F_ytran_aver =14.5887uN;
F_ztran=-sign(alpha).*V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
% F_ytran=V_nonr.^2.*C_N*Coeff_liftdragF_N*10^(-3); 
% F_ztran=V_nonr.^2.*C_T*Coeff_liftdragF_N*10^(-3);  
% figure(12)
% plot(t/T,F_ytran,'k-',t/T,F_ztran,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} & F_{z,tran} (uN)')
% legend('F_{y,tran}','F_{z,tran}')
% title('ƽ������������(���������)��ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_tran=-sign(alpha).*F_ytran.*cos(alpha);              % ��ֱ����vertical, %��λ��uN  %  
F_horiz_tran=F_ytran.*sin(abs(alpha));                           % ˮƽ����horizontal,%��λ��uN
% figure(13)                  % ����������ЩС  
% hold on
% plot(t/T,F_vert_tran,'k-',t/T,F_horiz_tran,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{vert,tran} &  F_{horiz,tran} (uN)')
% legend('\itF_{vert,tran}','\itF_{horiz,tran}')
% title('��ת�����������Ĵ�ֱ�����ˮƽ����ķ�����ʱ��ı仯����')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ������������� ������ϵ�£� ��ת�����������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��3��
% C_avereff;   R_wingeff;   F_ndRot;  �����������������������Գ���ò������
% omega_h;  omega_x;                    �����������������������Գ�2DOF�˶�����������, �������
% F_zrot=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*omega_x.*omega_h*F_ndRot*(10^(-12)*10^3);%��λ��mN����δ���Ƿ���
% F_y_rot=C_R*sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6;    % ��λ��uN
% x0_nd=0.356737;
C_R=pi*(0.75-x0_nd);
% C_R=1;  
% C_R=1.55;
% F_yrot=-C_R*sign(alpha).*omega_x.*V_nonr*F_yrotcoeff*10^(-3); % 
F_yrot=k_C_R*C_R*omega_x.*V_nonr*F_yrotcoeff*10^(-3);      % ��λ��rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
F_yrot_aver=trapz(t,abs(F_yrot))/(3*T);                             % F_yrot_aver =2.8944uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(14)                  % ����������ЩС  
% hold on
% plot(t/T,F_ytran,'k-',t/T,F_yrot,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{y,tran} &  F_{y,rot} (uN)')
% legend('\itF_{y,tran}','\itF_{y,rot}')
% title('ƽ������������(����)����ת����������(����)��ʱ��ı仯����')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_rot=-sign(alpha).*F_yrot.*cos(alpha);              % ��ֱ����vertical, %��λ��uN  %  
F_horiz_rot=F_yrot.*sin(abs(alpha));                          % ˮƽ����horizontal,%��λ��uN
% figure(15)                  % ����������ЩС  
% hold on
% plot(t/T,F_vert_rot,'k-',t/T,F_horiz_rot,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{vert,rot} &  F_{horiz,rot} (uN)')
% legend('\itF_{vert,rot}','\itF_{horiz,rot}')
% title('��ת�����������Ĵ�ֱ�����ˮƽ����ķ�����ʱ��ı仯����')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����������������������ϵ�£��������������������������õ㣺ƬԪ�ҳ��е�
%%% �Ĵ�ͱ���Ťת��Эͬ���ã������������������ο�SK Agrawal�������Ƶ����̣��Ա��������������
%% �������Ϊ��10��
% I5y;   I6y;     % ��λ�� mg.mm=10^(-9) kg.m �����������������������Գ���ò������
% domega_z; omega_x; omega_y;                      �����������������������Գ�2DOF�˶���Ǽ�����, ������
% dphi; dalpha; ddphi; ddalpha; alpha;               �����������������������Գ�2DOF�˶�������, �Ǽ�����, ����
%%  (1) 2010_JFM_Aeromechanics-��������������ʽ����% I5y��I6y��λ�� mg.mm=*10^(-9) kg.m
% F_yadd=sign(phi).*(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*(10^(-9)*10^6);   %��λ��uN
% F_yadd1=(I5y*(domega_z-omega_x.*omega_y)+I6y*domega_x)*10^(-3); % ԭ���Ƴ��Ĺ�ʽ������Ϊ����ЩŶ
% ע�⹫ʽǰ(-)���������ݹ�Ӭ�������ϵ�޸��˽��ٶȺͽǼ��ٶ�
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)-I6y*domega_x/4)*10^(-3); %��λ��rad^2*s^-2*kg*m=N=10^6uN % ��2001-JEBƥ��
% F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3); % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4  % ������һ����
% F_yadd1=k_add*(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);    % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4  % �����ڶ�����
F_yadd1=k_add*(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);  
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
% F_yadd1=k_add*(I5y*(ddphi.*cos(psi)+dphi.*dpsi.*sin(psi))+I6y*domega_x)*10^(-3);  
F_yadd1_aver=trapz(t,abs(F_yadd1))/(3*T);       % F_yadd1_aver=6.0929uN(old);  F_yadd1_aver=4.5472;
% (2)������2001-JEB-Sanjay SP-��������������ʽ
% F_yadd2=(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3);% ԭ���Ƴ��Ĺ�ʽ��������
% �����alpha1Ϊȫ�����ι��ǣ������и���dalpha1��ddalpha1
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))-(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))+(I6y/4)*ddalpha1/4)*10^(-3); 
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))-(I6y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% F_yadd2=-(I5y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))+(I6y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% �����alpha2�����и�, dalpha��ddalpha��Ҫ��
% F_yadd2=-(I5y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I6y/4)*ddalpha/4)*10^(-3); 
i=length(alpha);   % size(alpha)=(2000*1)
dalpha=[diff(alpha);alpha(i,1)];    
j=length(dalpha);
ddalpha=[diff(dalpha);dalpha(j,1)];  % size(ddalpha)
% F_yadd2=-(I5y*(ddphi.*sin(abs(alpha))+dphi.*dalpha.*cos(abs(alpha)))-(I6y/4)*ddalpha/4)*10^(-3);% ������һ����
F_yadd2=(I5y*(ddphi.*sin(abs(alpha))+dphi.*dalpha.*cos(abs(alpha)))-(I6y/4)*ddalpha/4)*10^(-3);
%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(16)  % ��λ�� mg.mm*rad^2.*s^-2=10^(-9) kg.m*s^-2=10^(-9)N=10^(6)uN;   %��λ��uN
% F_yadd1_tran=-(I5y*(domega_z+omega_x.*omega_y))*10^(-3);  % ����������������ƽ�������������
% F_yadd1_rot=(I6y*domega_x)*10^(-3);                                         % ����������������ת�������������
% plot(t/T,F_yadd1_tran,'r-.',t/T,F_yadd1_rot,'g:',t/T,F_yadd1,'k-','LineWidth',2) 
% ylabel('��������������F_{norm}(t) (uN)'); 
% legend('F_{y,add1,tran}(t)','F_{y,add1,rot}(t)','F_{y,add1}(t)');  
% title('��������������������F_{y,add1}(t))��ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%
F_vert_add=-sign(alpha).*F_yadd1.*cos(alpha);              % ��ֱ����vertical, %��λ��uN  %  
F_horiz_add=F_yadd1.*sin(abs(alpha));                          % ˮƽ����horizontal,%��λ��uN
% figure(17) 
% plot(t/T,F_vert_add,'r-.',t/T,F_horiz_add,'g:','LineWidth',2) 
% ylabel('��������������F_{norm}(t) (uN)'); 
% legend('F_{vert,add}(t)','F_{horiz,add}(t)');  
% title('�����������������Ĵ�ֱ�����ˮƽ����ķ�����ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(18)
% hold on
% plot(t/T,F_yadd1,'r-',t/T,F_yadd2,'b-','LineWidth',2)    % ע�������domega_x=ddpsi��ddalpha�ķ��Ŵ��ڲ�ͬŶ
% xlabel('\itNormalized time')
% ylabel('\itF_{y,add} (uN)')
% legend('F_{y,add1}','F_{y,add2}')
% title('������������(����)��ʱ��ı仯����')  
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����������������
% % �������ϵ�µĽ��ٶ�
% omega_x=dpsi;                     % չ��
% omega_y=dphi.*sin(psi);       % ����(��ʼ����)
% omega_z=dphi.*cos(psi);      % ������
% omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% % �������ϵ�µĽǼ��ٶȡ����������������ļ���
% domega_x=ddpsi;
% % domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_wing=2.4*10^-9;                    % kg
% x_com=1.920243385*10^-3;      %  m                           % ���ĵ�չ������
% % z_com=-0.149785466+0.636=0.486215*10^-3;      % ��Ťת����������
% z_com=0.149785466*10^-3;       %  m                          % ��Ťת����������
% F_inert_y=-k_inert*m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % ������������������
% % F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6; 
% F_inert_z=-k_inert*m_wing*(-domega_y*x_com-omega_y.^2*z_com+omega_x.*(omega_z*x_com-omega_x*z_com))*10^6; % ������������������
% % F_inert_y=-m_wing*(ddphi.*cos(psi)*x_com-(ddpsi-dphi.^2.*sin(2*psi)/2)*z_com)*10^6;       % ������������������
% % F_inert_z=-m_wing*(-ddphi.*sin(psi)*x_com-(dpsi.^2+dphi.^2.*(sin(psi)).^2)*z_com)*10^6;  % ������������������
% figure(19) 
% plot(t/T,F_inert_y,'r-',t/T,F_inert_z,'g-','LineWidth',2) 
% ylabel('��������������������F_{inert,y}(t) & F_{inert,z}(t) (uN)'); 
% legend('F_{inert,y}(t)','F_{inert,z}(t)');  
% title('��������������������F_{inert,y}(t) & F_{inert,z}(t)��ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_vert_inert=-sign(alpha).*F_inert_y.*cos(alpha);              % ��ֱ����vertical, %��λ��uN  %  
% F_horiz_inert=F_inert_y.*sin(abs(alpha));                          % ˮƽ����horizontal,%��λ��uN
% figure(20) 
% plot(t/T,F_vert_inert,'r-',t/T,F_horiz_inert,'g-','LineWidth',2) 
% ylabel('�������������������Ĵ�ֱ����F_{vert,inert}(t)��ˮƽ����F_{horiz,inert}(t)�ķ���  (uN)'); 
% legend('F_{vert,inert}(t)','F_{horiz,inert}(t)');  
% title('�������������������Ĵ�ֱ�����ˮƽ����ķ�����ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������ ������ϵ�£���(�ϳ�)������������ƽ������������+��ת����������+�����������������õ㣺ƬԪ�ҳ��е�
F_N=F_ytran+F_yrot+F_yadd1;    % ��λ��rad^2*s^-2*kg*m=N=10^6uN  ���������ǳ�������
% F_N=F_ytran+F_yrot+F_yadd1+F_inert_y;    % ��λ��rad^2*s^-2*kg*m=N=10^6uN�������ǳ�������
% figure(21)
% plot(t/T,F_N,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_N (uN)')
% legend('F_N')
% title('ƽ������������+��ת����������+������������(����)��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƬԪ��������ֽ⵽ƬԪ��ֱ�����ˮƽ�������õ㣺ƬԪ�ҳ��е�
F_vertical=-sign(alpha).*cos(alpha).*F_N;                % ��ֱ����vertical, %��λ��uN  %����ʵ���ϸ�ʽ������죬������ʵ�����Ա��Է���
% F_vertical=sign(alpha).*sin(psi).*F_N;                   % ��ֱ����vertical, %��λ��uN  % �������ݹ�ʽ�Ƶ��ò��ñ�ʽ
F_horizontal=sin(abs(alpha)).*F_N;                          % ˮƽ����horizontal,%��λ��uN
% % F_horizontal=-sign(alpha).*F_N.*sin(alpha);     % ˮƽ����horizontal,%��λ��uN
% F_vertical=abs(F_N.*cos(alpha));                          % ��ֱ����vertical, %��λ��uN
% F_horizontal=abs(sign(alpha).*F_N.*sin(alpha));  % ˮƽ����horizontal,%��λ��uN
%% ��������ϵ�£� �����(side force)F_horiz_X��ˮƽ������(drag)=F_horiz_y
% F_Z=F_vertical;
% F_horiz_X=-sin(phi).*F_horizontal;
% F_horiz_Y=cos(phi).*F_horizontal;  
F_horiz_X=-sin(phi).*cos(psi).*F_N; % �����(side force)F_horiz_X
F_horiz_Y=cos(phi).*cos(psi).*F_N;  % ˮƽ������(drag)=F_horiz_y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_uN=9.821;   % �����������ٶ�:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
F_vaver=trapz(t,F_vertical)/(3*T)/g_uN;               % F_vaver =1.2532; %��λ��mg      1.2532(��������)
F_haver=trapz(t,F_horiz_Y)/(3*T)/g_uN;          % F_haver =-0.1109; %��λ��mg     -0.1109(��������)
F_haverabs=trapz(t,abs(F_horiz_Y))/(3*T);    
F_v2haver=F_vaver/F_haverabs;                % ������������ֵ: F_v2haver =0.1150;  0.1045(��������)������������ϵ���ı�ֵ��û��ʵ�����壿
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;         % �����������ٶ�:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
M_body =1.8;  %��Ӭ���������(mg)��Science����M_body =1.8e-06;(kg)
W=M_body*g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����science��е��Ӭ��������������ݡ�����һ��֮�������
% %%����Ķ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % output_1=[t_NOfreq,Fx_norm_stroke_butterfilt_mean,Fy_norm_stroke_butterfilt_mean,Fz_norm_stroke_butterfilt_mean];
% force_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2\force_science.xlsx','A1:D1145'); % ��������
%�������������ȷ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output_1=[t_NOfreq,Fx_norm_all_butterfilt_steady,Fy_norm_all_butterfilt_steady,Fz_norm_all_butterfilt_steady];
% force_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\wing_parameter\datanalysis_science_fruitfly\robotForcesTorques\ForceModulations\aeroforce_for_steady_wingbeat.xlsx','A1:D1145');
force_science=xlsread('aeroforce_for_steady_wingbeat.xlsx',1,'A1:D1145'); % ��������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_NOfreq1=force_science(:,1);  % t_NOfreq(1,1)=0;
Fx_norm1=-force_science(:,2);
Fy_norm1=force_science(:,3);
Fz_norm1=-force_science(:,4); 
t_NOfreq=[t_NOfreq1+0.0052824335;t_NOfreq1+T+0.0052824335;t_NOfreq1+2*T+0.0052824335];
Fx_norm=[Fx_norm1;Fx_norm1;Fx_norm1];
Fy_norm=[Fy_norm1;Fy_norm1;Fy_norm1];
Fz_norm=[Fz_norm1;Fz_norm1;Fz_norm1];       % size(Fz_norm)  % (3435*1)
% figure(22)
% plot(t_NOfreq/T,Fz_norm,'r--',t_NOfreq/T,Fx_norm,'b--','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itF_{z,exp} & F_{x,exp} (uN)')
% legend('F_{z,exp}','F_{x,exp}')
% title('����������ʵ���õĴ�ֱ����������ˮƽ�������������2014-Science-MH Dickinson')   
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��������ֽ⵽��ֱ�����ˮƽ����ķ�������ʱ��ı仯������ʵ����Խ���ĶԱ�
figure(44)
plot(t_NOfreq/T,Fz_norm,'r-',t_NOfreq/T,Fx_norm,'b-','LineWidth',2.5)  % ʵ����Խ��
hold on
% 1.525�������ǳ�������  % 1.65���������ǳ��������������������ϸ������Ͻ�Ӧ������������������������
% plot(t/T,2*F_vertical/W,'r-.',t/T,2*F_horiz_Y/W,'b-.','LineWidth',2)       % ����Ԥ����
% plot(t/T,1.65*F_vertical/W,'r-.',t/T,1.65*F_horiz_Y/W,'b-.','LineWidth',2)  
% plot(t/T,1.5*F_vertical/W,'r-.',t/T,2.4*F_horiz_Y/W,'b-.','LineWidth',2)    % ������ʵ�����Ա�Ŷ
% plot(t/T,1.45*F_vertical/W,'r-.',t/T,2.3*F_horiz_Y/W,'b-.','LineWidth',2.5)  
% plot(t/T,1.65*F_vertical/W,'r-.',t/T,2.3*F_horiz_Y/W,'b-.','LineWidth',2.5)  % �������� % ������ʵ�����Ա�Ŷ
plot(t/T,1.45*F_vertical/W,'r-.',t/T,2.0*F_horiz_Y/W,'b-.','LineWidth',2.5)     % ���������� % ������ʵ�����Ա�Ŷ x=[1,1.8,0.35,1];
% plot(t/T,2*F_vertical/W,'r-.',t/T,2.0*F_horiz_Y/W,'b-.','LineWidth',2.5) 
% plot(t/T,1.425*F_vertical/W,'r-.',t/T,2.225*F_horiz_Y/W,'b-.','LineWidth',2.5)  
xlabel('\rmNormalized time','Fontsize',24,'FontName','Times','FontWeight','Bold')
% ylabel('\itF_{vertical} & F_{horiz,Y} (uN)')
% legend('F_{vertical}','F_{horiz,Y}')
ylabel('\rmF / m_{body}g','Fontsize',24,'FontName','Times','FontWeight','Bold')
legend('F_{vert,z,exp}','F_{horiz,y,exp}','F_{vert,z,cal}','F_{horiz,y,cal}')
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',20,'FontName','Times','FontWeight','Bold')
% title('��������ֽ⵽��ֱ�����ˮƽ����ķ�������ʱ��ı仯������ʵ����Խ���ĶԱ�')   
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%
axis([0.9,4.05,-4,4])
% set(gca,'XTick',(0.9:0.1:4.05))
hold on
L=length(t);
plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis
axis([min(t/T),min(t/T)+1,-4,3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% figure(45)
% plot(t_NOfreq/T,Fz_norm,'r-',t_NOfreq/T,Fx_norm,'b-',t_NOfreq/T,Fy_norm,'g-','LineWidth',2)  % ʵ����Խ��
% hold on
% % 1.525�������ǳ�������  % 1.65���������ǳ��������������������ϸ������Ͻ�Ӧ������������������������
% plot(t/T,2*F_vertical/W,'r-.',t/T,2*F_horiz_Y/W,'b-.',t/T,0*F_horiz_X/W,'g-.','LineWidth',2)       % ����Ԥ����
% xlabel('\rmNormalized time','Fontsize',20,'FontName','Times','FontWeight','Bold')
% ylabel('\rmF / m_{body}g','Fontsize',20,'FontName','Times','FontWeight','Bold')
% legend('F_{vert,z,exp}','F_{horiz,y,exp}','F_{lateral,x,exp}','F_{vert,z,cal}','F_{horiz,y,cal}','F_{lateral,x,cal}')
% set(gca,'FontSize',12,'FontName','Times','FontWeight','Bold')
% axis([0.9,4.05,-4,4])
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ������ϵ�£� �ɷ���������ϵ�����: ��ֱ�����ˮƽ����������ϵ��
% % Coeff_liftdragF_N=0.0068201������λ��:mg*mm^-3*mm*mm^3= 10^(-9) kg*m
% % Coeff_liftdragF_N=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);  % kg*m
% C_N_total=F_N*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));          % ��λ: N/(rad^2*s^-2*kg*m)=һ������
% % C_v=F_vertical*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9)); 
% % C_h=F_horizontal*10^(-6)./(V_ref_aver.^2*Coeff_liftdragF_N*10^(-9));    % *10^3); 
% C_v=-sign(alpha).*C_N_total.*cos(alpha); 
% C_h=C_N_total.*sin(abs(alpha)); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %��˲ʱ������ϵ��ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ����ֱ�����ˮƽ����������ϵ��
% C_vaver=trapz(t,C_v)/(3*T)                 % C_vaver =2.3852;       2.4003(��������)   
% C_haver=trapz(t,C_h)/(3*T)                 % C_haver =-0.2655;   -0.2125(��������)  
% C_habsaver=trapz(t,abs(C_h))/(3*T);     
% C_v2haver=C_vaver/C_habsaver           % ��ֱ����ϵ����ˮƽ����ϵ���ı�ֵ: C_v2haver = 1.1298;  1.026(��������)
% figure(24)
% plot(t/T,C_N_total,'g-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itC_{N,total}')
% legend('C_{N,total}')
% title('�����������ϵ����ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% figure(25)
% plot(t/T,C_v,'r-',t/T,C_h,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itC_v & C_h')
% legend('C_v','C_h')
% title('��ֱ�����ˮƽ����������ϵ����ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ƽ��ֵ�ĶԱ�
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_v_norm_aver=trapz(t,1.5*F_vertical/W)/(3*T)
% Fz_norm_aver=trapz(t,Fz_norm)/(3*T)
% F_haver_norm=trapz(t,2.4*F_horiz_Y/W)/(3*T)
% Fx_norm_aver=trapz(t,Fx_norm)/(3*T)
% F_vert_relaerror=abs(F_v_norm_aver-Fz_norm_aver)/Fz_norm_aver*100
% F_horiz_relaerror=abs(abs(F_haver_norm)-abs(Fx_norm_aver))/abs(Fx_norm_aver)*100
% F_horiz_relaerror=abs(abs(F_haver_norm)-abs(Fx_norm_aver))/abs(F_haver_norm)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ĸ�������ϵ���ĳͷ����Լ��
coeff_con=aeroF_coeff_constraint(x);
% ������������ϵ���ĳͷ����penaltyfun
s=2000;
penaltyfun=s*coeff_con;           % penaltyfun =;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stdev=std((2*F_vertical-Fz_norm),0,1);   % ����S1(0)������(1)Ԫ�صı�׼����
stdev_v=std((abs(2*F_vertical)-abs(Fz_norm)),0,1);   % ����S1(0)������(1)Ԫ�صı�׼����
stdev_h=std((abs(2*F_horiz_Y)-abs(Fx_norm)),0,1);   % ����S1(0)������(1)Ԫ�صı�׼����
obj_function=stdev_v+stdev_h+penaltyfun;  % obj_function =40.1316;
% toc  % Elapsed time is 3.401964 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
