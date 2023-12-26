%% ������: Solution of the aerodynamic force
clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ú�����ò��������������ϵ���ĺ���
% wing_para_output=[R_wingeff,C_avereff,F_nd,Y_rcpnd,Y_rnd,I_xxam,I_xyam,I3,I3y,I4y];
wing_para=wing_shape_fruitfly();      %���ú���wing_shape;  
% size(wing_para)
r2_nd=0.5801; 
R_wingeff=wing_para(1,1);  % R_wingeff =3.0040;        % ��λ�� mm
C_avereff=wing_para(1,2);  % C_avereff =0.8854;          % ��λ�� mm
F_nd=wing_para(1,3);         % F_nd = 0.4639;                % ���ٵ�λ��mm^4
Y_rnd=wing_para(1,4);        % Y_rnd =0.1402;               % ���ٵ�λ��mm
% I_xyam=wing_para(1,5);      % I_xyam =-0.0012;            % ��λ�� mg.mm^2;  
% I_xxam=wing_para(1,6);      % I_xxam =2.5080e-004;    % ��λ�� mg.mm^2;  
I3=wing_para(1,7);              % I3 =0.7485;                      % �����٣����ٻ���λ�� mm^4
I3y=wing_para(1,8);            % I3y = 0.0051;                    % ��λ�� mg.mm
I4y=wing_para(1,9);            % I4y =-5.3461e-004;          % ��λ�� mg.mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',dalpha',ddphi',ddpsi',ddalpha',C_L',C_D',C_N1',C_N2',C_N3',C_T',alpha2'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly();      %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)
size(wing_kenimatics)
t=wing_kenimatics(:,1);                % ��λ��ms
phi=wing_kenimatics(:,2);        % �Ĵ�ǡ�����λ��rad
psi=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
alpha1=wing_kenimatics(:,4);      % alpha1=pi/2+psi.*sign(dphi);    %���ι��ǡ���������             %���������ȫ�����ι���
alpha2=wing_kenimatics(:,5);      % alpha2=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi=wing_kenimatics(:,6);          % ��λ��rad/s
dpsi=wing_kenimatics(:,7);          % ��λ��rad/s�������ܻ�����⡪alpha��XXXXXX
% dalpha1=wing_kenimatics(:,8);      % ��λ��rad/s����%������������и�
ddphi=wing_kenimatics(:,9);       % ��λ��rad/s^2
ddpsi=wing_kenimatics(:,10);     % ��λ��rad/s^2�������ܻ�����⡪alpha��XXXXXX
% ddalpha1=wing_kenimatics(:,11);  % ��λ��rad/s^2����%������������и�
C_L=wing_kenimatics(:,12);          % ���ܻ�����⡪alpha��XXXXXX
C_D=wing_kenimatics(:,13);         % ���ܻ�����⡪alpha��XXXXXX
C_N1=wing_kenimatics(:,14);   % ���ܻ�����⡪alpha��XXXXXX
C_N2=wing_kenimatics(:,15);       % ���ܻ�����⡪alpha��XXXXXX
% C_N3=wing_kenimatics(:,16);   % ���ܻ�����⡪alpha��XXXXXX
% C_T=wing_kenimatics(:,17);     % ���ܻ�����⡪alpha��XXXXXX
alpha=alpha2;
% C_N=C_N2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
% �������ϵ�µĽ��ٶ�
omega_x=-dpsi; 
omega_y=dphi.*sin(psi); 
omega_z=dphi.*cos(psi);
omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% �������ϵ�µĽǼ��ٶȡ����������������ļ���
domega_x=-ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ������������������ϵ�£�ƽ�������������������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��4��
% C_avereff;  R_wingeff;    F_nd;                �����������������������Գ���ò������;% F_nd=0.50236;
% omega_h;  C_L(alpha);  C_D(alpha);       �����������������������Գ�2DOF�˶�����������, ����
%%
Rou=1.225;        %���ﵥλ��Kg/m^3�����ɻ���ɡ���10^6/(10^3)^3=10^-3mg/mm^3    
g=9.821/10^6;   % �����������ٶ�:g=9.821N/kg=9.821/10^6=9.821e-006   N/mg  ����g�Ĺ��ʵ�λ��m/s^2��N/kg
LD_aerocoeff=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);  % ��λ�ǣ�kg/m^3*mm*m^3= 10^(-12) kg.m
lift_wingframe=sign(alpha).*omega_h.^2.*C_L*LD_aerocoeff/g;      % ��λ��N / (N/mg) =mg
drag_wingframe=omega_h.^2.*C_D*LD_aerocoeff/g;
f=188.7;  T=1/f;          % Hz
Phi=max(phi);              % Phi*180/pi=73.4427
S=0.0266*100*10^(-6);      % m^2
U_ref=2*f*Phi*r2_nd*R_wingeff*10^(-3);   % m/s   % U_ref =0.8430
C_L=(lift_wingframe*g)/(Rou*S*U_ref^2);
figure(30)
plot(t/T,C_L,'r-');   
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��˲ʱ����������ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ��������
lift_aver=trapz(t,lift_wingframe)/T                                    % 
drag_aver=trapz(t,drag_wingframe)/T                             %���������Ҳ��ȥ�˾���ֵ�ģ�
drag_instabs_aver=trapz(t,abs(drag_wingframe))/T;        %���������Ҳ��ȥ�˾���ֵ�ģ�
lift2drag_aver=lift_aver/drag_instabs_aver                 %�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������:
% lift_aver =1.0913
% drag_aver =0.1573
% lift2drag_aver =1.1037
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13)
F_LD=plot(t/T,lift_wingframe,'r-',t/T,drag_wingframe,'k-');   
xlabel('\itNormalized time')
ylabel('\itAerodynamic force (mg)')                    %��λ���㵽mg��Ŷ
title('Aerodynamic lift and drag force \itvs. \itt \rm for flapping wing')
% legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\it|F_D|)')
legend('\itinstantaneous lift (\itF_L)','\itinstantaneous drag (\itF_D)')
set(F_LD,'LineWidth',2)
grid on
% axis([-inf,inf,-900,900])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƽ������������(1)���������������õ㣺ƬԪ�ҳ��е�
F_ny_coeff=(1/2)*Rou*C_avereff*R_wingeff^3*F_nd*10^(-12);        % ��λ��kg/m^3*m*m^3=kg.m
F_ny1=omega_h.^2.*C_N1*F_ny_coeff*10^6;  % ��λ��uN
F_ny2=omega_h.^2.*C_N2*F_ny_coeff*10^6;  % ��λ��uN�����������������
figure(14)
plot(t/T,F_ny1,'r-',t/T,F_ny2,'g-')
xlabel('\itNormalized time')
ylabel('\itF_{ny_1} and F_{ny_2} (uN)')
legend('F_{ny_1}','F_{ny_2}')
title('ƽ������������(����)��ʱ��ı仯����')   % ƽ������������(����)��ʱ��ı仯����
grid on
% ƽ������������/ƽ��������������������������ϵ��
% C_Nt=F_ny1./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %��λ��uN
% figure(15)                                                       
% plot(t/T,C_Nt,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt')
% legend('C_Nt')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƽ������������(2)�����������������������õ㣺ƬԪ�ҳ��е�
% % ����������������F_nd=R2nd2+2*xr_nd*R1nd1+xr_nd^2;   �ڳ�������������پ���xr_nd=0ʱ; F_nd=R2nd2
% R2_nd2=0.2350;      % ������������������F_nd=R2nd2�����������Ľ������׼ȷ
% F_ny_total=-1/2*Rou*R_wingeff^3*C_avereff*R2_nd2*C_N3.*dphi.*abs(dphi)*(10^-3*10^-9*10^3);                    %���򣬵�λ��mN
% F_Ttz=-1/2*Rou*R_wingeff^3*C_avereff*R2_nd2*C_T.*dphi.*abs(dphi).*psi.*abs(psi)*(10^-3*10^-9*10^3);  %���򣬵�λ��mN
% figure(16)
% F_LD2=plot(t/T,F_ny_total,'r-',t/T,F_Ttz,'k-');   
% xlabel('\itNormalized time')
% ylabel('\itAerodynamic force (mN)')
% title('Aerodynamic normal and tangential force \itvs. \itt \rm for flapping wing')
% legend('\itinstantaneous normal force (\itF_N)','\itinstantaneous tangential force (\itF_T)')
% set(F_LD2,'LineWidth',2)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ�������������������ϵ�£��������������������������õ㣺ƬԪ�ҳ��е�
%%% �Ĵ�ͱ���Ťת��Эͬ���ã������������������ο�SK Agrawal�������Ƶ����̣��Ա��������������
%% �������Ϊ��10��
% I3y;   I4y;                                                         �����������������������Գ���ò������
% domega_z; omega_x; omega_y;                     �����������������������Գ�2DOF�˶���Ǽ�����, ������
% dphi; dalpha; ddphi; ddalpha; alpha;            �����������������������Գ�2DOF�˶�������, �Ǽ�����, ����
%%  (1) 2010_JFM_Aeromechanics-��������������ʽ����% I3y��I4y��λ�� mg.mm-���㵽kg.m��Ҫ*10^(-9)
% F_yam=sign(phi).*(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x)*(10^(-9)*10^6);   %��λ��uN
% F_yam1=(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x)*(10^(-9)*10^6); % ԭ���Ƴ��Ĺ�ʽ
% ע�⹫ʽǰ(-)���������ݹ�Ӭ�������ϵ�޸��˽��ٶȺͽǼ��ٶ�
F_yam1=-(I3y*(domega_z-omega_x.*omega_y)+I4y*domega_x/4)*10^(-3);    %��λ��uN
%% (2)������2001-JEB-Sanjay SP-��������������ʽ
% F_yam2=(I3y*(ddphi.*sin(alpha)+dphi.*dalpha.*cos(alpha))-(I4y/4)*ddalpha/4)*10^(-3);% ԭ���Ƴ��Ĺ�ʽ��������
% �����alpha1Ϊȫ�����ι��ǣ������и���dalpha1��ddalpha1
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))-(I4y/4)*ddalpha1/4)*10^(-3); 
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dalpha1.*cos(alpha1))+(I4y/4)*ddalpha1/4)*10^(-3); 
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))-(I4y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% F_yam2=-(I3y*(ddphi.*sin(alpha1)+dphi.*dpsi.*cos(alpha1))+(I4y/4)*ddpsi/4)*10^(-3); %�����dpsi��ddpsi����仯��Ҫȷ��
% �����alpha2�����и�, dalpha2��ddalpha2��Ҫ��
% F_yam2=-(I3y*(ddphi.*sin(alpha2)+dphi.*dalpha2.*cos(alpha2))-(I4y/4)*ddalpha2/4)*10^(-3); 
i=length(alpha2);   % size(alpha2)=(2000*1)
dalpha2=[diff(alpha2);alpha2(i,1)];    
j=length(dalpha2);
ddalpha2=[diff(dalpha2);dalpha2(j,1)];  % size(ddalpha2)
F_yam2=-(I3y*(ddphi.*sin(abs(alpha2))+dphi.*dalpha2.*cos(abs(alpha2)))-(I4y/4)*ddalpha2/4)*10^(-3); 
figure(17)          % ����������ЩС
hold on
plot(t/T,F_yam1,'r-',t/T,F_yam2,'b-')                    % ע�������domega_x=ddpsi��ddalpha�ķ��Ŵ��ڲ�ͬŶ
xlabel('\itNormalized time')
ylabel('\itF_{am_y} (uN)')
legend('F_{am1_y}','F_{am2_y}')
title('������������(����)��ʱ��ı仯����')          % ������������(����)��ʱ��ı仯����
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-');     %��x-axis
grid on
%% ������������/ƽ��������������������������ϵ��
% C_Na=F_yam2./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %��λ��uN
% figure(18)                                                       
% plot(t/T,C_Na,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Na')
% legend('C_Na')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���������������� ������ϵ�£� ��ת�����������������������õ㣺ƬԪ�ҳ��е�
%% �������Ϊ��3��
% C_avereff;   R_wingeff;   I3;                              �����������������������Գ���ò������
% omega_h;  omega_x;                                       �����������������������Գ�2DOF�˶�����������, �������
%%
C_R=1.55;
% F_zrot=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*omega_x.*omega_h*I3*(10^(-12)*10^3);%��λ��mN����δ���Ƿ���
F_y_rotcoeff=(1/2)*Rou*R_wingeff^2*C_avereff^2*C_R*I3*10^(-12);          % ��λ��kg/m^3*mm^2*mm^2= kg.m 10^(-12)
% F_y_rot=sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6;    % ��λ��uN
F_y_rot=omega_x.*omega_h*F_y_rotcoeff*10^6;    % ��λ��uN
figure(19)                  % ����������ЩС  
plot(t/T,F_y_rot,'b-')
xlabel('\itNormalized time')
ylabel('\itF_{y,rot} (uN)')
legend('\itF_{y,rot}')
title('��ת����������(����)��ʱ��ı仯����')   % ��ת����������(����)��ʱ��ı仯����
grid on
%% ��ת����������/ƽ��������������������������ϵ��               %  ���ת������������ϵ�����岻���ֵ
% figure(20)                                                                       
% % C_Nr=F_y_rot./(-sign(dphi).*omega_h.^2*F_ny_coeff*(10^-6));  %��λ��uN
% C_Nr=F_y_rot./(sign(alpha2).*abs(omega_x).*abs(omega_h)*F_y_rotcoeff*10^6);  %��λ��mN
% plot(t/T,C_Nr,'r-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nr')
% legend('C_Nr')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ֻ��Ƶķ���������ϵ���Լ��ܺ�
% C_Ns=C_Nt+C_Na+C_Nr;  %��λ��mN
% figure(21)                                                                         %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Nt,'r-',t/T,C_Na,'g-',t/T,C_Nr,'b-',t/T,C_Ns,'k-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt & C_Na & C_Nr & C_Ns')
% legend('C_Nt','C_Na','C_Nr','C_Ns')
% grid on
%% ����������ϵ���ܺͷֽ⵽��ֱ����������(C_L)�Լ�ƽ������������(C_D)
% C_Lv=sign(dphi).*C_Ns.*sin(psi);     % ��ֱ����vertical������ϵ��
% C_Dh=sign(dphi).*C_Ns.*cos(psi);    % ˮƽ����horizontal������ϵ��
% figure(22)                                                                       %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Nt,'r-',t/T,C_Lv,'g-',t/T,C_Dh,'b-')
% xlabel('\itNormalized time')
% ylabel('\itC_Nt & C_Lv & C_Dh')
% legend('C_Nt','C_Lv','C_Dh')
% title('����������ϵ���ܺͷֽ⵽��ֱ����������(C_L)�Լ�ƽ������������(C_D)') 
% grid on
% figure(23)  % ����������ϵ���ܺͷֽ⵽��ֱ����������(C_L)                                                                       %   XXXXXXXXXXXXXXXXX
% plot(t/T,C_Lv,'g-')
% xlabel('\itNormalized time')
% ylabel('\itC_Lv')
% legend('C_Lv')
% title('����������ϵ���ܺͷֽ⵽��ֱ����������(C_L)')   
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������ ������ϵ�£���(�ϳ�)������������ƽ������������+��ת����������+�����������������õ㣺ƬԪ�ҳ��е�
F_N=F_ny1+F_y_rot+F_yam1;    % ����ƽ��������F_ny��ƽ��������ϵ���ϳɶ���,%��λ��mN
% F_N=F_ny_total+F_y_rot+F_yam;  % ����ƽ��������F_ny_total��ƽ������������ϵ���������,%��λ��mN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������ȶ������ڵ���(�ϳ�)�������������������Excel����mN
% �Ա�thr_dim_chord2����ĵ������ڻ������ͼ
B=F_N;    % size(B)  % (2000*1)     %����mN
xlswrite('Forcenormal_oneT.xlsx',B,'sheet1','A1:A2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(24)
plot(t/T,F_N,'r-')
xlabel('\itNormalized time')
ylabel('\itF_N (mN)')
legend('F_N')
title('ƽ������������+��ת����������+������������(����)��ʱ��ı仯����')
grid on
%% ������ϵ�£�ƬԪ��������ֽ⵽ƬԪ��ֱ�����ˮƽ�������õ㣺ƬԪ�ҳ��е�
F_vertical=abs(F_N.*cos(alpha));                          % ��ֱ����vertical, %��λ��mN
F_horizontal=sign(alpha).*F_N.*sin(alpha);   % ˮƽ����horizontal,%��λ��mN
% F_vertical=abs(F_N.*cos(alpha));                          % ��ֱ����vertical, %��λ��mN
% F_horizontal=abs(sign(alpha).*F_N.*sin(alpha));                                   % ˮƽ����horizontal,%��λ��mN
% F_v=(abs(F_N).*cos(alpha))*10^-3/g;     % ��ֱ����vertical,%��λ��mg
% F_h=(F_N.*sin(alpha))*10^-3/g;             % ˮƽ����horizontal,%��λ��mg
F_vaver=trapz(t,F_vertical)/T                       % T��������������Ŷ
F_haver=trapz(t,F_horizontal)/T                       %���������Ҳ��ȡ�˾���ֵ�ģ�
F_haverabs=trapz(t,abs(F_horizontal))/T;         %���������Ҳ��ȡ�˾���ֵ�ģ�
F_v2haver=F_vaver/F_haverabs          %�����
%%%output: Units: mg
% F_vaver =10.8146
% F_haver =1.5524
% F_v2haver =0.9314
figure(25)
plot(t/T,F_vertical,'r-',t/T,F_horizontal,'b-')
xlabel('\itNormalized time')
ylabel('\itF_{vertical} & F_{horizontal} (mN)')
legend('F_{vertical}','F_{Phorizontal}')
title('��������ֽ⵽��ֱ�����ˮƽ����ķ�������ʱ��ı仯����')   
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£� �ɷ���������ϵ�����: ��ֱ�����ˮƽ����������ϵ���������������ⲻ��ȷ
% C_N=F_N./((1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));   % ��λ������ŶXXXX
C_N=F_N./(-(1/2*Rou*omega_h.^2.*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));   % ��λ������ŶXXXX
C_v=F_vertical./(-sign(dphi).*(1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3)); 
C_h=F_horizontal./(-sign(dphi).*(1/2*Rou*sign(dphi).*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3));
% % C_h=F_horizontal./((1/2*Rou*omega_h.^2*C_avereff*R_wingeff^3*F_nd)*(10^-3*10^-9*10^3)); 
% C_v=C_N.*cos(alpha); 
% C_h=C_N.*sin(alpha); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��˲ʱ������ϵ��ʹ�����λ��ֺ���trapz������ֵ���֣���������, ���ƽ����ֱ�����ˮƽ����������ϵ��
C_vaver=trapz(t,C_v)/T                     
C_haver=trapz(t,C_h)/T                      %���������Ҳ��ȡ�˾���ֵ�ģ�
C_habsaver=trapz(t,abs(C_h))/T;        %���������Ҳ��ȡ�˾���ֵ�ģ�
C_v2haver=C_vaver/C_habsaver              %�����
% ���
% C_vaver = 3.7772e+005
% C_haver =4.5096e+006
% C_v2haver =0.0540
% figure(26)
% plot(t/T,C_N,'g-')
% xlabel('\itNormalized time')
% ylabel('\itC_N')
% legend('C_N')
% title('�����������ϵ����ʱ��ı仯����')
% grid on
% figure(27)
% plot(t/T,abs(C_v),'r-',t/T,C_h,'b-')
% xlabel('\itNormalized time')
% ylabel('\itC_v & C_h')
% legend('C_v','C_h')
% title('��ֱ�����ˮƽ����������ϵ����ʱ��ı仯����')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

