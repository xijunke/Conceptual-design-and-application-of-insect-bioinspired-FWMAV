function F_M=Aero_F_M_fruitfly(x)
% function F_M=Aero_F_M_fruitfly()
% R_wing=3.004;
% C_aver=0.8854;
% xr0=0.3289;
% C_maxyaxis=0.25;
% R_wing=x(1);
% C_aver=x(2);
% xr0=x(3);
% C_maxyaxis=x(4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=x(1);    T=1/f;          % Hz
% phi_m=x(2);
% epsilon=x(3);
% psi_m=x(4);
% zeta=x(5);
% psi_0=x(6);
% %%%%%%%%%%%%%%%%%%
R_wing=x(7);
C_aver=x(8);
xr0=x(9);
C_maxyaxis=x(10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing_para=wing_shape_fruitfly_sixteen_good2(R_wing,C_aver,xr0,C_maxyaxis);% ע����õĺ�����wing_shape_fruitfly_sixteen_good2
% (1) ƽ�����������������ز���
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46391;  ������, ���ٻ���λΪmm^4
Coeff_liftdragF_N=wing_para(1,2);  % Coeff_liftdragF_N=0.00682;  %��λ��mg*mm  %ƽ������������
M_xaercoeff=wing_para(1,3);              % M_xaercoeff=0.006038;   %��λ��: mg.mm^2   % ƽ�������������ز��������Ƴ�ƽ���µ�չ����
% C_aver1=M_xaercoeff/Coeff_liftdragF_N;  % C_aver1=0.8854;
I1z=wing_para(1,4);                             % I1y=0.016158   % ��λ�� mg.mm^2            % ƽ�������������ز��������Ƴ�ƽ���µ�������
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  ������, ���ٻ���λΪmm
M_xrdcoeff=wing_para(1,6);                % M_xrdcoeff=0.0001839; % ��λ��mg.mm^2 %ת�������������ز������Ƴ�ƽ���µ�չ����
% (2) ת�����������������ز���
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74847;  ������, ���ٻ���λΪmm^4
F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.003243;  % ��λ�� mg.mm   % ת������������
M_xRotcoeff=wing_para(1,9);             % M_xRotcoeff=0.002871;   % ��λ�� mg.mm^2 % ת��������������ϵ�������Ƴ�ƽ���µ�չ����
% C_aver2=M_xRotcoeff/F_yrotcoeff    % C_aver2=0.8854;
I2z=wing_para(1,10);                          % I2y=0.006943;        % ��λ�� mg.mm^2            % ת�������������ز��������Ƴ�ƽ���µ�������
% (3) �����������������ز���
I_xzam=wing_para(1,11);                    % I_xzam =0.001424  % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
I_xxam=wing_para(1,12);                    % I_xxam =0.000338  % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
I5y=wing_para(1,13);                          % I5z=0.0050926   % ��λ�� mg.mm       % ������������������������
I6y=wing_para(1,14);                          % I6z=0.00077164    % ��λ�� mg.mm       % ������������������������
I7z=wing_para(1,15);                          % I7y=0.0109056;      % ��λ�� mg.mm^2   % �������������ز��������Ƴ�ƽ���µ�������
M_zrdcoeff=wing_para(1,16);             % M_zrdcoeff=0.001169; % ��λ�� mg.mm^2 % ת�������������ز������Ƴ�ƽ���µ�������
% C_max_LtoT=wing_para(1,17);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
%%%%%%%%%%%%%%%%%%%%%%%%%%
f=x(1);  T=1/f;   % Hz
phi_m=x(2);
K=x(3);                  %  epsilon=x(3);    
eta_m=x(4);          % phi_0=x(4);       
C_eta=x(5);           %  psi_m=x(5);
Phi_eta=x(6);        % zeta=x(6); 
% eta_0=x(7);       % psi_0=x(7);       
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,K,eta_m,C_eta,Phi_eta); %  6���������˶�ѧ
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,epsilon,psi_m,zeta,psi_0);  %  6���������˶�ѧ
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_sim();
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
% f=188.7;  % f=234;
% T=1/f;
t=wing_kenimatics(:,1);                % ��λ��ms
phi=wing_kenimatics(:,2);        % �Ĵ�ǡ�����λ��rad
psi=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
alpha2=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi=wing_kenimatics(:,5);          % ��λ��rad/s
dpsi=wing_kenimatics(:,6);          % ��λ��rad/s
ddphi=wing_kenimatics(:,7);       % ��λ��rad/s^2
ddpsi=wing_kenimatics(:,8);     % ��λ��rad/s^2
C_N=wing_kenimatics(:,9);   
% C_L=wing_kenimatics(:,10);          
% C_D=wing_kenimatics(:,11);   
% C_T=wing_kenimatics(:,12);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% (9) 2014-JRSI-Nabawy�������ǳ�򼸺���ò����Ӱ���������ϵ��
% % C_La2d=2*pi;    % ��λ��rad^(-1), ����2*pi (rad^(-1))= 0.11 (deg^(-1));   
% C_La2d=5.15662;  % ��λ��rad^(-1), ��Ե���������ŵ����ƽ��, ����ʹ�� 0.09 (deg^(-1)) % ����ǰԵ�������ⷨ�����Ĳ�����ɭ����ϵ��ģ��C_L_polhamus1��ȫ�غ�
% % The parameter E is the edge correction proposed by Jones [42] for the lifting line theory 
% % and is evaluated as the quotient of the wing semi-perimeter to its length[42,43].
% % lambda=0.755; % ��Ϊ����һ����������, �Ա���ʵ����Զ��õ�������ϵ�����жԱ�, �� C_La2d=2*pi;ʱ��ʵ����Ͻ���Ͻ�
% % lambda=0.61; % ��Ϊ����һ����������, �Ա���ʵ����Զ��õ�������ϵ�����жԱ�, �� C_La2d=5.15662; ;ʱ��ʵ����Ͻ���Ͻ�
% lambda=0.70; % ��������������������Dickinson��ϵ�ϵ��������õ����رȺ����������ܶȽϽӽ�,��P_asterisk=29.7503;  L =1.2761;
% C_periEdge=Cfcoeff(R_wing,C_aver);  % ���ú��������İ��ܳ�
% E=lambda*C_periEdge;   % E=C/2/R; ��Գ�ʼ��Ӭ��� E_C_peri=1.146071666; 
% % The parameter k is the so-called 'k-factor' included to correct for the difference in efficiency between the assumed ideal
% % uniform downwash distribution and real downwash distribution [35,40,44]. For this work, the k-factor required within
% % equation (2.11) will be estimated using the induced power factor expression of hovering actuator disc models [45].
% k=1.51;  % ��Թ�Ӭ
% % for the span of a single  wing 
% % AR=3.40158;     % ע�������չ�ұ�AR��(3,5)�ǱȽϺõġ�% ��Թ�Ӭ
% AR=(R_wing+xr0)/C_aver; 
% C_La_Ny=C_La2d/(E+k*C_La2d/(pi*AR));     % �����������������ۻ�õ���ά�����������б��  % ���: C_La_Ny =2.7506;
% % C_L_nabawy=(0.5*C_La2d/(E+k*C_La2d/(pi*AR)))*sin(2*alpha_4);
% C_L_nabawy=0.5*C_La_Ny*sin(2*abs(alpha2));   % ��2009-JEB-Rotational accelerations ��Fig.7ͼʵ�����ݶԱȷ��Ϻܺ�
% % C_D_nabawy=C_Lnabawy.*tan(alpha_4);   % C_Dnabawy=2*C_T*(sin(alpha))^2=C_La_Ny*(sin(alpha))^2;  
% C_D_nabawy=C_La_Ny*(sin(abs(alpha2))).^2;  
% C_L=C_L_nabawy;
% C_D=C_D_nabawy;
% C_N=cos(abs(alpha2)).*C_L+sin(abs(alpha2)).*C_D; %��������ϵ���ϳɡ���2010-JFM-RJ Wood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N=length(t);
% C_N=zeros(N,1);
% for i=1:1:N
%       C_N(i,1)=Cfcoeff(R_wing,C_aver,xr0,alpha2(i,1));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
% �������ϵ�µĽ��ٶ�
omega_x=dpsi;                     % չ��
omega_y=dphi.*sin(psi);       % ����(��ʼ����)
omega_z=dphi.*cos(psi);      % ������
% omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% �������ϵ�µĽǼ��ٶȡ����������������ļ���
domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������Ǻ͵��������ٶȵļ��㡪��ȡ�����ٶ�����ڸ�������ٶȣ���������Ӹ��ţ�
v_y_nonr=-omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;    % v_z=-r*dphi*sin(psi)
% alpha2=atan2(-v_y_nonr,v_z_nonr);   % ��ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ
% % ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % ���������ٶ�V_nonr=omega_h=dphi;   % ��λ�� rad/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ������ϵ�£�����������
% ��һ������������������ϵ�£�ƽ�������������������������������õ㣺ƬԪ�ҳ��е�
F_ytran=-sign(alpha2).*V_nonr.^2.*abs(C_N)*Coeff_liftdragF_N*10^(-3);   % ��λ��rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
% �ڶ������������� ������ϵ�£� ��ת�����������������������õ㣺ƬԪ�ҳ��е�
% C_R=1.55;
x_rd0=C_maxyaxis;
C_R=pi*(0.75-x_rd0);
F_yrot=C_R*omega_x.*V_nonr*F_yrotcoeff*10^(-3);      % ��λ��rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
% ����������������������ϵ�£��������������������������õ㣺ƬԪ�ҳ��е�
F_yadd1=(I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);  % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4  % �����ڶ�����
% �����ֳ���������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m_wing=2.4*10^-9;                    % kg
% x_com=1.920243385*10^-3;      %  m                           % ���ĵ�չ������
% % z_com=-0.149785466+0.636=0.486215*10^-3;      % ��Ťת����������
% z_com=0.149785466*10^-3;       %  m                          % ��Ťת����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis);        %���ú���inertia_moment;
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];
x_com=I_moment(1,1)*10^(-3);    % ���ĵ�չ������ % m;
z_com=I_moment(1,2)*10^(-3);    % ��Ťת���������� % m;
m_wing=I_moment(1,3)*10^(-3);  % �����òѧ�����仯�仯������� % kg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % ������������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(9)
% plot(t/T,F_ytran,'r-',t/T,F_yrot,'b-.',t/T,F_yadd1,'c-.',t/T,F_inert_y,'m-','LineWidth',2.5)
% legend('^{rw}\itF_{y,tran}','^{rw}\itF_{y,rot}','^{rw}\itF_{y,add}','^{rw}\itF_{y,inertia}')
% xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
% ylabel('Aerodynamic forces and inertia force ( \muN )','FontSize',20,'FontName','Times','FontWeight','Bold')
% % title('The nomal aerodynamic forces and inertial force in the wing planform fixed frame','FontSize',20,'FontWeight','Bold')
% title('��ƽ�������µķ����������͹�����','FontSize',20,'FontWeight','Bold')
% grid on
% % axis([t(1,1)/T,t(1,1)/T+1,min(F_ytran)-0.5,max(F_ytran)+0.5])
% axis([t(1,1)/T,t(1,1)/T+1,-inf,inf])
% box on
% set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ������������ ������ϵ�£���(�ϳ�)������������ƽ������������+��ת����������+�����������������õ㣺ƬԪ�ҳ��е�
% F_N=F_ytran+F_yadd1+F_yrot+F_inert_y;    % ��λ��rad^2*s^-2*kg*m=N=10^6uN
F_N=F_ytran+F_yadd1+F_yrot;    % ��λ��rad^2*s^-2*kg*m=N=10^6uN
B=[F_ytran,F_yadd1,F_yrot,F_inert_y];    % size(B)  % (1000*4)     %����mN
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal\optimal_wing_para_motion\Forcenormal_3T.xlsx',B,'sheet1','A1:D1000');
xlswrite('Forcenormal_3T.xlsx',B,'sheet1','A1:D1000'); % �����������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�ƬԪ��������ֽ⵽ƬԪ��ֱ�����ˮƽ�������õ㣺ƬԪ�ҳ��е�
F_vertical=-sign(alpha2).*F_N.*cos(alpha2);       % ��ֱ����vertical, %��λ��uN  % 
F_verticalaver=trapz(t,F_vertical)/(3*T);              % F_verticalaver=12.3074uN; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_horizontal=F_N.*sin(abs(alpha2));                   % ˮƽ����horizontal,%��λ��uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��������ϵ�£�����(drag)=F_horiz_y �Ͳ����(side force)F_horiz_X
% F_Z=F_vertical;
F_horiz_Y=cos(phi).*F_horizontal;                       % ˮƽ����_Y_axis����
% F_horiz_X=sin(phi).*F_horizontal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N1=length(F_horiz_Y);
for i=1:N1
    if F_horiz_Y(i)==0
        F_horiz_Y(i)=0.0000001;
    end
end
F_horiz_Y_aver=trapz(t,abs(F_horiz_Y))/(3*T);
Rario_F_vertical_horiz_aver=F_verticalaver/F_horiz_Y_aver    % ���Ʊ�
F_vertical_to_horiz_rms=sqrt(trapz(t,(abs(F_vertical)./(abs(F_horiz_Y))-Rario_F_vertical_horiz_aver).^2)/(3*T)) % ���Ʊ�rms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��������ֽ⵽��ֱ�����ˮƽ����ķ�������ʱ��ı仯������ʵ����Խ���ĶԱ�
% figure(10)
% % plot(t/T,F_vertical,'r-',t/T,-sign(alpha2).*F_horiz_Y,'b-.','LineWidth',2.5)
% plot(t/T,F_vertical,'r-',t/T,F_horiz_Y,'b-.','LineWidth',2.5)
% xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
% ylabel('Aerodynamic forces  ( \muN )','FontSize',20,'FontName','Times','FontWeight','Bold')
% legend('^{rr}\itF_{\itvertical,Z} ','^{rr}\itF_{\ithorizontal,Y}')
% % title('��������µĴ�ֱ�����ˮƽ����ķ�������ʱ��ı仯����','FontSize',20,'FontWeight','Bold')
% grid on
% t=t';
% % axis([min(t/T),min(t/T)+1,-12.5,17.5])
% % axis([t(1,1)/T,t(1,length(t))/T,min(F_horiz_Y)-0.5,max(F_vertical)+0.5])
% axis([t(1,1)/T,t(1,1)/T+1,min(F_horiz_Y)-0.5,max(F_vertical)+0.5])
% % set(gca,'XTick',(0.9:0.1:4.05))
% box on
% set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
% plot(t/T,F_ytran,'r-.',t/T,F_yrot,'b-.',t/T,F_yadd1,'g-.',t/T,F_inert_y,'c-.',t/T,F_vertical,'k-',t/T,F_horiz_Y,'k:','LineWidth',2.5)
% legend('^{rw}\itF_{y,tran}','^{rw}\itF_{y,rot}','^{rw}\itF_{y,add}','^{rw}\itF_{y,inertia}','^{rr}\itF_{\itvertical,Z} ','^{rr}\itF_{\ithorizontal,Y}')
plot(t/T,F_ytran,'b-.',t/T,F_yrot,'b:',t/T,F_yadd1,'g-.',t/T,F_vertical,'k-',t/T,F_horiz_Y,'k:','LineWidth',2.5)
legend('^{rw}\itF_{y,tran}','^{rw}\itF_{y,rot}','^{rw}\itF_{y,add}','^{rr}\itF_{\itvertical,Z} ','^{rr}\itF_{\ithorizontal,Y}')
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('Forces ( \muN )','FontSize',20,'FontName','Times','FontWeight','Bold')
% title('Aerodynamic forces and lift and thrust') % aerodynamic_forces_and_lift_and_thrust
grid on
% axis([t(1,1)/T,t(1,1)/T+1,min(F_horiz_Y)-0.5,max(F_vertical)+0.5])
% axis([t(1,1)/T,t(1,1)/T+1,min(F_ytran)-0.5,max(F_vertical)+0.5])
axis([t(1,1)/T,t(1,1)/T+1,-inf,inf])
box on
set(gca,'YGrid','off','XTick',(min(t/T):0.025:min(t/T)+1),'LineStyle','-','LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ʵļ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һģ�顪��չ��Ťת�ᡪ���������ط���%%%%
%%%%���߲���, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡���ƽ��������������������%%%%%%
%%%%�ڶ����֡���ת��������������������%%%%%%
%%%%�������֡���ת��������������%%%%%%%%%
%%%%���Ĳ��֡���ת������������%%%%%%%%%%
%%%%���岿�֡���Ťת�����ĵ��Իظ�����%%%%%%
%%%%�������֡������������%%%%%%%%%%%%
%%%%���߲��֡�����������������%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡���ƽ�������������������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
%% (3) ����3���������ת���������ء�����ѹ��λ��
k_xaero=1; 
N=length(t);
M_xtrans=zeros(N,1);
Y_rcpnd_trans=zeros(N,1);
for i=1:1:N
    Y_rcpnd_trans(i,1)=COP_Ycpnd2_TransCirc(alpha2(i,1),R_wing,xr0,C_maxyaxis); %���ú���COP_Ycpnd2_TransCirc��⾻ѹ�ĵ�������λ��Y_rcpnd; % ��������
    %  Y_rcpnd_trans=abs(Y_rcpnd_trans(i,1));  % ��
    % ����ĵ�λ�� (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
    %  M_xaercoeff=0.0060;   %��λ��: mg.mm^2
    % ƽ��������������ת����������һ��ʼ����ʱ���
     % M_xtrans(i,1)=-k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);      
     M_xtrans(i,1)=k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);  
end
%% �ڶ����֡���ת�������������������ء���Ťת��
k_xRot=1;  
% C_R=1.55;    % ��ϵ�������޸�
% x_rd0=C_maxyaxis;
% C_R=pi*(0.75-x_rd0);
M_xRotcoeff=k_xRot*C_R*M_xRotcoeff;
N=length(t);
M_xrotcirc=zeros(N,1);
Y_rcpnd_rot=zeros(N,1);
for i=1:1:N
    % ת��������Ťת���������ء������ú���COP_Ycpnd2_RotCirc��⾻ѹ�ĵ�������λ��Y_rcpnd_rot
    % ת�����������ġ�ѹ�ķֲ�����Dickinson���� or ѹ�������ҵ� or ѹ����c(r)/4����Ťת������
    Y_rcpnd_rot(i,1)=COP_Ycpnd2_RotCirc(alpha2(i,1),R_wing,xr0,C_maxyaxis); % ѹ�ķֲ�����Dickinson����: 
   % Y_rcpnd_rot=abs(Y_rcpnd_rot(i,1));
   % ����ĵ�λ�� (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
   %  M_xrotcirc(i,1)=k_xRot*C_R*sign(alpha2(i,1)).*omega_x(i,1).*abs(V_nonr(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3); % ת��������Ťת����������  
   M_xrotcirc(i,1)=-k_xRot*C_R*omega_x(i,1).*V_nonr(i,1).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3);    % ת��������Ťת����������
end
%% �������֡���ת�������������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
C_RD=1;  
M_xrd=-C_RD*omega_x.*abs(omega_x)*M_xrdcoeff*10^(-3);   % ��������һ��ʼ��˳ʱ��� % ע���������ȡ����omega_x
%% ���Ĳ��֡���ת������������+ƽ�����������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
k_am=1;   % ����������ϵ��; 2.35�ǽ������ޣ� psi_max =32.9386
% M_xam=-k_am*(-I_xzam*(domega_z+omega_x.*omega_y)+I_xxam*domega_x)*10^(-3);   % ��ʼ��ʱ��(-)
M_xam=k_am*(-I_xzam*(domega_z+omega_x.*omega_y)-I_xxam*domega_x)*10^(-3);   % ��ʼ��ʱ��(-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(101)
plot(t/T,M_xtrans,'r--',t/T,M_xrotcirc,'b-',t/T,M_xrd,'g--',t/T,M_xam,'c-','LineWidth',2)
legend('^{rw}\itM_{x,trans}','^{rw}\itM_{x,rot}','^{rw}\itM_{x,rd}','^{rw}\itM_{x,add}')
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('Aerodynamic moments  ( \muN.mm )','FontSize',20,'FontName','Times','FontWeight','Bold')
title('����Ťת��^{rw}\itx-axis\rm����������','FontSize',20,'FontWeight','Bold')
grid on
% axis([t(1,1)/T,t(1,1)/T+1,min(M_xtrans)-1,max(M_xtrans)+1])
axis([min(t/T),min(t/T)+1,-inf,inf])
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��Ťת��Ĺ��ʼ���
% ƽ���������ع���  % uW
P_xtrans=M_xtrans.*omega_x*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% ת���������ع���  % uW
P_xrd=M_xrd.*omega_x*10^-3;  
% ת���������ع���  % uW  
P_xrotcirc=M_xrotcirc.*omega_x*10^-3; 
% ���������ع��ʺͳ������������ع���  % uW
P_xam=M_xam.*omega_x*10^-3;
P_aerox=P_xtrans+P_xrd+P_xrotcirc+P_xam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(102)
% plot(t/T,P_xtrans,'r--',t/T,P_xrotcirc,'b-',t/T,P_xrd,'g--',t/T,P_xam,'c-',t/T,P_aerox,'k-','LineWidth',2)
% legend('ƽ������','ת������','ת�����Ṧ��','����������','Ťת���ܹ���')
% title('����Ťת��x-axis����������')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ�ģ�顪������ת���ᡪ���������ط���%%%%%%
%%%%���߲���, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡���ƽ��������������������%%%%%%%
%%%%�ڶ����֡���ת��������������%%%%%%%%%%%
%%%%�������֡���ת������������%%%%%%%%%%%%
%%%%���Ĳ��֡���ת����������%%%%%%%%%%%%%
%%%%���岿�֡����Ĵ���ת�������ĵ��Իظ�����%%%%
%%%%�������֡����Ĵ�����������%%%%%%%%%%%%
%%%%���߲��֡�����������������%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���� ƽ�������������������������ء����Ƴ�ƽ���µ�������
k_ztrans=1;  
M_ztrans=k_ztrans*sign(alpha2).*(I1z.*abs(C_N).*V_nonr.^2)*10^(-3);   % I1z=0.0162; % ��λ�� mg.mm^2=10^-3uN*mm=10^-6 mN*mm
%% �ڶ����� ת�������������������������ء����Ƴ�ƽ���µ�������
C_R=1.55;    % ��ϵ��������
k_zrot=1;
M_zrot=-k_zrot*I2z*C_R*omega_x.*V_nonr*10^(-3);  % I2z=0.0069; ��λ��: mg.mm^2*(rad*s^-1)^2=10^-6mN.mm=10^-3uN.mm
%% �������� �������������������������ء����Ƴ�ƽ���µ������� % I5y��I6y��λ�� mg.mm=*10^(-9) kg.m  %��λ��10^6uN
% k_za=0.35; 
k_za=1;  % ��ϵ��������ô��
% M_zadd=-k_za*(-I7z*(domega_z+omega_x.*omega_y)+I_xzam*domega_x)*10^(-3); % ԭ���Ƴ��Ĺ�ʽ������Ϊ����ЩŶ   % ��ʼ��ʱ��(-)
M_zadd=k_za*(-I7z*(domega_z+omega_x.*omega_y)-I_xzam*domega_x)*10^(-3); 
%% ���Ĳ��� ת���������ء����Ƴ�ƽ���µ�������
% C_RD2=0.05*C_RD;
C_RD2=1*C_RD;  % ��ϵ��������ô��
M_zrd=-C_RD2*omega_x.*abs(omega_x)*M_zrdcoeff*10^(-3);   % ��������һ��ʼ��˳ʱ��� % ע���������ȡ����omega_x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(103)
plot(t/T,M_ztrans,'r--',t/T,M_zrot,'b-',t/T,M_zrd,'g--',t/T,M_zadd,'c-','LineWidth',2)
legend('^{rw}\itM_{Z,trans}','^{rw}\itM_{Z,rot}','^{rw}\itM_{Z,rd}','^{rw}\itM_{Z,add}')
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('Aerodynamic moments  ( \muN.mm )','FontSize',20,'FontName','Times','FontWeight','Bold')
title('�����Ĵ���^{rr}\itZ-axis\rm����������','FontSize',20,'FontWeight','Bold')
grid on
% axis([t(1,1)/T,t(1,1)/T+1,min(M_ztrans)-1,max(M_ztrans)+1])
axis([min(t/T),min(t/T)+1,-inf,inf])
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �Ĵ���Ĺ���
% % ƽ���������ع���  % uW
% P_ztrans=M_ztrans.*omega_z*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ת���������ع���  % uW  
% P_zrotcirc=M_zrot.*omega_z*10^-3; 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ���������ع��ʺͳ������������ع���  % uW
% P_zam=M_zadd.*omega_z*10^-3;
% % P_inert_z=M_inert_z.*omega_z*10^-3;
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ת���������ع���  % uW
% P_zrd=M_zrd.*omega_z*10^-3;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ƽ���������ع���  % uW
P_ztrans=cos(psi).*M_ztrans.*dphi*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
%%%%%%%%%%%%%%%%%%%%%%%%%
% ת���������ع���  % uW  
P_zrotcirc=cos(psi).*M_zrot.*dphi*10^-3; 
%%%%%%%%%%%%%%%%%%%%%%%%%
% ���������ع��ʺͳ������������ع���  % uW
P_zam=cos(psi).*M_zadd.*dphi*10^-3;
% P_inert_z=M_inert_z.*omega_z*10^-3;
%%%%%%%%%%%%%%%%%%%%%%%%%
% ת���������ع���  % uW
P_zrd=cos(psi).*M_zrd.*dphi*10^-3;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_aeroz=P_ztrans+P_zrd+P_zrotcirc+P_zam;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(104)
% plot(t/T,P_ztrans,'r--',t/T,P_zrotcirc,'b-',t/T,P_zrd,'g--',t/T,P_zam,'c-',t/T,P_aeroz,'k-','LineWidth',2)
% legend('ƽ������','ת������','ת�����Ṧ��','����������','�Ĵ����ܹ���')
% title('�����Ĵ���z-axis����������')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��������ɶ��˶��Ĺ��ʼ��㡪��ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aerox=P_xtrans+P_xrd+P_xrotcirc+P_xam;
% P_totalx=P_xtrans+P_xrd+P_xrotcirc+P_xam+P_inert_x;
P_psix_total=P_aerox;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aeroz=P_ztrans+P_zrd+P_zrotcirc+P_zam;
% P_totalz=P_ztrans+P_zrd+P_zrotcirc+P_zam+P_inert_z;
P_phiz_total=P_aeroz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I_moment=inertia_moment(R_wing,C_aver,xr0,C_maxyaxis);        %���ú���inertia_moment;
% uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% I_moment=[XC_tran,YC_tran,M_big,Ix_inertia,Iz_inertia];
Ix_inertia=I_moment(1,4)*10^(-3);     % ��ƽ���ᶨ����  % g.mm^2=10^(-9)kg.m^2����for 10^-3uW
Iz_inertia=I_moment(1,5)*10^(-3);     % ��ƽ���ᶨ����  % g.mm^2=10^(-9)kg.m^2����for 10^-3uW
% ���桪��rad*s^-1*10^(-9)kg.m^2*rad*s^-2=10^(-9)(N=kg*m.s^-2)*(m.s^-1)=10^-9W=10^-6mW=10^-3uW
P_Ix_inertia=omega_x.*Ix_inertia.*domega_x;      % uW
P_Iz_inertia=dphi.*Iz_inertia.*ddphi;   % uW     % cos(psi).*
P_totalx=P_psix_total+P_Ix_inertia;
P_totalz=P_phiz_total+P_Iz_inertia;
% P_total=[P_totalx,P_totalz];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_I_inertia=P_Ix_inertia+P_Iz_inertia;
figure(105)
plot(t/T,P_Ix_inertia,'r--',t/T,P_Iz_inertia,'b-.',t/T,P_I_inertia,'k-','LineWidth',2)
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('Power output ( \muW )','FontSize',20,'FontName','Times','FontWeight','Bold')
legend('^{rw}\itP_{x,inertia}','^{rr}\itP_{Z,inertia}','\itP_{total,inertia}')
title('����Ťת��^{rw}\itx-axis\rm���Ĵ���^{rr}\itZ-axis\rm�Ĺ��Թ���','FontSize',20,'FontWeight','Bold')
grid on
axis([min(t/T),min(t/T)+1,-inf,inf])
box on
set(gca,'LineStyle','-','LineWidth',1.5,'FontSize',16,'FontName','Times','FontWeight','Bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aero_aver=trapz(t,(P_psix_total+P_phiz_total))/(3*T)
% P_inertia_aver=trapz(t,(P_Ix_inertia+P_Iz_inertia))/(3*T)
% P_aero_aver_to_inertia_aver=P_aero_aver/P_inertia_aver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_total=Aero_M_fruitfly2_exp(x);      % ���ú������������ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
% % size(P_total)              % (1000*2)
% P_totalx=P_total(:,1);    % (1000*1)
% P_totalz=P_total(:,2);    % (1000*1)
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total.xlsx',P_total,'sheet1','A1:B2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_psix_total_posi=P_psix_total;% ����Ťת�����ܹ���
N=length(P_psix_total_posi);
for i=1:N
    if P_psix_total_posi(i)<=0
        P_psix_total_posi(i)=0;
    end
end
%%%%%%%%%%%%%%%%%%%
P_Ix_inertia_posi=P_Ix_inertia;% ����Ťת���Թ���
for i=1:N % N=length(P_Ix_inertia_posi);
    if P_Ix_inertia_posi(i)<=0
        P_Ix_inertia_posi(i)=0;
    end
end
%%%%%%%%%%%%%%%%%%%
P_totalx_posi=P_totalx; % ����Ťת�ܹ���=�����ܹ���+���Թ���
for i=1:N %N=length(P_totalx_posi);
    if P_totalx_posi(i)<=0
        P_totalx_posi(i)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_phiz_total_posi=P_phiz_total;% �����Ĵ������ܹ���
for j=1:N  % N=length(P_totalx_posi);
    if P_phiz_total_posi(j)<=0
        P_phiz_total_posi(j)=0;
    end
end
%%%%%%%%%%%%%%%%%%%
P_Iz_inertia_posi=P_Iz_inertia;% �����Ĵ���Թ���
for j=1:N  % N=length(P_totalx_posi);
    if P_Iz_inertia_posi(j)<=0
        P_Iz_inertia_posi(j)=0;
    end
end
%%%%%%%%%%%%%%%%%%%
P_totalz_posi=P_totalz;% �����Ĵ��ܹ���=�����ܹ���+���Թ���
for j=1:N  % N=length(P_totalx_posi);
    if P_totalz_posi(j)<=0
        P_totalz_posi(j)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P_aero_aver_posi=trapz(t,(P_psix_total_posi+P_phiz_total_posi))/(3*T);
P_inertia_aver_posi=trapz(t,(P_Ix_inertia_posi+P_Iz_inertia_posi))/(3*T);
P_aero_aver_to_inertia_aver=P_aero_aver_posi/P_inertia_aver_posi       % ȫ������������ȫ�����Թ��ʵı�ֵ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(106)
P_x=P_aerox+P_Ix_inertia;
plot(t/T,P_xtrans,'r-.',t/T,P_xrotcirc,'b-.',t/T,P_xrd,'g-.',t/T,P_xam,'c-.',t/T,P_aerox,'r-',t/T,P_Ix_inertia,'m-',t/T,P_x,'k-',t/T,P_totalx_posi,'g-','LineWidth',2)
% xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('\rmPitch power output ( \muW )','FontSize',20,'FontName','Times','FontWeight','Bold')
legend('\itP_{x,trans}','\itP_{x,rot}','\itP_{x,rd}','\itP_{x,add}','\itP_{x,aero}','\itP_{x,inertia}','\itP_{x,total}','\itP_{x,total,posi}')
% title('����Ťת��^{rw}\itx-axis\rm���������ʺ͹��Թ���','FontSize',20,'FontWeight','Bold')
grid on
% axis([min(t/T),min(t/T)+1,min(P_aerox)-1,max(P_Ix_inertia)+1])
axis([min(t/T),min(t/T)+1,-inf,inf])
box on
set(gca,'YGrid','off','XTick',(min(t/T):0.025:min(t/T)+1),'LineStyle','-','LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(107)
P_z=P_aeroz+P_Iz_inertia;
plot(t/T,P_ztrans,'r-.',t/T,P_zrotcirc,'b-.',t/T,P_zrd,'g-.',t/T,P_zam,'c-.',t/T,P_aeroz,'r-',t/T,P_Iz_inertia,'m-',t/T,P_z,'k-',t/T,P_totalz_posi,'g-','LineWidth',2)
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('\rmFlapping power output ( \muW )','FontSize',20,'FontName','Times','FontWeight','Bold')
legend('\itP_{Z,trans}','\itP_{Z,rot}','\itP_{Z,rd}','\itP_{Z,add}','\itP_{Z,aero}','\itP_{Z,inertia}','\itP_{Z,total}','\itP_{Z,total,posi}')
% title('�����Ĵ���^{rr}\itZ-axis\rm���������ʺ͹��Թ���','FontSize',20,'FontWeight','Bold')
grid on
axis([min(t/T),min(t/T)+1,-inf,inf])
box on
set(gca,'YGrid','off','XTick',(min(t/T):0.025:min(t/T)+1),'LineStyle','-','LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ų����ò�����ĳ�Ťת���ʡ��Ĵ��ʺ��ܹ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(108)
plot(t/T,P_totalx_posi,'r--',t/T,P_totalz_posi,'b-','LineWidth',2)
hold on
P_total=P_totalx_posi+P_totalz_posi;
plot(t/T,P_total,'k-','LineWidth',2.5)
xlabel('Normalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('Power output ( \muW )','FontSize',20,'FontName','Times','FontWeight','Bold')
% legend('\itP_{\itx,\psi,posi}','\itP_{\itZ,\phi,posi}','\itP_{\ittotal,,posi}')
% legend('\it^{\itrw}P_{\itx,posi}','\it^{\itrr}P_{\itZ,posi}','\itP_{\ittotal,,posi}')
legend('\itP_{\itx,posi}','\itP_{\itZ,posi}','\itP_{\ittotal,,posi}')
title('������ų����ò�����ĳ�Ťת���ʡ��Ĵ��ʺ��ܹ������','FontSize',20,'FontWeight','Bold')   
grid on
% axis([t(1,1)/T,t(1,length(t))/T,0,max(P_total)+0.5])
axis([t(1,1)/T,t(1,1)/T+1,0,max(P_total)+0.5])
box on
% set(gca,'XTick',(0.9:0.1:4.05))
set(gca,'YGrid','off','XTick',(min(t/T):0.025:min(t/T)+1),'LineStyle','-','LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_total_posi=[P_totalx_posi,P_totalz_posi]; 
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\optimal_wing_para_motion\P_total_positive.xlsx',P_total_posi,'sheet1','A1:B1000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=188.7; 
% f=x(1);
% T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
% t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1����XXXXXX
% t_00= 1/(4*f);    % t_00= -1/(4*f);
% t=linspace(t_00,t_00+3*T,1000); 
P_totalx_aver=trapz(t,P_totalx_posi)/(3*T);  % ʱ�������ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
P_totalz_aver=trapz(t,P_totalz_posi)/(3*T);  % ʱ�������ʡ���ƽ�����ʺ�Ťת���ʡ���Ťת�Ṧ�ʺ��Ĵ��Ṧ��
P_total_aver=P_totalx_aver+P_totalz_aver;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_M=[F_verticalaver,P_total_aver];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%