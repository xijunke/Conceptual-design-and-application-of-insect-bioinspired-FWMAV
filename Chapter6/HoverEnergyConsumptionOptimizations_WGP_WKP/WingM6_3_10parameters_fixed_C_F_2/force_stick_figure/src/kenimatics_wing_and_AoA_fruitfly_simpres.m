function wing_m_output=kenimatics_wing_and_AoA_fruitfly_simpres()
%% kenimatics_wing_and_AoA
% ���˶����ɺͼ��ι���(AOA)
% (ע�⣺���ռ��������Ťת��,��Ҫʱ���������ʵ��ʱ������׼)��������Ҫ���ͺ˶�
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ź������󵼡�����ý����ʡ������������ڼ���������������ʱ����Ҫ���¼�����õ�Ťת��
%% (1) ��Ϊ�趨�Ĵ�Ǻ���
syms t         % t_range=[0.0052824335,0.0052824335+5*T];
w =1185.6;     % ��Ƶ��   %  f=188.7; T=1/f;  %����Ƶ�� (Hz)������ 
% f_var=188.7;  % Scienceʵ���Ĵ���Ľ׸���Ҷ���������Ļ�Ƶ
% w =2*pi*f_var;   PHI=1.1487;
% phi_pres =PHI*cos(w*t);     
% ��Ϊ�趨�Ĵ��(��ֵΪPHI, Ƶ��ΪScienceʵ���Ĵ���Ľ׸���Ҷ���������Ļ�Ƶ188.7Hz) % ��Ϊ�趨�Ĵ��: prescribed=�涨��
phi_pres =sym('1.149*sin(w*t+1.571)');  %���ź���
dphi_pres =diff(phi_pres ,t,1);
ddphi_pres =diff(phi_pres ,t,2); 
dphi=inline(vectorize(dphi_pres ),'w','t');                    % ��ֵ����
ddphi=inline(vectorize(ddphi_pres ),'w','t');  
% dphi_pres =-w*PHI*sin(w*t);
% ddphi_pres =-w^2*PHI*cos(w*t);   
%% (2) ���������ģ�����õ���Ťת�Ǻͽ��ٶ��Լ��Ǽ��ٶ�    
psi_sim=sym('-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)-0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180');  %���ź���
dpsi_sim=diff(psi_sim,t,1);   
ddpsi_sim=diff(psi_sim,t,2);    
dpsi=inline(vectorize(dpsi_sim),'w','t');        %��ֵ����
ddpsi=inline(vectorize(ddpsi_sim),'w','t');    %��ֵ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������ֵ
f=188.7; T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1   
dphi=dphi(w,t);         %(1*100)�������������Ĵ������                                                           % ���
ddphi=ddphi(w,t);    %(1*100)�������������Ĵ�Ǽ�����                                                         % ���
dpsi=dpsi(w,t);         %(1*100)������������Ťת������                                                             % ���
ddpsi=ddpsi(w,t);     %(1*100)������������Ťת�Ǽ�����                                                         % ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �Ĵ�Ǻ�Ťת��
% �Ĵ�ǡ��� (1*100)�������������Ĵ�ǡ���������            %���                                     
phi_pres=1.149*sin(w*t+1.571);
% Ťת�ǡ���(1*100)������������Ťת�ǡ���������            %���
psi_sim=-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+...
              0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)...
             -0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180;   
% [phi_min,k]=min(phi);  % ���:  phi_min =-1.0157;  k =10756;      
% t_0=t(k);                       % ���: t0 =0.0028;  
% figure(1)                                                          % ͼ1�����Ĵ�Ǻ�Ťת��
% subplot(311)
% plot(t*f,psi_sim*180/pi,'g:',t*f,phi_pres*180/pi,'k-.','LineWidth',2) %ת��Ϊms �� ����degree   *10^3   *180/pi
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) vs \phi_{pres}(t)');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)');  
% title('����ת����\psi_{sim}(t)���Ĵ��phi_{pres}(t)��ʱ��ı仯')  
% grid on  % ����ת����psi_sim(t)��Ťת��phi_pres(t)��ʱ��ı仯
% axis([0.9,4.05,-105,105])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-�Ĵ�����ʺ�Ťת������
% subplot(312)
% plot(t/T,dphi,'r-',t/T,dpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('������ (rad/s)')
% legend('\itd\phi(t)','\itd\psi(t)')
% title('�Ĵ�����ʺ�Ťת��������ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
% Plot-�Ĵ�Ǽ����ʺ�Ťת�ǽǼ�����
% subplot(313)
% plot(t/T,ddphi,'r-',t/T,ddpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('�Ǽ����� (rad/s^-2)')
% legend('\itdd\phi(t)','\itdd\psi(t)')
% title('�Ĵ�Ǽ����ʺ�Ťת�Ǽ�������ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-Ťת�Ǻͼ��ι���AOA
alpha1=pi/2-psi_sim.*sign(dphi);        %(1*200)���������������ι��ǡ���������        %�������ȫ�����ι��� 
% Y = sign(X) returns an array Y the same size as X, where each element of Y is:
% *1 if the corresponding element of X is greater than zero
% * 0 if the corresponding element of X equals zero
% *-1 if the corresponding element of X is less than zero
% figure(2)                                                              % ͼ2���� ע�����Ｘ�ι���ʼ��ȡ��ֵ
% % % x_interval=[0,1/2,1/2,0];
% % % y_interval=[-100,-100,100,100];
% % % fill(x_interval,y_interval,'y');
% % % hold on
% % % legend('','\it\psi(t)','\it\alpha_1(t)')
% plot(t/T,psi_sim*180/pi,'b-',t/T,alpha1*180/pi,'g-')      %Ťת�Ǻͼ��ι���AOA��ʱ��ı仯����
% xlabel('\itNormalized time')
% ylabel('\itAngle (��)')
% legend('\it\psi_{sim}(t)','\it\alpha_1(t)')
% % hold on
% % L=length(t);
% % plot([0,t(L)/T],[0,0],'k-');     %��x-axis
% title('Ťת�Ǻͼ��ι���alpha_1��ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-105,105])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ��ٶ�
% omega_x=dpsi;   % Ťת��:   x(1)=psi;  x(2)=dpsi;   
omega_y=dphi.*sin(psi_sim);
omega_z=dphi.*cos(psi_sim);
% omega_h=dphi;    % �������ϵ�½���������    % omega_h=sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
%% �������ϵ�µĽǼ��ٶȡ����������������ļ���
% domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi); 
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% �������Ǽ���
v_y_nonr=omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
alpha2=atan2(v_y_nonr,-v_z_nonr);   % ��ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ
% ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)                                               % ͼ3����ʹ�ó������ϵ�µĽ��ٶ���⼸�ι���
% plot(t/T,alpha1*180/pi,'r-',t/T,alpha2*180/pi,'b-','LineWidth',2)    
% xlabel('\itNormalized time')
% ylabel('\it\alpha_1 & \alpha_2 (deg)')
% legend('\alpha_1 (t)','\alpha_2 (t)')
% title('������ʱ��ı仯����')   % ������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������ϵ�£�������ϵ��
%% lift and drag coefficients with Alpha�����湥�Ǳ仯��������ϵ���ͷ���������ϵ��
% alpha=45;     %�ٶ����Ǻ㶨����;���湫ʽ��������ϵ��Ϊ:C_L=1.8046;C_D=1.7037.
alpha=alpha2;
% (1)����ľ��鹫ʽ����1999-science-MH Dickinson�����ס�����һ��
% C_L =0.225+1.58*sin((2.13*alpha*180/pi-7.2)*pi/180);  % ���ι�����������(alpha),����ʹ�þ���ֵ����abs(alpha);
% C_D =1.92-1.55*cos((2.04*alpha*180/pi-9.82)*pi/180);  % ���ι�����������(alpha),����ʹ�þ���ֵ����abs(alpha);
% C_N=cos(alpha).*C_L2+sin(alpha).*C_D2;  % ���ι�����������(alpha2),����ʹ�þ���ֵ����abs(alpha);
% (2) lift and drag coefficients with Alpha�����湥�Ǳ仯��������ϵ���ͷ���������ϵ�����ڶ��֡������仯��ϵ��
% (a) 1999-Science-MH Dickinson�����Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�
C_L =(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180)); 
C_D =(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));
% (b) 2004-JEB-Wang ZJ�����Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�
% C_D =(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));
% C_D =1.92+1.55*cos((2.04*alpha1-9.82)*pi/180);
% (c) ƽ��������ϵ�������ֵ
% C_L_aver=trapz(t,C_L)/(3*T); % C_L_aver =1.5088;
% C_D_aver=trapz(t,C_D)/(3*T) % C_D_aver =1.9123;
% C_L2C_D_aver=C_L_aver/C_D_aver; % C_L2C_D_aver =0.7890;
C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %��������ϵ���ϳɡ���2010-JFM-RJ Wood
% (3) lift and drag coefficients with Alpha�����湥�Ǳ仯��������ϵ���ͷ���������ϵ�����������֡���ȫ��ϵ��
C_L2 =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  % ���ι���ȫ��abs(alpha);��������ϵ��ȫ��
C_D2 =1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180); % ���ι���ȫ��abs(alpha);��������ϵ��ȫ��
% C_N2=cos(abs(alpha)).*C_L2+sin(abs(alpha)).*C_D2;                % ���ι���Ϊ��abs(alpha)����������ϵ��ȫ��
% ��������-��������ϵ���ϳɡ��� 2012-IEEE ASME TM-Veaceslav Arabagi ��2009-Science-Deng Xinyan
C_N2=sign(alpha2).*sqrt(C_L2.^2+C_D2.^2); 
% (4) Normal and tangential coefficients with Alpha�����湥�Ǳ仯�ķ��������������ϵ�����������֡���ȫ��ϵ��
C_N3=3.4*sin(abs(alpha));      % 2006-IEEE TR-Deng Xinyan �� 2010-EAE-RJ Wood-����������ϵ�����
% ����(tangential)������ϵ��C_Tֻ��alpha��(-pi/4,pi/4)ʱ��Ϊ�㣬�������Ϊ��
% if alpha<pi/4 & alpha>-pi/4
%     C_T=0.4*cos(2*alpha).^2
% else
%     C_T=0;
% end
C_T=0.4*(cos(2*alpha)).^2.*(alpha>-pi/4 & alpha<pi/4);       %����������ϵ�����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)        % ͼ6��������ϵ�¡�������ϵ��
% plot(t/T,C_L,'r-',t/T,C_D,'b-','LineWidth',2);      %ʱ�����򻯣�����Ƶ��f, ���߳�������T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of lift and drag \itvs. t \rm for flapping wing')
% title('����������ϵ������ʱ��ı仯����')
% legend('\itC_L','\itC_D')
% grid on
% axis([0.9,4.05,-inf,inf])
% figure(5)     % ͼ7��������ϵ�¡����������������ϵ��
% plot(t/T,C_N3,'r-',t/T,C_T,'b-')                                % �õ��ķ��������������ϵ����Ҫ��ʵ
% xlabel('Normalized time')
% ylabel('C_N3(\alpha(t)) & C_T(\alpha (t))')
% legend('C_{N3}(\alpha)','C_T(\alpha)')
% title('���������������ϵ����ʱ��ı仯����')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %���ַ���������ϵ��
% figure(6)     % ͼ8��������ϵ�¡�����������ϵ��
% plot(t/T,C_N1,'r-',t/T,C_N2,'b-',t/T,C_N3,'g-','LineWidth',2);      %ʱ�����򻯣�����Ƶ��f, ���߳�������T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of normal force \itvs. t \rm for flapping wing')
% title('����������ϵ����ʱ��ı仯����')
% legend('\itC_{N1}','\itC_{N2}','\itC_{N3}')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���˶����������������ϵ�����
% wing_m_output=zeros(1000,12);
wing_m_output=[t',phi_pres',psi_sim',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

