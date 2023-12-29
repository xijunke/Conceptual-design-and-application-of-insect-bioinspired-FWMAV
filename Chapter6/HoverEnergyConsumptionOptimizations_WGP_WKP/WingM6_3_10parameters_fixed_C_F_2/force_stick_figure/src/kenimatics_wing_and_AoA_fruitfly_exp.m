function wing_m_output=kenimatics_wing_and_AoA_fruitfly_exp()
%% kenimatics_wing_and_AoA
% ���˶����ɺͼ��ι���(AOA)
% (ע�⣺���ռ��������Ťת��,��Ҫʱ���������ʵ��ʱ������׼)��������Ҫ���ͺ˶�
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ź������󵼡�����ý����ʡ������������ڼ���������������ʱ����Ҫ���¼�����õ�Ťת��
% (1) 2014-Science-ʵ����Ĵ�ǡ��� ����4�׸���Ҷ������� 
syms t         % t  �Լ�����ʱ����
w =1185.6;     % ��Ƶ��   %  f=188.7; T=1/f;  %����Ƶ�� (Hz)������  
phi1=sym('(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180');  %���ź���: ԭʼ���ݵĶ�������ת��Ϊ������
dphi1=diff(phi1,t,1);
ddphi1=diff(phi1,t,2); 
dphi=inline(vectorize(dphi1),'w','t');                    % ��ֵ����
ddphi=inline(vectorize(ddphi1),'w','t');   
% (2) 2014-Science-ʵ���Ťת�ǡ�������8�׸���Ҷ�������    
psi1=sym('-(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180');  %���ź���:ԭʼ���ݵĶ�������ת��Ϊ������
dpsi1=diff(psi1,t,1);   
ddpsi1=diff(psi1,t,2);    
dpsi=inline(vectorize(dpsi1),'w','t');  %��ֵ����
ddpsi=inline(vectorize(ddpsi1),'w','t');  %��ֵ����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
phi=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+...
        0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180; % ԭʼ���ݵĶ�������ת��Ϊ������
% Ťת�ǡ���(1*100)������������Ťת�ǡ���������            %���
psi=-(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+... 
        3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+...
        0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180;  % ԭʼ���ݵĶ�������ת��Ϊ������
% [phi_min,k]=min(phi);  % ���:  phi_min =-1.0157;  k =10756;      
% t_0=t(k);                       % ���: t0 =0.0028;  
% figure(1)       % ͼ1�����Ĵ�Ǻ�Ťת��
% % subplot(311)
% plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-','LineWidth',2)  %ת��Ϊms �� ����degree   *10^3   *180/pi
% xlabel('\itNormalized time')
% ylabel('\itAngle (��)')
% legend('\it\phi(t)','\it\psi(t)')
% title('�Ĵ�Ǻ�Ťת����ʱ��ı仯����')   % �Ĵ�Ǻ�Ťת����ʱ��ı仯����
% grid on
% axis([0.9,4.05,-105,105])
% set(gca,'XTick',(0.9:0.1:4.05))
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
% set(gca,'XTick',(0.9:0.1:4.05))
% % Plot-�Ĵ�Ǽ����ʺ�Ťת�ǽǼ�����
% subplot(313)
% plot(t/T,ddphi,'r-',t/T,ddpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('�Ǽ����� (rad/s^-2)')
% legend('\itdd\phi(t)','\itdd\psi(t)')
% title('�Ĵ�Ǽ����ʺ�Ťת�Ǽ�������ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-Ťת�Ǻͼ��ι���AOA
% % alpha1=pi/2+psi.*sign(dphi);        %(1*200)���������������ι��ǡ���������        %�������ȫ�����ι��� 
% alpha1=pi/2-psi.*sign(dphi);        %(1*200)���������������ι��ǡ���������        %�������ȫ�����ι���  % ʵ��Ťת��ǰ����
% % Y = sign(X) returns an array Y the same size as X, where each element of Y is:
% % *1 if the corresponding element of X is greater than zero
% % * 0 if the corresponding element of X equals zero
% % *-1 if the corresponding element of X is less than zero
% figure(2)                                                              % ͼ2���� ע�����Ｘ�ι���ʼ��ȡ��ֵ
% % % x_interval=[0,1/2,1/2,0];
% % % y_interval=[-100,-100,100,100];
% % % fill(x_interval,y_interval,'y');
% % % hold on
% % % legend('','\it\psi(t)','\it\alpha_1(t)')
% plot(t/T,psi*180/pi,'b-',t/T,alpha1*180/pi,'g-','LineWidth',2)      %Ťת�Ǻͼ��ι���AOA��ʱ��ı仯����
% xlabel('\itNormalized time')
% ylabel('\itAngle (��)')
% legend('\it\psi(t)','\it\alpha_1(t)')
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-');     %��x-axis
% title('Ťת�Ǻͼ��ι���alpha_1��ʱ��ı仯����') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ��ٶ�
% omega_x=dpsi;   % Ťת��:   x(1)=psi;  x(2)=dpsi;   
omega_y=dphi.*sin(psi);
omega_z=dphi.*cos(psi);
% omega_h=dphi;    % �������ϵ�½���������    % omega_h=sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
%% �������ϵ�µĽǼ��ٶȡ����������������ļ���
% domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi); 
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% �������Ǽ���
v_y_nonr=omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
alpha2=atan2(v_y_nonr,-v_z_nonr);   % ��ȷ����ע�������ĵ�alpha=atan2(omega_z,-omega_y)*180/pi; ��ͬ  % ʵ��Ťת��ǰ����
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
% C_L2 =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  % ���ι���ȫ��abs(alpha);��������ϵ��ȫ��
% C_D2 =1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180); % ���ι���ȫ��abs(alpha);��������ϵ��ȫ��
% % C_N2=cos(abs(alpha)).*C_L2+sin(abs(alpha)).*C_D2;                % ���ι���Ϊ��abs(alpha)����������ϵ��ȫ��
% ��������-��������ϵ���ϳɡ��� 2012-IEEE ASME TM-Veaceslav Arabagi ��2009-Science-Deng Xinyan
% C_N2=sign(alpha2).*sqrt(C_L2.^2+C_D2.^2); 
% (4) Normal and tangential coefficients with Alpha�����湥�Ǳ仯�ķ��������������ϵ�����������֡���ȫ��ϵ��
% C_N3=3.4*sin(abs(alpha));      % 2006-IEEE TR-Deng Xinyan �� 2010-EAE-RJ Wood-����������ϵ�����
% ����(tangential)������ϵ��C_Tֻ��alpha��(-pi/4,pi/4)ʱ��Ϊ�㣬�������Ϊ��
% if alpha<pi/4 & alpha>-pi/4
%     C_T=0.4*cos(2*alpha).^2
% else
%     C_T=0;
% end
% C_T=0.4*(cos(2*alpha)).^2.*(alpha>-pi/4 & alpha<pi/4);       %����������ϵ�����
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
% wing_m_output=zeros(200,12);
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_N1',C_L',C_D',C_T'];
wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_N1'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




