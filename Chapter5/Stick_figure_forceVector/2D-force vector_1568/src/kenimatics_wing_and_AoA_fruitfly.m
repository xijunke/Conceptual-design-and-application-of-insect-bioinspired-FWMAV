function wing_m_output=kenimatics_wing_and_AoA_fruitfly()
% %% kenimatics_wing_and_AoA
% % ���˶����ɺͼ��ι���(AOA)
% % (ע�⣺���ռ��������Ťת��,��Ҫʱ���������ʵ��ʱ������׼)��������Ҫ������˶�
% % clear all;clc;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% �������ź������󵼡�����ý�����
% % (1) 2014-Science-ʵ���Ťת�ǡ��� ����4�׸���Ҷ������� 
% syms t         % t  �Լ�����ʱ����
% w =1185.6;     % ��Ƶ��   %  f=188.7; T=1/f;  %����Ƶ�� (Hz)������  
% phi1=sym('(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180');  %���ź���
% dphi1=diff(phi1,t,1);
% ddphi1=diff(phi1,t,2); 
% dphi=inline(vectorize(dphi1),'w','t');                    % ��ֵ����
% ddphi=inline(vectorize(ddphi1),'w','t');   
% % (2) 2014-Science-ʵ���Ťת�ǡ�������8�׸���Ҷ�������    
% psi1=sym('(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180');  %���ź���
% dpsi1=diff(psi1,t,1);   
% ddpsi1=diff(psi1,t,2);    
% dpsi=inline(vectorize(dpsi1),'w','t');  %��ֵ����
% ddpsi=inline(vectorize(ddpsi1),'w','t');  %��ֵ����
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ������������ֵ
% f=188.7; T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
% % t=linspace(0.0028,0.0028+2*T,200);  % t_steady1                                                                %���
% t=linspace(0.0028,0.0028+T,200);  % t_steady1  
% dphi=dphi(w,t);               %(1*200)�������������Ĵ������                                                     %���
% ddphi=ddphi(w,t);          %(1*200)�������������Ĵ�Ǽ�����                                                  %���
% dpsi=dpsi(w,t);          %(1*200)������������Ťת������                                                           %���
% ddpsi=ddpsi(w,t);     %(1*200)������������Ťת�Ǽ�����                                                        %���
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% �Ĵ�Ǻ�Ťת��
% % �Ĵ�ǡ��� (1*200)�������������Ĵ�ǡ���������            %���                                     
% phi=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+...
%         0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180;
% % Ťת�ǡ���(1*200)������������Ťת�ǡ���������            %���
% psi=(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+... 
%         3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+...
%         0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180;   
% % [phi_min,k]=min(phi);  % ���:  phi_min =-1.0157;  k =10756;      
% % t_0=t(k);                       % ���: t0 =0.0028;  
% figure(1)                                                          % ͼ1�����Ĵ�Ǻ�Ťת��
% plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-')  %ת��Ϊms �� ����degree   *10^3   *180/pi
% xlabel('\itNormalized time')
% ylabel('\itAngle (��)')
% legend('\it\phi(t)','\it\psi(t)')
% title('�Ĵ�Ǻ�Ťת����ʱ��ı仯����')   % �Ĵ�Ǻ�Ťת����ʱ��ı仯����
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) ����ȡ�ĵ����ȶ������ڡ����Թ��Ǳ仯����Ϊ������˶�ѧ���ݶ��롪��������
wingmotion_oneT=xlsread('wingmotion_oneT.xlsx','A1:G2000');
% f=188.7;  %T=1/f;  %����Ƶ�� (Hz)������  % w =1185.6; 
t=wingmotion_oneT(:,1);             % ��λ��ms
psi=wingmotion_oneT(:,2);            % Ťת�ǡ�����λ��rad
dpsi=wingmotion_oneT(:,3);            % Ťת�ǡ�����λ��rad
ddpsi=wingmotion_oneT(:,4);             % Ťת�ǡ�����λ��rad
phi=wingmotion_oneT(:,5);                   % �Ĵ�ǡ�����λ��rad
dphi=wingmotion_oneT(:,6);                     % �Ĵ�ǡ�����λ��rad
ddphi=wingmotion_oneT(:,7);                      % �Ĵ�ǡ�����λ��rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-Ťת�Ǻͼ��ι���AOA
alpha1=pi/2+psi.*sign(dphi);        %(1*200)���������������ι��ǡ���������        %�������ȫ�����ι��� 
dalpha1=dpsi.*sign(dphi);             %(1*200)���������������ι���-������             %������������и�
ddalpha1=ddpsi.*sign(dphi);         %(1*200)���������������ι���-������             %������������и�
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
% plot(t/T,psi*180/pi,'b-',t/T,alpha1*180/pi,'g-')      %Ťת�Ǻͼ��ι���AOA��ʱ��ı仯����
% xlabel('\itNormalized time')
% ylabel('\itAngle (��)')
% legend('\it\psi(t)','\it\alpha_1(t)')
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-');     %��x-axis
% title('Ťת�Ǻͼ��ι���alpha_1��ʱ��ı仯����') 
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ��ٶ�
% omega_x=-dpsi;   % Ťת��:   x(1)=psi;  x(2)=dpsi;   
omega_y=dphi.*sin(psi); 
omega_z=dphi.*cos(psi);
% omega_h=dphi;     % �������ϵ�½���������    % omega_h=sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
%% �������ϵ�µĽǼ��ٶȡ����������������ļ���
% domega_x=-ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi); 
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ǵ����
% ����alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
alpha2=atan2(omega_z,-omega_y);      % �������: �������ι��ǡ����ᳫʹ��XXX���ι��ǡ�����ȡ����ֵ����ȷ
% figure(3)                                               % ͼ3����ʹ�ó������ϵ�µĽ��ٶ���⼸�ι���
% plot(t/T,alpha1*180/pi,'r-',t/T,alpha2*180/pi,'b-')  
% % plot(t/T,alpha*180/pi,'r-',t/T,abs(alpha2)*180/pi,'b-')  
% xlabel('\itNormalized time')
% ylabel('\it\alpha_1 & \alpha_2 (deg)')
% legend('\alpha_1 (t)','\alpha_2 (t)')
% title('������ʱ��ı仯����')   % ������ʱ��ı仯����
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-�Ĵ�����ʺ�Ťת������
%% �Ĵ������dphi��Ťת���ǽ�����dalpha��ǰ������
% figure(4)                                                              % ͼ4�����Ĵ�����ʺ�Ťת������
% % plot(t/T,abs(dphi),'r-',t/T,abs(dalpha1),'b-')           % �Ĵ�����ʺ�Ťת��������ʱ��ı仯����
% plot(t/T,dphi,'r-',t/T,dalpha1,'b-')  
% xlabel('\itNormalized time')
% ylabel('������ (rad/s)')
% legend('\itd\phi(t)','\itd\alpha_1(t)')
% title('�Ĵ�����ʺ�Ťת��������ʱ��ı仯����') 
% grid on
%% Plot-�Ĵ�Ǽ����ʣ�Ťת�Ǽ����ʣ����ǽǼ�����
%�Ĵ�Ǽ�����ddphi��Ťת���ǽǼ�����ddalpha1��ǰ������
% figure(5)        % ͼ5���Ĵ�Ǽ�����, Ťת�Ǽ�����, ���ǽǼ�������ʱ��ı仯����
% plot(t/T,abs(ddphi),'r-',t/T,ddpsi,'b-',t/T,ddalpha1,'g-')
% % plot(t/T,abs(ddphi),'r-',t/T,abs(ddpsi),'b-',t/T,abs(ddalpha1),'g-') 
% xlabel('\itNormalized time')
% ylabel('�Ǽ����� (rad/s^-2)')
% legend('\itdd\phi(t)','\itdd\psi(t)','\itdd\alpha_1(t)')
% title('�Ĵ�Ǽ����ʣ�Ťת�Ǽ����ʣ����ǽǼ�������ʱ��ı仯����') 
% grid on
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
C_L =sign(alpha).*(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180));  %�Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�Ŷ
C_D =sign(alpha).*(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));  %�Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�Ŷ��1999-Science-MH Dickinson
% C_D =1.92+1.55*cos((2.04*alpha1-9.82)*pi/180);      %�Ƕ�alpha�Զ�����ʾ,�������Ǻ�������Ի����Ƶ�Ŷ��2004-JEB-Wang ZJ
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
% % figure(6)        % ͼ6��������ϵ�¡�������ϵ��
% % C_LD=plot(t/T,C_L,'r-',t/T,C_D,'b-');      %ʱ�����򻯣�����Ƶ��f, ���߳�������T;
% % xlabel('\itt (Normalized time with flapping period)')
% % ylabel('\itForce coefficients')
% % % title('Coefficients of lift and drag \itvs. t \rm for flapping wing')
% % title('����������ϵ������ʱ��ı仯����')
% % legend('\itC_L','\itC_D')
% % set(C_LD,'LineWidth',2)
% % grid on
% figure(7)     % ͼ7��������ϵ�¡����������������ϵ��
% plot(t/T,C_N3,'r-',t/T,C_T,'b-')                                % �õ��ķ��������������ϵ����Ҫ��ʵ
% xlabel('Normalized time')
% ylabel('C_N3(\alpha(t)) & C_T(\alpha (t))')
% legend('C_{N3}(\alpha)','C_T(\alpha)')
% title('���������������ϵ����ʱ��ı仯����')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %���ַ���������ϵ��
% figure(8)     % ͼ8��������ϵ�¡�����������ϵ��
% C_N=plot(t/T,C_N1,'r-',t/T,C_N2,'b-',t/T,C_N3,'g-');      %ʱ�����򻯣�����Ƶ��f, ���߳�������T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of normal force \itvs. t \rm for flapping wing')
% title('����������ϵ����ʱ��ı仯����')
% legend('\itC_{N1}','\itC_{N2}','\itC_{N3}')
% set(C_N,'LineWidth',2)
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���˶����������������ϵ�����
% wing_m_output=zeros(200,17);
% wing_m_output=[t',phi',psi',alpha1',alpha2',dphi',dpsi',dalpha1',ddphi',ddpsi',ddalpha1',C_L',C_D',C_N1',C_N2',C_N3',C_T'];
wing_m_output=[t,phi,psi,alpha1,alpha2,dphi,dpsi,dalpha1,ddphi,ddpsi,ddalpha1,C_L,C_D,C_N1,C_N2,C_N3,C_T];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




