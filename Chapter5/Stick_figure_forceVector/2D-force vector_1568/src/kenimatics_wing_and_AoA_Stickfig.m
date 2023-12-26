function wing_m_output=kenimatics_wing_and_AoA_Stickfig()
%% kenimatics_wing_and_AoA
% ��Դ�Ʒ䡪��Bumblebee
% clear all; clc;
syms   phi_m  K      theta_m   Phi_theta   theta_0   N        eta_m   C_eta   Phi_eta    eta_0      f  t  
%% ��̽Ǻͽ��ٶȡ��Լ��Ǽ��ٶ�
phi=sym('phi_m*asin(K*sin(2*pi*f*t))/asin(K)');          % ���ź���
dphi1=diff(phi,t,1);
ddphi1=diff(phi,t,2);
dphi=inline(vectorize(dphi1),'phi_m','K','f','t');            % ��ֵ����
ddphi=inline(vectorize(ddphi1),'phi_m','K','f','t');
%% �ڶ��Ǻͽ��ٶȡ��Լ��Ǽ��ٶ�
theta=sym('theta_m*cos(2*pi*N*f*t+Phi_theta)+theta_0');        %���ź���
dtheta1=diff(theta,t,1);
ddtheta1=diff(theta,t,2);
dtheta=inline(vectorize(dtheta1),'theta_m','N','Phi_theta','theta_0','f','t');            %��ֵ����
ddtheta=inline(vectorize(ddtheta1),'theta_m','N','Phi_theta','theta_0','f','t');
%% Ťת�Ǻͽ��ٶȡ��Լ��Ǽ��ٶ�
% eta=eta_m*tanh(C_eta*sin(2*pi*f*ti+Phi_eta))./tanh(C_eta);
% eta=eta_m*tanh(C_eta*sin(2*pi*f*ti+Phi_eta))./tanh(C_eta)+pi/2;
eta=sym('eta_m*tanh(C_eta*sin(2*pi*f*t+Phi_eta))/tanh(C_eta)+eta_0');        %���ź���
deta1=diff(eta,t,1);
ddeta1=diff(eta,t,2);
deta=inline(vectorize(deta1),'eta_m','C_eta','Phi_eta','eta_0','f','t');            %��ֵ����
ddeta=inline(vectorize(ddeta1),'eta_m','C_eta','Phi_eta','eta_0','f','t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ������������ֵ
f=122;  T=1/f;                        % Hz
phi_m=pi/2;                          % rad.
theta_m=12.3*pi/180;          % rad.
eta_m=87*pi/180;                % rad.
theta_0=1.83*pi/180;           % rad.
eta_0=-91.8*pi/180;             % rad.
K=0.925; 
C_eta=1.223;
N=2;
Phi_theta=-102.2*pi/180;     % rad.
Phi_eta=-91.8*pi/180;          % rad.
%% ���˶���
t=linspace(0.25*T,1.25*T,32);  %���
phi=phi_m*asin(K*sin(2*pi*f*t))/asin(K);
theta=theta_m*cos(2*pi*N*f*t+Phi_theta)+theta_0;
eta=eta_m*tanh(C_eta*sin(2*pi*f*t+Phi_eta))/tanh(C_eta)+eta_0+pi/2;   % ע������+pi/2
%% ���˶����ٶȺͽǼ��ٶ�
dphi=dphi(phi_m,K,f,t);        %(1*100)�������������Ĵ������                                                    %���
ddphi=ddphi(phi_m,K,f,t);    %(1*100)�������������Ĵ�Ǽ�����                                                %���
dtheta=dtheta(theta_m,N,Phi_theta,theta_0,f,t);          %(1*100)������������Ťת������             %���
ddtheta=ddtheta(theta_m,N,Phi_theta,theta_0,f,t);     %(1*100)������������Ťת�Ǽ�����           %���
deta=deta(eta_m,C_eta,Phi_eta,eta_0,f,t);                %(1*100)������������Ťת������                  %���
ddeta=ddeta(eta_m,C_eta,Phi_eta,eta_0,f,t);            %(1*100)������������Ťת�Ǽ�����              %���

%% ���������˶����ɡ������ӻ�
figure(1)
plot(t*f,phi*180/pi,'r-',t*f,theta*180/pi,'g-',t*f,eta*180/pi,'b-');
grid on
xlabel('������ʱ��')
ylabel('��̽� (\phi(t)), �ڶ��� (\theta(t)), Ťת�� (\eta(t))  (deg.)')
title('���������˶�����')
legend('\phi(t)','\theta(t)','\eta(t)')

%% ���ٶ�
% v_x(r)=r.*(dphi(t).*cos(theta(t)).*cos(eta(t))+dtheta(t).*sin(eta(t)));
% v_y(r)=r.*(dtheta(t).*cos(eta(t))-dphi(t).*cos(theta(t)).*sin(eta(t)));
v_x=dphi.*cos(theta).*cos(eta)+dtheta.*sin(eta);
v_y=dtheta.*cos(eta)-dphi.*cos(theta).*sin(eta);
%% �߼��ٶ�
% a_x(r)=r.*((ddphi(t).*cos(theta(t))+dtheta(t).*(deta(t)-dphi(t).*sin(theta(t)))).*cos(eta(t))+(ddtheta(t)-deta(t).*dphi(t).*cos(theta(t))).*sin(eta(t)));
% a_y(r)=r.*((-ddphi(t).*cos(theta(t))+dtheta(t).*(deta(t)-dphi(t).*sin(theta(t)))).*sin(eta(t))+(ddtheta(t)-deta(t).*dphi(t).*cos(theta(t))).*cos(eta(t)));
a_x=(ddphi.*cos(theta)+dtheta.*(deta-dphi.*sin(theta))).*cos(eta)+(ddtheta-deta.*dphi.*cos(theta)).*sin(eta);
a_y=(-ddphi.*cos(theta)+dtheta.*(deta-dphi.*sin(theta))).*sin(eta)+(ddtheta-deta.*dphi.*cos(theta)).*cos(eta);
%% ���ι���
% alpha=atan(sign(v_x(R/2)./v_y(R/2)).*v_x(R/2)./v_y(R/2));
% alpha=atan(v_x(R/2)./v_y(R/2));
% alpha=atan2(v0_x(R/2),v0_y(R/2));   %
% ȡ����1/2�᳤�������ٶ��������λ�õĹ��ǡ�������Ŷ���ڳ�����ϵ������ͷ�������ٶȵı�ֵ�໥������r;
alpha=atan2(v_x,v_y);
%% ������ϵ�����ٶȿ��ӻ�
figure(2)
subplot(211)
% plot(t,v_x(R/2),'r-',t*f,v_y(R/2),'b-');
plot(t*f,v_x,'r-',t*f,v_y,'b-');
grid on
xlabel('������ʱ��')
ylabel('v_x(t) �� v_y(t) (m/s)')
title('�������µ����ٶȱ仯����')
legend('v_x(t)','v_y(t)')
%% ���ι��ǵĿ��ӻ�
subplot(212)
plot(t*f,alpha*180/pi,'r-');  
grid on
xlabel('������ʱ��')
ylabel('���ι��� (\alpha (t)) (deg.)')
title('�������¼��ι����˶�����')
legend('\alpha(t)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���˶����������������ϵ�����
wing_m_output=zeros(1000,15);
wing_m_output=[t',phi',theta',eta',alpha',dphi',dtheta',deta',ddphi',ddtheta',ddeta',v_x',v_y',a_x',a_y'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
