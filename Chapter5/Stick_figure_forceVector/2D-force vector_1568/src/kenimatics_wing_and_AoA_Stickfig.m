function wing_m_output=kenimatics_wing_and_AoA_Stickfig()
%% kenimatics_wing_and_AoA
% 针对大黄蜂——Bumblebee
% clear all; clc;
syms   phi_m  K      theta_m   Phi_theta   theta_0   N        eta_m   C_eta   Phi_eta    eta_0      f  t  
%% 冲程角和角速度、以及角加速度
phi=sym('phi_m*asin(K*sin(2*pi*f*t))/asin(K)');          % 符号函数
dphi1=diff(phi,t,1);
ddphi1=diff(phi,t,2);
dphi=inline(vectorize(dphi1),'phi_m','K','f','t');            % 数值函数
ddphi=inline(vectorize(ddphi1),'phi_m','K','f','t');
%% 摆动角和角速度、以及角加速度
theta=sym('theta_m*cos(2*pi*N*f*t+Phi_theta)+theta_0');        %符号函数
dtheta1=diff(theta,t,1);
ddtheta1=diff(theta,t,2);
dtheta=inline(vectorize(dtheta1),'theta_m','N','Phi_theta','theta_0','f','t');            %数值函数
ddtheta=inline(vectorize(ddtheta1),'theta_m','N','Phi_theta','theta_0','f','t');
%% 扭转角和角速度、以及角加速度
% eta=eta_m*tanh(C_eta*sin(2*pi*f*ti+Phi_eta))./tanh(C_eta);
% eta=eta_m*tanh(C_eta*sin(2*pi*f*ti+Phi_eta))./tanh(C_eta)+pi/2;
eta=sym('eta_m*tanh(C_eta*sin(2*pi*f*t+Phi_eta))/tanh(C_eta)+eta_0');        %符号函数
deta1=diff(eta,t,1);
ddeta1=diff(eta,t,2);
deta=inline(vectorize(deta1),'eta_m','C_eta','Phi_eta','eta_0','f','t');            %数值函数
ddeta=inline(vectorize(ddeta1),'eta_m','C_eta','Phi_eta','eta_0','f','t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 给各个函数赋值
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
%% 翅运动角
t=linspace(0.25*T,1.25*T,32);  %输出
phi=phi_m*asin(K*sin(2*pi*f*t))/asin(K);
theta=theta_m*cos(2*pi*N*f*t+Phi_theta)+theta_0;
eta=eta_m*tanh(C_eta*sin(2*pi*f*t+Phi_eta))/tanh(C_eta)+eta_0+pi/2;   % 注意这里+pi/2
%% 翅运动角速度和角加速度
dphi=dphi(phi_m,K,f,t);        %(1*100)的行向量——拍打角速率                                                    %输出
ddphi=ddphi(phi_m,K,f,t);    %(1*100)的行向量——拍打角加速率                                                %输出
dtheta=dtheta(theta_m,N,Phi_theta,theta_0,f,t);          %(1*100)的行向量——扭转角速率             %输出
ddtheta=ddtheta(theta_m,N,Phi_theta,theta_0,f,t);     %(1*100)的行向量——扭转角加速率           %输出
deta=deta(eta_m,C_eta,Phi_eta,eta_0,f,t);                %(1*100)的行向量——扭转角速率                  %输出
ddeta=ddeta(eta_m,C_eta,Phi_eta,eta_0,f,t);            %(1*100)的行向量——扭转角加速率              %输出

%% 参数化翅运动规律——可视化
figure(1)
plot(t*f,phi*180/pi,'r-',t*f,theta*180/pi,'g-',t*f,eta*180/pi,'b-');
grid on
xlabel('无量纲时间')
ylabel('冲程角 (\phi(t)), 摆动角 (\theta(t)), 扭转角 (\eta(t))  (deg.)')
title('参数化翅运动规律')
legend('\phi(t)','\theta(t)','\eta(t)')

%% 线速度
% v_x(r)=r.*(dphi(t).*cos(theta(t)).*cos(eta(t))+dtheta(t).*sin(eta(t)));
% v_y(r)=r.*(dtheta(t).*cos(eta(t))-dphi(t).*cos(theta(t)).*sin(eta(t)));
v_x=dphi.*cos(theta).*cos(eta)+dtheta.*sin(eta);
v_y=dtheta.*cos(eta)-dphi.*cos(theta).*sin(eta);
%% 线加速度
% a_x(r)=r.*((ddphi(t).*cos(theta(t))+dtheta(t).*(deta(t)-dphi(t).*sin(theta(t)))).*cos(eta(t))+(ddtheta(t)-deta(t).*dphi(t).*cos(theta(t))).*sin(eta(t)));
% a_y(r)=r.*((-ddphi(t).*cos(theta(t))+dtheta(t).*(deta(t)-dphi(t).*sin(theta(t)))).*sin(eta(t))+(ddtheta(t)-deta(t).*dphi(t).*cos(theta(t))).*cos(eta(t)));
a_x=(ddphi.*cos(theta)+dtheta.*(deta-dphi.*sin(theta))).*cos(eta)+(ddtheta-deta.*dphi.*cos(theta)).*sin(eta);
a_y=(-ddphi.*cos(theta)+dtheta.*(deta-dphi.*sin(theta))).*sin(eta)+(ddtheta-deta.*dphi.*cos(theta)).*cos(eta);
%% 几何攻角
% alpha=atan(sign(v_x(R/2)./v_y(R/2)).*v_x(R/2)./v_y(R/2));
% alpha=atan(v_x(R/2)./v_y(R/2));
% alpha=atan2(v0_x(R/2),v0_y(R/2));   %
% 取得是1/2翅长处的线速度来计算该位置的攻角——不对哦，在翅坐标系下弦向和法向的线速度的比值相互抵消了r;
alpha=atan2(v_x,v_y);
%% 翅坐标系的线速度可视化
figure(2)
subplot(211)
% plot(t,v_x(R/2),'r-',t*f,v_y(R/2),'b-');
plot(t*f,v_x,'r-',t*f,v_y,'b-');
grid on
xlabel('无量纲时间')
ylabel('v_x(t) 和 v_y(t) (m/s)')
title('翅坐标下的线速度变化规律')
legend('v_x(t)','v_y(t)')
%% 几何攻角的可视化
subplot(212)
plot(t*f,alpha*180/pi,'r-');  
grid on
xlabel('无量纲时间')
ylabel('几何攻角 (\alpha (t)) (deg.)')
title('翅坐标下几何攻角运动规律')
legend('\alpha(t)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅运动规律输出和气动力系数输出
wing_m_output=zeros(1000,15);
wing_m_output=[t',phi',theta',eta',alpha',dphi',dtheta',deta',ddphi',ddtheta',ddeta',v_x',v_y',a_x',a_y'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
