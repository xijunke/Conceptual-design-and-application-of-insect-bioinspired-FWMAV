% function wing_m_output=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,K,eta_m,C_eta,Phi_eta) % ,eta_0
% function wing_m_output=kenimatics_wing_and_AoA_fruitfly_sim(f,phi_m,epsilon,psi_m,zeta,psi_0)
%% kenimatics_wing_and_AoA
% 翅运动规律和几何攻角(AOA)
% (注意：最终计算出来的扭转角,需要时间轴区间和实验时间轴配准)———需要求解和核对
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%% %翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% x =[46.6104,1.3076,0.1891,1.2715,2.5103,-1.5169,3.9813,1.96,1.7845,0.0001]; % P_asterisk =4.7867;L =1.0000; delta = -2.9524e-008; AR =2.9417; Re =186.1461;% 下一行为精确值
x =[4.6610447e+001,1.3075639e+0,1.8911207e-001,1.2714657e+0,2.5103035e+0, -1.5168925e+0,3.9812634e+0,1.9600006e+0,1.7844515e+0,7.9074669e-005];
f=x(1);  T=1/f;   % Hz 
phi_m=x(2);
K=x(3);                  %  epsilon=x(3);    
eta_m=x(4);          % phi_0=x(4);       
C_eta=x(5);           %  psi_m=x(5);
Phi_eta=x(6);        % zeta=x(6); 
% eta_0=x(7);       % psi_0=x(7);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 创建符号函数，求导——获得角速率—————用于计算气动力和力矩时，需要更新计算而得的扭转角
% % (1) 2014-Science-实测翅拍打角—— 按照4阶傅里叶级数拟合 
% syms t         % t  自己设置时间轴
% w =1185.6;     % 角频率   %  f=188.7; T=1/f;  %翅拍频率 (Hz)和周期  
% phi1=sym('(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180');  %符号函数: 原始数据的度数——转化为弧度制
% dphi1=diff(phi1,t,1);
% ddphi1=diff(phi1,t,2); 
% dphi=inline(vectorize(dphi1),'w','t');                    % 数值函数
% ddphi=inline(vectorize(ddphi1),'w','t');   
% % (2) 2014-Science-实测翅扭转角——按照8阶傅里叶级数拟合    
% psi1=sym('-(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180');  %符号函数:原始数据的度数——转化为弧度制
% dpsi1=diff(psi1,t,1);   
% ddpsi1=diff(psi1,t,2);    
% dpsi=inline(vectorize(dpsi1),'w','t');  %数值函数
% ddpsi=inline(vectorize(ddpsi1),'w','t');  %数值函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 给各个函数赋值
% 七个变量——(phi_m,K,eta_m,C_eta,Phi_eta,eta_0,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 给各个函数赋值——雄蜂——最优的翅膀运动学参数
% f=122;                            % Hz——最优值
% phi_m=pi/2;                   % rad.——最优值
% eta_m=87.0*pi/180;       % rad.——最优值
% % eta_0=-pi/2;                  % rad.——最优值
% eta_0=0;                        % rad.——————初始值符合
% K=0.925;                        %——最优值
% C_eta=1.223;                  %——最优值
% Phi_eta=-91.8*pi/180;    % rad.——最优值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 给各个函数赋值——果蝇——最优的翅膀运动学参数
% f=234;  % Hz——f∈[0,inf];     % T=1/f;                            %—Frequency—最优值            
% phi_m=pi/2;                % rad.——phi_m∈[0,pi/2];            %—Azimuthal amplitude—最优值
% eta_m=72.7*pi/180;    % rad.——eta_m∈[0,pi];                %—Pitching amplitude—最优值
% eta_0=pi/2;             % rad.——eta_0∈[eta_m-pi,pi-eta_m];     %—Pitching offset—最优值
% % eta_0=0;                     % rad.——————————————————Pitching offset———初始值符合
% K=0.704;                     % K∈[0,1]; % —Affects the shape of phi(t )—最优值
% C_eta=2.375;               % C_eta∈[0,inf];                           %—Affects the duration of wing rotation—最优值
% Phi_eta=-72.4*pi/180;  % rad.——Phi_eta∈[-pi,pi];          %—Pitching phase offset—最优值
% % P_asterisk=49.6805; L=3.2297;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eta_m=0;                       % rad.—Pitching amplitude—Min
% eta_m=pi/2;  
% eta_m=pi;                     % rad.—Pitching amplitude—Max
% %%%%%%%%%%%%%%
% eta_0=eta_m-pi;            % rad.—Pitching offset—Min
% eta_0=pi-eta_m;            % rad.—Pitching offset—Max
% % %%%%%%%%%%%%%%
% K=0.0001;       % K∈[0,1];——K=0时phi为正弦曲线; K=1时phi为三角方波; %K可被看作是翅膀多快地逆反其方向的度量，即逆反速度的度量？
% C_eta=0.0001; % C_eta∈[0,inf];——C_eta=0时eta为正弦曲线; C_eta=inf时eta倾向为阶跃函数; % C_eta的值与翅膀俯仰逆反的时长呈负相关？？？
% Phi_eta=-pi/2;     % rad.——Phi_eta∈[-pi,pi];——Pitching phase offset——比较好
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w =1185.6;     % 角频率   %  f=188.7; T=1/f;  %翅拍频率 (Hz)和周期 
% f=188.7; 
T=1/f;  %翅拍频率 (Hz)和周期  % w =1185.6; 
% t_00=solve(dphi(phi_m,K,f,t)=0,t); 
t_00= 1/(4*f);    % t_00= -1/(4*f); % for Wang ZJ wingbeat pattern
% t_00= -epsilon/(2*pi*f)+T;        % for 谐波函数wingbeat  pattern
t=linspace(t_00,t_00+3*T,1000); 
% t=linspace(0.0052824335,0.0052824335+3*T,100000); % t_steady1——针对2014-Science-实测翅拍打角和扭转角
% t=linspace(0.0102,0.0102+3*T,100000);                           % t_steady1——针对2007-JFM-Wang ZJ-人为设计运动学规律——雄蜂
% t=linspace(0.0139,0.0139+3*T,1000);                           % t_steady1——针对2007-JFM-Wang ZJ-人为设计运动学规律——果蝇
% dphi=dphi(w,t);         %(1*100)的行向量——拍打角速率                                                           % 输出
% ddphi=ddphi(w,t);    %(1*100)的行向量——拍打角加速率                                                         % 输出
% dpsi=dpsi(w,t);         %(1*100)的行向量——扭转角速率                                                             % 输出
% ddpsi=ddpsi(w,t);     %(1*100)的行向量——扭转角加速率                                                         % 输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 拍打角和扭转角
% % 拍打角—— (1*100)的行向量——拍打角——弧度制            %输出                                     
% phi=(3.9008+65.0445*cos(w*t)+4.2642*sin(w*t)+3.5806*cos(2*w*t)-2.9492*sin(2*w*t)+...
%         0.1319*cos(3*w*t)+0.3639*sin(3*w*t)+0.7844*cos(4*w*t)+0.2098*sin(4*w*t))*pi/180; % 原始数据的度数——转化为弧度制
% % 扭转角——(1*100)的行向量——扭转角——弧度制            %输出
% psi=-(4.2936+1.8536*cos(t*w)+59.6529*sin(t*w)+5.1852*cos(2*t*w)+1.6095*sin(2*t*w)-8.4569*cos(3*t*w)+9.8057*sin(3*t*w)+... 
%         3.9562*cos(4*t*w)+5.8064*sin(4*t*w)-3.0337*cos(5*t*w)-2.8749*sin(5*t*w)-2.8771*cos(6*t*w)+...
%         0.6686*sin(6*t*w)+0.8649*cos(7*t*w)-0.6137*sin(7*t*w)+0.0771*cos(8*t*w)-1.0007*sin(8*t*w))*pi/180;  % 原始数据的度数——转化为弧度制
% % [phi_min,k]=min(phi);  % 输出:  phi_min =-1.0157;  k =10756;      
% % t_0=t(k);                       % 输出: t0 =0.0028;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% syms  phi_m   K        eta_m  C_eta  Phi_eta  eta_0          f    t
% phi=phi_m.*asin(K.*sin(2.*pi.*f.*t))./asin(K);
% dphi=diff(phi,t,1);
% ddphi=diff(phi,t,2);
% eta=eta_m.*tanh(C_eta.*sin(2.*pi.*f.*t+Phi_eta))./tanh(C_eta)+eta_0;
% deta=diff(eta,t,1);
% ddeta=diff(eta,t,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第二方案————7个变量——(f,phi_m,K,eta_m,C_eta,Phi_eta,eta_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi=phi_m.*asin(K.*sin(2.*pi.*f.*t))./asin(K);    % 拍打角、角速度和角加速度
dphi=(2.*K.*pi.*f.*phi_m.*cos(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(1./2));
ddphi=(4.*K.^3.*pi.^2.*f.^2.*phi_m.*cos(2.*pi.*f.*t).^2.*sin(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(3./2))...
            - (4.*K.*pi.^2.*f.^2.*phi_m.*sin(2.*pi.*f.*t))./(asin(K).*(1 - K.^2.*sin(2.*f.*pi.*t).^2).^(1./2));
% %%%%%%%%%%%%%%%
% % [phimax,index]=max(phi);
% % phi_max =phimax*180/pi
% % t_0=t(index)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psi=-(eta_m.*tanh(C_eta.*sin(2.*pi.*f.*t+Phi_eta))./tanh(C_eta)+eta_0); 
psi=-(eta_m.*tanh(C_eta.*sin(2.*pi.*f.*t+Phi_eta))./tanh(C_eta)); % 扭转角、角速度和角加速度——原始公式没有负号;
dpsi=(2.*C_eta.*pi.*eta_m.*f.*cos(Phi_eta + 2.*pi.*f.*t).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta);% 扭转角速度
ddpsi=- (4.*C_eta.*pi.^2.*eta_m.*f.^2.*sin(Phi_eta + 2.*pi.*f.*t).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta)...
      - (8.*C_eta.^2.*pi.^2.*eta_m.*f.^2.*cos(Phi_eta + 2.*pi.*f.*t).^2.*tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).*(tanh(C_eta.*sin(Phi_eta + 2.*pi.*f.*t)).^2 - 1))./tanh(C_eta);% 扭转角加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 第四方案———— 人为设计谐波运动(harmonic_motion)——6个变量——(f,phi_m,epsilon,psi_m,zeta,psi_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phi=phi_m.*cos(2.*pi.*f.*t+epsilon);                        % 拍打角
% dphi=-2.*pi.*f.*phi_m.*sin(epsilon + 2.*pi.*f.*t);                 % 拍打角速度
% ddphi=-4.*pi.^2.*f.^2.*phi_m.*cos(epsilon + 2.*pi.*f.*t);   % 拍打角加速度
% psi=psi_m.*sin(2.*pi.*f.*t+zeta)+psi_0;                         % 扭转角
% dpsi =2.*pi.*f.*psi_m.*cos(zeta + 2.*pi.*f.*t);                % 扭转角速度
% ddpsi =-4.*pi.^2.*f.^2.*psi_m.*sin(zeta + 2.*pi.*f.*t);   % 扭转角加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi_min,k]=min(phi); 
t_0=t(k); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)       % 前后冲程背景图
hold on
% subplot(313)
% % plot(t_0/T,phi_min*180/pi,'dr','LineWidth',2)
% % hold on
%%%%%%%%%%%%%%%%%%%%%%%%%
juxingx1=[t(1,1)/T,(t_0)/T,(t_0)/T,t(1,1)/T];
juxingy1=[-105,-105,105,105];
% fill(juxingx1,juxingy1,'y');
fill(juxingx1,juxingy1,[0.95,0.95,0.95]);
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%
juxingx2=[(t_0)/T,(t(1,1)+T)/T,(t(1,1)+T)/T,(t_0)/T];
juxingy2=[-105,-105,105,105];
% fill(juxingx2,juxingy2,'c');
fill(juxingx2,juxingy2,[0.7,0.7,0.7]);
hold on
axis([min(t/T),min(t/T)+1,-105,105])
% set(gca,'XTick','off','YGrid','off')
set(gca,'XTick',(0:0:0),'YTick',(0:0:0))
%%%%%%%%%%%%%%%%%%%%%%%%%
% plot([(t_0)/T,(t_0)/T],[-105,105],'k-.','LineWidth',2)
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)       % 图1——拍打角和扭转角
plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-.','LineWidth',2)  %转换为ms 和 度数degree   *10^3   *180/pi
% plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-.','Color',[0,0,1],'LineWidth',2)
% plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:','LineWidth',2)  %转换为ms 和 度数degree   *10^3   *180/pi
% xlabel('\rmNormalized time','FontSize',20,'FontName','Times','FontWeight','Bold')
ylabel('\rmFlapping angle and pitch angle (°)','FontSize',20,'FontName','Times','FontWeight','Bold')
legend('\it\phi(t)','\it\psi(t)')
title('Flapping angle and pitch angle','FontSize',20,'FontWeight','Bold')   % 拍打角和扭转角随时间的变化规律
grid on
hold on
% plot([min(t/T),max(t/T)],[0,0],'k-.','LineWidth',2)
% axis([-0.45,2.6,-105,105]) 
% axis([min(t/T),max(t/T),-105,105])
axis([min(t/T),min(t/T)+1,-105,105])
% set(gca,'XTick',(0.9:0.1:4.05))
set(gca,'YGrid','off','XTick',(min(t/T):0.025:min(t/T)+1),'LineStyle','-','LineWidth',1.5,'FontSize',14,'FontName','Times','FontWeight','Bold') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)       % 图1——拍打角和扭转角
% subplot(311)
% plot(t/T,phi*180/pi,'r-',t/T,psi*180/pi,'b-','LineWidth',2)  %转换为ms 和 度数degree   *10^3   *180/pi
% xlabel('\itNormalized time')
% ylabel('\itAngle (°)')
% legend('\it\phi(t)','\it\psi(t)')
% title('拍打角和扭转角随时间的变化规律')   % 拍打角和扭转角随时间的变化规律
% grid on
% % axis([1.2444,4.2444,-105,105])            % ——雄蜂  t_0(0.0102)*f(122)=1.2444;
% % set(gca,'XTick',(1.2444:0.2:4.2444))     % ——雄蜂  t_0(0.0102)*f(122)=1.2444;
% % axis([3.2526,6.2526,-105,105])                % ——果蝇  t_0(0.0139)*f(234)=3.2526;
% % set(gca,'XTick',(3.2526:0.2:6.2526))         %——果蝇  t_0(0.0139)*f(234)=3.2526;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% t_phi_max=-epsilon/(2.*pi.*f)+T;                             % t_phi_max =-0.0109;
% t_psi_0=(asin(-psi_0/psi_m)-zeta)/(2.*pi.*f)+T;        % t_psi_0 =1.7425e-004;
% plot((t_phi_max+T/2)/T,-phi_m*180/pi,'rs',t_psi_0/T,0,'bs','LineWidth',3)
% hold on
% plot([t_psi_0/T,t_psi_0/T],[-100,100],'g-.',[(t_phi_max+T/2)/T,(t_phi_max+T/2)/T],[-100,100],'g-.','LineWidth',1.5)
% % hold on
% % plot((t_phi_max+T)/T,phi_m*180/pi,'rs',(t_psi_0+T/2)/T,0,'bs','LineWidth',4)
% delta_t1=t_phi_max+T-(t_psi_0+T/2);  % delta_t1 = 0.0018;
% delta_t=(pi-((epsilon+asin(-psi_0/psi_m)-zeta)))/(2.*pi.*f); % delta_t =0.0018;
% delta=pi-(epsilon+asin(-psi_0/psi_m)-zeta); % delta =0.4495;
% delta_deg=delta*180/pi;   % delta_deg =25.7563;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold on
% k=0;
% t_phi_0=(pi/2+2*k*pi-epsilon)/(2.*pi.*f)+T;   % t_phi_0 =-0.0044;
% psi_specific=psi_m.*sin(2.*pi.*f.*t_phi_0+zeta)+psi_0;  % psi_specific =-0.8767;
% alpha_specific=(pi/2+psi_specific)*180/pi;    % alpha_specific =39.7662;
% plot(t_phi_0/T,0,'rd',t_phi_0/T,psi_specific*180/pi,'bd','LineWidth',3)
% hold on
% plot([t_phi_0/T,t_phi_0/T],[0,psi_specific*180/pi],'k:','LineWidth',1.5)
% hold on
% t_phi_01=(pi/2+2*k*pi-epsilon)/(2.*pi.*f)+T/2+T;  % t_phi_01 =0.0085;
% psi_specific1=psi_m.*sin(2.*pi.*f.*t_phi_01+zeta)+psi_0;  % psi_specific1 =0.8797;
% alpha_specific1=(pi/2-psi_specific1)*180/pi;   % alpha_specific1 =39.5943;
% plot(t_phi_01/T,0,'rd',t_phi_01/T,psi_specific1*180/pi,'bd','LineWidth',4)
% hold on
% plot([t_phi_01/T,t_phi_01/T],[0,psi_specific1*180/pi],'k:','LineWidth',1.5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-拍打角速率和扭转角速率
% subplot(312)                                                       
% plot(t/T,dphi,'r-',t/T,dpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('角速率 (rad/s)')
% legend('\itd\phi(t)','\itd\psi(t)')
% title('拍打角速率和扭转角速率随时间的变化规律') 
% grid on
% % axis([1.2444,4.2444,-inf,inf]) % ——雄蜂
% % set(gca,'XTick',(1.2444:0.2:4.2444)) % ——雄蜂
% % axis([3.2526,6.2526,-inf,inf])  % ——果蝇
% % set(gca,'XTick',(3.2526:0.2:6.2526)) %——果蝇
% % Plot-拍打角加速率和扭转角角加速率
% subplot(313)
% plot(t/T,ddphi,'r-',t/T,ddpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('角加速率 (rad/s^-2)')
% legend('\itdd\phi(t)','\itdd\psi(t)')
% title('拍打角加速率和扭转角加速率随时间的变化规律') 
% grid on
% % axis([1.2444,4.2444,-inf,inf]) % ——雄蜂
% % set(gca,'XTick',(1.2444:0.2:4.2444)) % ——雄蜂
% % axis([3.2526,6.2526,-inf,inf])  % ——果蝇
% % set(gca,'XTick',(3.2526:0.2:6.2526)) %——果蝇
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-扭转角和几何攻角AOA
% % alpha1=pi/2+psi.*sign(dphi);        %(1*200)的行向量——几何攻角——弧度制        %输出——全正几何攻角 
% alpha1=pi/2-psi.*sign(dphi);        %(1*200)的行向量——几何攻角——弧度制        %输出——全正几何攻角  % 实测扭转角前负号
% % Y = sign(X) returns an array Y the same size as X, where each element of Y is:
% % *1 if the corresponding element of X is greater than zero
% % * 0 if the corresponding element of X equals zero
% % *-1 if the corresponding element of X is less than zero
% figure(4)                                                              % 图2—— 注意这里几何攻角始终取正值
% % % x_interval=[0,1/2,1/2,0];
% % % y_interval=[-100,-100,100,100];
% % % fill(x_interval,y_interval,'y');
% % % hold on
% % % legend('','\it\psi(t)','\it\alpha_1(t)')
% plot(t/T,psi*180/pi,'b-',t/T,alpha1*180/pi,'g-','LineWidth',2)      %扭转角和几何攻角AOA随时间的变化规律
% xlabel('\itNormalized time')
% ylabel('\itAngle (°)')
% legend('\it\psi(t)','\it\alpha_1(t)')
% hold on
% L=length(t);
% plot([0,t(L)/T],[0,0],'k-');     %画x-axis
% title('扭转角和几何攻角alpha_1随时间的变化规律') 
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀坐标系下的角速度
% omega_x=dpsi;   % 扭转角:   x(1)=psi;  x(2)=dpsi;   
omega_y=dphi.*sin(psi);
omega_z=dphi.*cos(psi);
% omega_h=dphi;    % 翅膀坐标系下铰链角速率    % omega_h=sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
%% 翅膀坐标系下的角加速度——用于虚质量力的计算
% domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi); 
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% 气动攻角计算——取流产速度相对于刚体好了速度，所以整体加负号；
v_y_nonr=-omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;    % v_z=-r*dphi*sin(psi)
warning('off')
alpha2=atan2(-v_y_nonr,v_z_nonr);   % 正确——注意与下文的alpha=atan2(omega_z,-omega_y)*180/pi; 不同  % 实测扭转角前负号
% 由于alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %这里atan2给出象限正负值，尤其是alpha>pi/2时
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(5)                                               % 图3——使用翅膀坐标系下的角速度求解几何攻角
% plot(t/T,alpha1*180/pi,'r-',t/T,alpha2*180/pi,'b-','LineWidth',2)    
% xlabel('\itNormalized time')
% ylabel('\it\alpha_1 & \alpha_2 (deg)')
% legend('\alpha_1 (t)','\alpha_2 (t)')
% title('攻角随时间的变化规律')   % 攻角随时间的变化规律
% grid on
% % axis([1.2444,4.2444,-inf,inf]) %——雄蜂
% % set(gca,'XTick',(1.2444:0.2:4.2444)) %——雄蜂
% % axis([3.2526,6.2526,-inf,inf])  % ——果蝇
% % set(gca,'XTick',(3.2526:0.2:6.2526)) %——果蝇
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅坐标系下：气动力系数
%% lift and drag coefficients with Alpha——随攻角变化的升阻力系数和法向气动力系数
% alpha=45;     %假定攻角恒定不变;下面公式的升阻力系数为:C_L=1.8046;C_D=1.7037.
alpha=alpha2;
% (1)下面的经验公式来自1999-science-MH Dickinson的文献——第一种
% C_L =0.225+1.58*sin((2.13*alpha*180/pi-7.2)*pi/180);  % 几何攻角正负交替(alpha),无需使用绝对值符号abs(alpha);
% C_D =1.92-1.55*cos((2.04*alpha*180/pi-9.82)*pi/180);  % 几何攻角正负交替(alpha),无需使用绝对值符号abs(alpha);
% C_N=cos(alpha).*C_L2+sin(alpha).*C_D2;  % 几何攻角正负交替(alpha2),无需使用绝对值符号abs(alpha);
% (2) lift and drag coefficients with Alpha——随攻角变化的升阻力系数和法向气动力系数—第二种—正负变化的系数
% (a) 1999-Science-MH Dickinson——角度alpha以度数表示,但是三角函数是针对弧度制的
C_L =(0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180)); 
C_D =(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));
% (b) 2004-JEB-Wang ZJ——角度alpha以度数表示,但是三角函数是针对弧度制的
% C_D =(1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180));
% C_D =1.92+1.55*cos((2.04*alpha1-9.82)*pi/180);
% (c) 平均升阻力系数和其比值
% C_L_aver=trapz(t,C_L)/(3*T); % C_L_aver =1.5088;
% C_D_aver=trapz(t,C_D)/(3*T) % C_D_aver =1.9123;
% C_L2C_D_aver=C_L_aver/C_D_aver; % C_L2C_D_aver =0.7890;
C_N1=cos(abs(alpha)).*C_L+sin(abs(alpha)).*C_D; %由升阻力系数合成——2010-JFM-RJ Wood
% % (3) lift and drag coefficients with Alpha——随攻角变化的升阻力系数和法向气动力系数——第三种——全正系数
% C_L2 =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  % 几何攻角全正abs(alpha);——升力系数全正
% C_D2 =1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180); % 几何攻角全正abs(alpha);——阻力系数全正
% C_N2=cos(abs(alpha)).*C_L2+sin(abs(alpha)).*C_D2;                % 几何攻角为正abs(alpha)——法向力系数全正
% % 正负交替-由升阻力系数合成—— 2012-IEEE ASME TM-Veaceslav Arabagi 或2009-Science-Deng Xinyan
% C_N2=sign(alpha2).*sqrt(C_L2.^2+C_D2.^2); 
% % (4) Normal and tangential coefficients with Alpha——随攻角变化的法向和切向气动力系数——第四种——全正系数
% C_N3=3.4*sin(abs(alpha));      % 2006-IEEE TR-Deng Xinyan 或 2010-EAE-RJ Wood-法向气动力系数输出
% 切向(tangential)气动力系数C_T只在alpha∈(-pi/4,pi/4)时不为零，其他情况为零
% if alpha<pi/4 & alpha>-pi/4
%     C_T=0.4*cos(2*alpha).^2
% else
%     C_T=0;
% end
% C_T=0.4*(cos(2*alpha)).^2.*(alpha>-pi/4 & alpha<pi/4);       %切向气动力系数输出
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(6)        % 图6—翅坐标系下—升阻力系数
% plot(t/T,C_L,'r-',t/T,C_D,'b-','LineWidth',2);      %时间正则化，乘以频率f, 或者除以周期T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of lift and drag \itvs. t \rm for flapping wing')
% title('气动升阻力系数随着时间的变化规律')
% legend('\itC_L','\itC_D')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7)
% plot(t/T,C_N1,'r-','LineWidth',2); 
% % axis([1.2444,4.2444,-inf,inf])  % ——雄蜂
% % set(gca,'XTick',(1.2444:0.2:4.2444))  % ——雄蜂
% % axis([3.2526,6.2526,-inf,inf])  % ——果蝇
% % set(gca,'XTick',(3.2526:0.2:6.2526)) %——果蝇
% figure(8)     % 图8—翅坐标系下—法向和切向气动力系数
% plot(t/T,C_N3,'r-',t/T,C_T,'b-','LineWidth',2)                                % 得到的法向和切向气动力系数需要核实
% xlabel('Normalized time')
% ylabel('C_N3(\alpha(t)) & C_T(\alpha (t))')
% legend('C_{N3}(\alpha)','C_T(\alpha)')
% title('法向和切向气动力系数随时间的变化规律')
% grid on
% axis([1.2444,4.2444,-inf,inf])  % ——雄蜂
% set(gca,'XTick',(1.2444:0.2:4.2444))  % ——雄蜂
% axis([3.2526,6.2526,-inf,inf])  % ——果蝇
% set(gca,'XTick',(3.2526:0.2:6.2526)) %——果蝇
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %三种法向气动力系数
% figure(9)     % 图9—翅坐标系下—法向气动力系数
% plot(t/T,C_N1,'r-',t/T,C_N2,'b-',t/T,C_N3,'g-','LineWidth',2);      %时间正则化，乘以频率f, 或者除以周期T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of normal force \itvs. t \rm for flapping wing')
% title('法向气动力系数随时间的变化规律')
% legend('\itC_{N1}','\itC_{N2}','\itC_{N3}')
% grid on
% axis([1.2444,4.2444,-inf,inf])            % ——雄蜂
% set(gca,'XTick',(1.2444:0.2:4.2444))  % ——雄蜂
% axis([3.2526,6.2526,-inf,inf])            % ——果蝇
% set(gca,'XTick',(3.2526:0.2:6.2526))  %——果蝇
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅运动规律输出和气动力系数输出
% wing_m_output=zeros(200,12);
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_N1',C_L',C_D',C_T'];
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_N1'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




