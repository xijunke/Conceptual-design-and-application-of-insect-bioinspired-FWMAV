function wing_m_output=kenimatics_wing_and_AoA_fruitfly_simpres()
%% kenimatics_wing_and_AoA
% 翅运动规律和几何攻角(AOA)
% (注意：最终计算出来的扭转角,需要时间轴区间和实验时间轴配准)———需要求解和核对
% clear all;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 创建符号函数，求导——获得角速率—————用于计算气动力和力矩时，需要更新计算而得的扭转角
%% (1) 人为设定拍打角函数
syms t         % t_range=[0.0052824335,0.0052824335+5*T];
w =1185.6;     % 角频率   %  f=188.7; T=1/f;  %翅拍频率 (Hz)和周期 
% f_var=188.7;  % Science实测拍打角四阶傅里叶级数函数的基频
% w =2*pi*f_var;   PHI=1.1487;
% phi_pres =PHI*cos(w*t);     
% 人为设定拍打角(幅值为PHI, 频率为Science实测拍打角四阶傅里叶级数函数的基频188.7Hz) % 人为设定拍打角: prescribed=规定的
phi_pres =sym('1.149*sin(w*t+1.571)');  %符号函数
dphi_pres =diff(phi_pres ,t,1);
ddphi_pres =diff(phi_pres ,t,2); 
dphi=inline(vectorize(dphi_pres ),'w','t');                    % 数值函数
ddphi=inline(vectorize(ddphi_pres ),'w','t');  
% dphi_pres =-w*PHI*sin(w*t);
% ddphi_pres =-w^2*PHI*cos(w*t);   
%% (2) 下面是求解模拟计算得到的扭转角和角速度以及角加速度    
psi_sim=sym('-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)-0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180');  %符号函数
dpsi_sim=diff(psi_sim,t,1);   
ddpsi_sim=diff(psi_sim,t,2);    
dpsi=inline(vectorize(dpsi_sim),'w','t');        %数值函数
ddpsi=inline(vectorize(ddpsi_sim),'w','t');    %数值函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 给各个函数赋值
f=188.7; T=1/f;  %翅拍频率 (Hz)和周期  % w =1185.6; 
t=linspace(0.0052824335,0.0052824335+3*T,1000);  % t_steady1   
dphi=dphi(w,t);         %(1*100)的行向量——拍打角速率                                                           % 输出
ddphi=ddphi(w,t);    %(1*100)的行向量——拍打角加速率                                                         % 输出
dpsi=dpsi(w,t);         %(1*100)的行向量——扭转角速率                                                             % 输出
ddpsi=ddpsi(w,t);     %(1*100)的行向量——扭转角加速率                                                         % 输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 拍打角和扭转角
% 拍打角—— (1*100)的行向量——拍打角——弧度制            %输出                                     
phi_pres=1.149*sin(w*t+1.571);
% 扭转角——(1*100)的行向量——扭转角——弧度制            %输出
psi_sim=-(0.005045+8.152*cos(t*w)+65.33*sin(t*w)+0.01177*cos(2*t*w)-0.006058*sin(2*t*w)-8.152*cos(3*t*w)+10.67*sin(3*t*w)+...
              0.03641*cos(4*t*w)-0.008574*sin(4*t*w)-1.666*cos(5*t*w)-0.8538*sin(5*t*w)+0.0234*cos(6*t*w)...
             -0.003305*sin(6*t*w)-0.1157*cos(7*t*w)+0.1072*sin(7*t*w)+0.02172*cos(8*t*w)-0.001624*sin(8*t*w))*pi/180;   
% [phi_min,k]=min(phi);  % 输出:  phi_min =-1.0157;  k =10756;      
% t_0=t(k);                       % 输出: t0 =0.0028;  
% figure(1)                                                          % 图1——拍打角和扭转角
% subplot(311)
% plot(t*f,psi_sim*180/pi,'g:',t*f,phi_pres*180/pi,'k-.','LineWidth',2) %转换为ms 和 度数degree   *10^3   *180/pi
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) vs \phi_{pres}(t)');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)');  
% title('被动转动角\psi_{sim}(t)和拍打角phi_{pres}(t)随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)随时间的变化
% axis([0.9,4.05,-105,105])
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
% axis([0.9,4.05,-inf,inf])
% Plot-拍打角加速率和扭转角角加速率
% subplot(313)
% plot(t/T,ddphi,'r-',t/T,ddpsi,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('角加速率 (rad/s^-2)')
% legend('\itdd\phi(t)','\itdd\psi(t)')
% title('拍打角加速率和扭转角加速率随时间的变化规律') 
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot-扭转角和几何攻角AOA
alpha1=pi/2-psi_sim.*sign(dphi);        %(1*200)的行向量——几何攻角——弧度制        %输出——全正几何攻角 
% Y = sign(X) returns an array Y the same size as X, where each element of Y is:
% *1 if the corresponding element of X is greater than zero
% * 0 if the corresponding element of X equals zero
% *-1 if the corresponding element of X is less than zero
% figure(2)                                                              % 图2—— 注意这里几何攻角始终取正值
% % % x_interval=[0,1/2,1/2,0];
% % % y_interval=[-100,-100,100,100];
% % % fill(x_interval,y_interval,'y');
% % % hold on
% % % legend('','\it\psi(t)','\it\alpha_1(t)')
% plot(t/T,psi_sim*180/pi,'b-',t/T,alpha1*180/pi,'g-')      %扭转角和几何攻角AOA随时间的变化规律
% xlabel('\itNormalized time')
% ylabel('\itAngle (°)')
% legend('\it\psi_{sim}(t)','\it\alpha_1(t)')
% % hold on
% % L=length(t);
% % plot([0,t(L)/T],[0,0],'k-');     %画x-axis
% title('扭转角和几何攻角alpha_1随时间的变化规律') 
% grid on
% axis([0.9,4.05,-105,105])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅膀坐标系下的角速度
% omega_x=dpsi;   % 扭转角:   x(1)=psi;  x(2)=dpsi;   
omega_y=dphi.*sin(psi_sim);
omega_z=dphi.*cos(psi_sim);
% omega_h=dphi;    % 翅膀坐标系下铰链角速率    % omega_h=sqrt(omega_y^2+omega_z^2);  % 右手法则顺时针
%% 翅膀坐标系下的角加速度——用于虚质量力的计算
% domega_x=ddpsi;
% domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi); 
% domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%% 气动攻角计算
v_y_nonr=omega_z;   % v_y=r*dphi*cos(psi)
v_z_nonr=-omega_y;    % v_z=-r*dphi*sin(psi)
alpha2=atan2(v_y_nonr,-v_z_nonr);   % 正确——注意与下文的alpha=atan2(omega_z,-omega_y)*180/pi; 不同
% 由于alpha=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %这里atan2给出象限正负值，尤其是alpha>pi/2时
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)                                               % 图3——使用翅膀坐标系下的角速度求解几何攻角
% plot(t/T,alpha1*180/pi,'r-',t/T,alpha2*180/pi,'b-','LineWidth',2)    
% xlabel('\itNormalized time')
% ylabel('\it\alpha_1 & \alpha_2 (deg)')
% legend('\alpha_1 (t)','\alpha_2 (t)')
% title('攻角随时间的变化规律')   % 攻角随时间的变化规律
% grid on
% axis([0.9,4.05,-inf,inf])
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
% (3) lift and drag coefficients with Alpha——随攻角变化的升阻力系数和法向气动力系数——第三种——全正系数
C_L2 =0.225+1.58*sin((2.13*abs(alpha)*180/pi-7.2)*pi/180);  % 几何攻角全正abs(alpha);——升力系数全正
C_D2 =1.92-1.55*cos((2.04*abs(alpha)*180/pi-9.82)*pi/180); % 几何攻角全正abs(alpha);——阻力系数全正
% C_N2=cos(abs(alpha)).*C_L2+sin(abs(alpha)).*C_D2;                % 几何攻角为正abs(alpha)——法向力系数全正
% 正负交替-由升阻力系数合成—— 2012-IEEE ASME TM-Veaceslav Arabagi 或2009-Science-Deng Xinyan
C_N2=sign(alpha2).*sqrt(C_L2.^2+C_D2.^2); 
% (4) Normal and tangential coefficients with Alpha——随攻角变化的法向和切向气动力系数——第四种——全正系数
C_N3=3.4*sin(abs(alpha));      % 2006-IEEE TR-Deng Xinyan 或 2010-EAE-RJ Wood-法向气动力系数输出
% 切向(tangential)气动力系数C_T只在alpha∈(-pi/4,pi/4)时不为零，其他情况为零
% if alpha<pi/4 & alpha>-pi/4
%     C_T=0.4*cos(2*alpha).^2
% else
%     C_T=0;
% end
C_T=0.4*(cos(2*alpha)).^2.*(alpha>-pi/4 & alpha<pi/4);       %切向气动力系数输出
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4)        % 图6—翅坐标系下—升阻力系数
% plot(t/T,C_L,'r-',t/T,C_D,'b-','LineWidth',2);      %时间正则化，乘以频率f, 或者除以周期T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of lift and drag \itvs. t \rm for flapping wing')
% title('气动升阻力系数随着时间的变化规律')
% legend('\itC_L','\itC_D')
% grid on
% axis([0.9,4.05,-inf,inf])
% figure(5)     % 图7—翅坐标系下—法向和切向气动力系数
% plot(t/T,C_N3,'r-',t/T,C_T,'b-')                                % 得到的法向和切向气动力系数需要核实
% xlabel('Normalized time')
% ylabel('C_N3(\alpha(t)) & C_T(\alpha (t))')
% legend('C_{N3}(\alpha)','C_T(\alpha)')
% title('法向和切向气动力系数随时间的变化规律')
% grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %三种法向气动力系数
% figure(6)     % 图8—翅坐标系下—法向气动力系数
% plot(t/T,C_N1,'r-',t/T,C_N2,'b-',t/T,C_N3,'g-','LineWidth',2);      %时间正则化，乘以频率f, 或者除以周期T;
% xlabel('\itt (Normalized time with flapping period)')
% ylabel('\itForce coefficients')
% % title('Coefficients of normal force \itvs. t \rm for flapping wing')
% title('法向气动力系数随时间的变化规律')
% legend('\itC_{N1}','\itC_{N2}','\itC_{N3}')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 翅运动规律输出和气动力系数输出
% wing_m_output=zeros(1000,12);
wing_m_output=[t',phi_pres',psi_sim',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

