%% Stick figure force
clear all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第一部分——调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres(); %调用函数kenimatics_wing_and_AoA; % size(wing_kenimatics) % (1000,11)
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                        % 单位是ms
phi_pres=wing_kenimatics(:,2);            % 拍打角——单位是rad
psi_sim=wing_kenimatics(:,3);              % 拍打角——单位是rad
alpha_sim=wing_kenimatics(:,4);          % alpha=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi_pres=wing_kenimatics(:,5);          % 单位是rad/s
dpsi_sim=wing_kenimatics(:,6);            % 单位是rad/s
ddphi_pres=wing_kenimatics(:,7);        % 单位是rad/s^2
ddpsi_sim=wing_kenimatics(:,8);          % 单位是rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
f=188.7;                % Hz——翅拍频率 % T=1/f;  
% figure(1)
% plot(t*f,psi_sim*180/pi,'r-',t*f,phi_pres*180/pi,'g-',t*f,alpha_sim*180/pi,'b:','LineWidth',2)
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t)');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)');  
% title('被动转动角\psi_{sim}(t)和拍打角phi_{pres}(t)随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)随时间的变化
% axis([0.9,4.05,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 截取单个稳定的周期——以攻角变化周期为区间
% 采用攻角为90度对应时刻为上下冲程的结束时间点
[psi_sim_0, locat_0]=min(abs(psi_sim));   % psi_sim_0 =6.2058e-004; locat_0 =3;
T=1/f;
t_T=t(locat_0)+T;
index=find(t<=t_T);          % 336
num=length(index)+4;     % num=336+2;
% [alpha_sim_max, locat_0]=min(abs(alpha_sim));   % alpha_sim_max =0.3581; locat_max =219;
% T=1/f;
% t_T=t(locat_0)+T;
% index=find(t<=t_T);      % 552
% num=length(index);     % num=552;
% 刚好104个数据点=275-172+1;
t1=t((locat_0:num),1);                             % 单位是ms
psi1=psi_sim((locat_0:num),1);               % 扭转角——单位是rad
alpha1=alpha_sim((locat_0:num),1);       % 扭转角——单位是rad
dpsi1=dpsi_sim((locat_0:num),1);           % 扭转角——单位是rad
% ddpsi1=ddpsi_sim((locat_0:num),1);  % 扭转角——单位是rad
ddpsi1=dpsi_sim((locat_0:num),1);     % 扭转角——单位是rad—注意这里的扭转角加速度取了扭转角速度ddpsi1=dpsi_sim, 因为ODE方程没有输出ddpsi
phi1=phi_pres((locat_0:num),1);         % 拍打角——单位是rad
dphi1=dphi_pres((locat_0:num),1);      % 拍打角——单位是rad
ddphi1=ddphi_pres((locat_0:num),1);  % 拍打角——单位是rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) 将单个稳定的周期的总(合成)法向气动力数据读入——mN
% 以备thr_dim_chord2程序的调用用于绘制球棍图——翅坐标系下的法向力失终点y坐标
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %——mN
F_normal=xlsread('Forcenormal_oneT_simpres.xlsx','A1:C1000');
% Forcenormal=F_normal((locat_0:num),1);   %  被缩小了1/5=0.2倍显示quiver3
F_ytran=F_normal((locat_0:num),1);                % 10^(-1)*
F_yadd=F_normal((locat_0:num),2);
% F_yrot=F_normal((locat_0:num),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-',t1*f,2*abs(F_ytran/0.2),'b-','LineWidth',2)  % 转换为度数
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 2*abs(Force_{normal}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','2*abs(Force_{normal}(t))');  
% title('被动转动角\psi_{sim}(t), 拍打角phi_{pres}(t) and 2*abs(Force_{normal}(t))随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)随时间的变化
% % axis([1.6,2.7,-110,110])
% axis([0.9,2.1,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 第二部分—— 绘制单片翅弦微元的二维空间位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % R_wingeff =3.0040;        % 单位是 mm   
% % xr=0.3289;                     % x-root offset  \mm
% % R_eff=R_wingeff+xr;
C_avereff =0.8854;          % 单位是 mm
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  分解上下冲程
N_skip=12;  % 一个周期取28个点，半个周期取14个点
num=338;  % num=336+2;
k_index=(1:N_skip:num); %  length=length(k_index)  % (1*28)  
Nstep=num/2;                        % Nstep=169;
% k_index1=(1:12:Nstep);       % 1...13..25......169  % k_index1(1,1:15);  
% N_up=length(k_index1);      % (15)
% k_index2=(Nstep:12:num);  % 169...181...193......337  % k_index2(1,15:29);
% N_down=length(k_index2); % (15)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 上冲程
k_index1=(1:12:Nstep);         % 
N_up=length(k_index1);        % (15)
% t_up=t1(k_index1);    % size(t_up)  %  (15*1)
psi_up=psi1(k_index1);  % 针对一个周期之内的数据，每隔8个点取一个数据形成向量
alpha_up=alpha1(k_index1);
% figure(3)
% plot(t_up*f,psi_up*180/pi,'rd','LineWidth',2)      % 转换为度数
% Forcenormalup=Forcenormal(k_index1);
F_ytranup=F_ytran(k_index1);  % size(F_ytranup)   % (15*1)     %  10^(-5)*
F_yaddup=F_yadd(k_index1);   % size(F_yaddup)   % (15*1)
% 扭转轴点坐标——人为设定等时间距离
Y_axisup=linspace(-2.5,2.5,15); 
Z_axisup=zeros(1,length(Y_axisup));
% Z_axisup=zeros(1,length(Y_axisup))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 翅拍上冲程——不考虑扭转
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:1:N_up   % N_up=15
%     % (1) 前缘点坐标——实心圆点绘制
%    y1=Y_axisup(1,i);
%    z1=C_025+Z_axisup(1,i);
%     % (2) 弦向压心点————实心圆点绘制 
%     alpha=pi/2;             % 扭转角——单位是rad
%     d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;        
%     if d_cp>C_025
%         y3=Y_axisup(1,i);
%         z3=-(d_cp-C_025)+Z_axisup(1,i);
%     elseif d_cp<=C_025;
%         y3=Y_axisup(1,i);
%         z3=(C_025-d_cp)+Z_axisup(1,i);
%     end
%     % (3) 中弦点坐标—力失的起点—实心圆绘制
%     y4=Y_axisup(1,i);
%     z4=-C_025+Z_axisup(1,i);
%     % (4) 后缘点坐标——实心圆点绘制
%     y5=Y_axisup(1,i);
%     z5=-C_075+Z_axisup(1,i);
%     % 绘制片元弦向的特殊点阵的分布
%     figure(4)
%     plot([y1,y5],[z1,z5],'k-','LineWidth',1.5)          %  绘制直线——连接前后缘
%     hold on
%     plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
%     hold on
%     plot(y1,z1,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);   % 前缘——实心圆点绘制
%     hold on
%     plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
%     hold on
%     % plot(y5,z5,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',1.5)   % 后缘——实心圆点绘制
%     % hold on
%     % (4) 弦向压心点—绘制力失       
%        % 弦向压心点力失的起点
%     plot(y3,z3,'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);  % 弦向压心点————实心圆点绘制
%     hold on
% %     % p_cop=p_coprot+displacement;  
% %     % displacement=[Y_axisup(1,i);Z_axisup(1,i)];  % (2*1)
% %     y_cop=-F_ytranup(i,1)+Y_axisup(1,i);    %——注意符号
% % %     z_cop=z3+Z_axisup(1,i);
% %      z_cop=0;
% %     F_ytran_cop0=[y_cop;z_cop];  % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
% % %     Fytran_coprot= R_w2s*F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
% %     Fytran_coprot= F_ytran_cop0;        % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
% %     Fytran_cop=Fytran_coprot; 
% %     quiver(y3,z3,...                                                                 % 起点
% %                Fytran_cop(1,1),Fytran_cop(2,1),0.5,'b-','LineWidth',2.0); % 比例因子0.5        % 终点
% %     hold on
% end
% xlabel('纵向(y)')
% ylabel('垂直方向向(z)')
% grid on
% axis([-4,4,-3,3])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% animation———————— 动画片段程序开始段——(1)
N_step=N_up;         % N_up=15; 
animation=1;         %  1=save animation file --------------------------- for change
np=1;                    %  调节总的帧数和时间步长
Nframe=np*15;     % number of movie frames -------------------------- for change
Nskip1=floor(N_step/Nframe);
Nspeed=3;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 5 seconds @Nframe=1*15
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
%输出如下: duratin of movie will be 5 seconds 
if animation == 1
    aviobj = avifile('Stick_figure_force_upstroke.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_up   % N_up=15
    % 下面一段是坐标系之间的变换
    psi=-psi_up(i,1);             % 扭转角——单位是rad    %——注意符号
    Pitch = [cos(psi)   -sin(psi); ...
                 sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向
    R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    displacement=[Y_axisup(1,i);Z_axisup(1,i)];  % (2*1)
   % (1) 前缘点坐标——实心圆点绘制
    y1=0;
    z1=C_025;
    p_lead0=[y1;z1];                   % (2*1)
    p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)
    p_lead=p_leadrot+displacement;          % (2*1)
   % (4) 后缘点坐标——实心圆点绘制
    y5=0;
    z5=-C_075;
    p_tail0=[y5;z5];                   % (2*1)
    p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
    p_tail=p_tailrot+displacement;          % (2*1)   
    % (2) 弦向压心点————实心圆点绘制 
    % alpha_up=alpha1(i,1);             % 扭转角——单位是rad
    d_cp=(0.82*abs(alpha_up(i,1))/pi+0.05)*C_avereff;        
    if d_cp>C_025
        y3=0;
        z3=-(d_cp-C_025);
    elseif d_cp<=C_025;
        y3=0;
        z3=C_025-d_cp;
    end
    p_cop0=[y3;z3];
    p_coprot=R_w2s*p_cop0;
    p_cop=p_coprot+displacement; 
    % (3) 中弦点坐标—力失的起点—实心圆绘制
    y4=0;
    z4=-C_025;
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % 片元弦向的特殊点阵
%     y=[p_lead(1,1),Y_axisup(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisup(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(5)
    plot([Y_axisup(1,1),Y_axisup(1,N_up)],[Z_axisup(1,1),Z_axisup(1,N_up)],'r:','LineWidth',1.5) % 扭转轴心——实心圆点绘制
    hold on
    xlabel('纵向(y)')
    ylabel('垂直方向向(z)')
    grid on
    axis([-4,4,-1,3.1])
    hold on
    % (1) 绘制直线——连接前后缘
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',2) % 绘制前缘后缘实线
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (2) 扭转轴心——实心圆点绘制
    plot(Y_axisup(1,i),Z_axisup(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (3) 绘制前缘实心圆点
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (4) 弦向压心点—绘制力失
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % 弦向压心点————实心圆点绘制
    hold on   
    % p_cop=p_coprot+displacement;                 % 弦向压心点力失的起点
    % displacement=[Y_axisup(1,i);Z_axisup(1,i)];          % (2*1)
    y_cop=-F_ytranup(i,1)+Y_axisup(1,i);                     %——注意符号
    z_cop=z3+Z_axisup(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % 起点
               Fytran_cop(1,1),Fytran_cop(2,1),0.14,'b-','LineWidth',1.5);    % 终点   % 比例因子0.05——————注意力失被缩小了0.14倍 
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % (5) 中弦点—绘制力失
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
   % hold on  
    y_mid=-F_yaddup(i,1)+Y_axisup(1,i);  % 中弦点力失的起点——>终点
    z_mid=z4+Z_axisup(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % 起点
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.135,'g-','LineWidth',1.5);   % 终点   % 比例因子0.05——————注意力失被缩小了0.135倍 
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 动画片段程序结束段——(2)
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  动画片段程序结束段——(3)
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) 下冲程
k_index2=(Nstep:12:num);     % 169...181...193......337
N_down=length(k_index2);    % (15)
% t_down=t1(k_index2);    %  size(t_down)  %  (15*1)
psi_down=psi1(k_index2);  % 针对一个周期之内的数据，每隔个8个点取一个数据形成向量
alpha_down=alpha1(k_index2);
% figure(6)
% plot(t_down*f,psi_down*180/pi,'rd','LineWidth',2)      % 转换为度数
% Forcenormaldown=Forcenormal(k_index2);
F_ytrandown=F_ytran(k_index2);  % size(F_ytrandown)   % (15*1)     %  10^(-5)*
F_yadddown=F_yadd(k_index2);   % size(F_yadddown)   % (15*1)
% 扭转轴点坐标——人为设定等时间距离
% Y_axisdown=linspace(-2.5,2.5,15); 
% Z_axisdown=zeros(1,length(Y_axisdown));
Y_axisdown=fliplr(Y_axisup); 
Z_axisdown=fliplr(Z_axisup);
%%  Z_axisdown=zeros(1,length(Y_axisdown))+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % animation———————— 动画片段程序开始段——(1)
Nstep=N_down;    % N_down=15; 
animation=1;       %  1=save animation file --------------------------- for change
np=1;                     %  调节总的帧数和时间步长
Nframe=np*15;   % number of movie frames -------------------------- for change
Nskip2=floor(Nstep/Nframe);
Nspeed=3;  % number of frames per second ------------------------ for change
duration=Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
%输出如下: duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('Stick_figure_force_downstroke.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:1:N_down  % N_down=15
    % 下面一段是坐标系之间的变换
    psi=psi_down(i,1);             % 扭转角——单位是rad
    Pitch = [cos(psi)   sin(psi); ...
                 -sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向
    R_w2s =Pitch;     % from wing frame to spar frame——从翅坐标系到体坐标系
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];  % (2*1)
    % (1) 前缘点坐标——实心圆点绘制
    y1=0;
    z1=C_025;
    p_lead0=[y1;z1];                   % (2*1)
    p_leadrot=R_w2s*p_lead0;    % (2*2)*(2*1)=(2*1)
    p_lead=p_leadrot+displacement;          % (2*1)
    % (4) 后缘点坐标——实心圆点绘制
    y5=0;
    z5=-C_075;
    p_tail0=[y5;z5];                   % (2*1)
    p_tailrot=R_w2s*p_tail0;      % (2*2)*(2*1)=(2*1)
    p_tail=p_tailrot+displacement;          % (2*1)
    % (2) 弦向压心点————实心圆点绘制 
    %  alpha_down=alpha1(i,1);             % 扭转角——单位是rad
    d_cp=(0.82*abs(alpha_down(i,1))/pi+0.05)*C_avereff;        
    if d_cp>C_025
        y3=0;
        z3=-(d_cp-C_025);
    elseif d_cp<=C_025;
        y3=0;
        z3=C_025-d_cp;
    end
    p_cop0=[y3;z3];
    p_coprot=R_w2s*p_cop0;
    p_cop=p_coprot+displacement; 
    % (3) 中弦点坐标—力失的起点—实心圆绘制
    y4=0;
    z4=-C_025;
    p_mid0=[y4;z4];
    p_midrot=R_w2s*p_mid0;
    p_mid=p_midrot+displacement; 
    % 片元弦向的特殊点阵
%     y=[p_lead(1,1),Y_axisdown(1,i),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
%     z=[p_lead(2,1),Z_axisdown(1,i),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
%     wing_C=[y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(7)
    plot([Y_axisdown(1,1),Y_axisdown(1,N_down-1)],[Z_axisdown(1,1),Z_axisdown(1,N_down-1)],'r:','LineWidth',1.5) % 扭转轴心——实心圆点绘制
    hold on
    xlabel('纵向(y)')
    ylabel('垂直方向向(z)')
    grid on
    axis([-4,4,-1,3.1])
    hold on
    % (1) 绘制直线——连接前后缘
    plot([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],'k-','LineWidth',1.5) % 绘制前缘后缘实线
    hold on
    % (2) 扭转轴心——实心圆点绘制
    plot(Y_axisdown(1,i),Z_axisdown(1,i),'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3) % 扭转轴心——实心圆点绘制
    hold on
     % (3) 绘制前缘实心圆点
    plot(p_lead(1,1),p_lead(2,1),'ko-','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
    hold on
   % (4) 弦向压心点—绘制力失
    plot(p_cop(1,1),p_cop(2,1),'bo','LineWidth',1.7,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  % 弦向压心点————实心圆点绘制
    hold on
    % p_cop=p_coprot+displacement;                 % 弦向压心点力失的起点
    % displacement=[Y_axisdown(1,i);Z_axisdown(1,i)];          % (2*1)
    y_cop=-F_ytrandown(i,1)+Y_axisdown(1,i);                     %——注意符号
    z_cop=z3+Z_axisdown(1,i);
    F_ytran_cop0=[y_cop;z_cop];                           % size(wing_chord)   % (2*2)    % 弦向压心点力失的终点
    Fytran_coprot= R_w2s*F_ytran_cop0;              % wing_chord in body frame  % (2*2)     %—R_w2s'—注意符号
    Fytran_cop=Fytran_coprot; 
    quiver(p_cop(1,1),p_cop(2,1),...                        % 起点
               Fytran_cop(1,1),Fytran_cop(2,1),0.14,'b-','LineWidth',1.5);    % 终点   % 比例因子0.05——————注意力失被缩小了0.14倍 
    hold on
   % (5) 中弦点—绘制力失
   % plot(y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);  % 中弦点—力失的起点—实心圆绘制
   % hold on  
    y_mid=-F_yadddown(i,1)+Y_axisdown(1,i);  % 中弦点力失的起点——>终点
    z_mid=z4+Z_axisdown(1,i);
    F_yadd_mid0=[y_mid;z_mid];  % size(wing_chord)   % (2*2)
    Fyadd_midrot=R_w2s*F_yadd_mid0;        % wing_chord in body frame  % (2*2)
    Fyadd_mid=Fyadd_midrot; 
    quiver(p_mid(1,1),p_mid(2,1),...                                                      % 起点
               Fyadd_mid(1,1),Fyadd_mid(2,1),0.135,'g-','LineWidth',1.5);   % 终点   % 比例因子0.05——————注意力失被缩小了0.135倍 
    hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 动画片段程序结束段——(2)
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  动画片段程序结束段——(3)
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
