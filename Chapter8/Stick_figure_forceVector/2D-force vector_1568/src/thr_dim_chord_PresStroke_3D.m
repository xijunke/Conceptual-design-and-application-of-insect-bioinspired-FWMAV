%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 绘制三维空间球棍模型图——受力图
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) 读取翅膀运动学-翅运动规律和几何攻角(AOA)等数据
clear all;clc;
%% 调用翅膀运动学-翅运动规律和几何攻角(AOA)等数据
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_simpres();      %调用函数kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,11)
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();
t=wing_kenimatics(:,1);                % 单位是ms
phi_pres=wing_kenimatics(:,2);        % 拍打角——单位是rad
psi_sim=wing_kenimatics(:,3);            % 拍打角——单位是rad
alpha_sim=wing_kenimatics(:,4);      % alpha=atan2(omega_z,-omega_y);  % 几何攻角——弧度制   %输出——有正有负
dphi_pres=wing_kenimatics(:,5);          % 单位是rad/s
dpsi_sim=wing_kenimatics(:,6);          % 单位是rad/s
ddphi_pres=wing_kenimatics(:,7);       % 单位是rad/s^2
ddpsi_sim=wing_kenimatics(:,8);     % 单位是rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
% C_N1=wing_kenimatics(:,11);   
% C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
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

%% (2) 截取单个稳定的周期——以攻角变化周期为区间
% 采用攻角为90度对应时刻为上下冲程的结束时间点
[phi_pres_max, locat_max]=max(phi_pres);   % phi_pres_max =1.149; locat_max =2;
T=1/f;
t_T=t(locat_max)+T;
index=find(t<=t_T);      % 2005
num=length(index);     % num=335;
% [alpha_sim_max, locat_max]=min(abs(alpha_sim));   % alpha_sim_max =0.3581; locat_max =219;
% T=1/f;
% t_T=t(locat_max)+T;
% index=find(t<=t_T);      % 552
% num=length(index);     % num=552;
% 刚好104个数据点=275-172+1;
t1=t((locat_max:num),1);                     % 单位是ms
psi1=psi_sim((locat_max:num),1);           % 扭转角——单位是rad
alpha1=alpha_sim((locat_max:num),1);         % 扭转角——单位是rad
dpsi1=dpsi_sim((locat_max:num),1);           % 扭转角——单位是rad
% ddpsi1=ddpsi_sim((locat_max:num),1);           % 扭转角——单位是rad
ddpsi1=dpsi_sim((locat_max:num),1);          % 扭转角——单位是rad—注意这里的扭转角加速度取了扭转角速度ddpsi1=dpsi_sim, 因为ODE方程没有输出ddpsi
phi1=phi_pres((locat_max:num),1);                        % 拍打角——单位是rad
dphi1=dphi_pres((locat_max:num),1);                        % 拍打角——单位是rad
ddphi1=ddphi_pres((locat_max:num),1);                        % 拍打角——单位是rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 将单个稳定的周期的总(合成)法向气动力数据读入——mN
% 以备thr_dim_chord2程序的调用用于绘制球棍图——翅坐标系下的法向力失终点y坐标
% B=[F_ytran,F_yadd1,F_yrot];    % size(B)  % (2000*1)     %——mN
F_normal=xlsread('Forcenormal_oneT_simpres.xlsx','A1:C1000');
% Forcenormal=0.2*F_normal((locat_max:num),1);   %  被缩小了1/5倍显示quiver3
F_ytran=0.1*F_normal((locat_max:num),1);
F_yadd=0.5*F_normal((locat_max:num),2);
F_yrot=0.2*F_normal((locat_max:num),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2)
% plot(t1*f,psi1*180/pi,'r-',t1*f,phi1*180/pi,'g-',t1*f,alpha1*180/pi,'k-',t1*f,2*abs(F_ytran/0.2),'b-','LineWidth',2)  % 转换为度数
% xlabel('\itNormalized time');  
% ylabel('\psi_{sim}(t) and \phi_{pres}(t) and \alpha_{sim}(t) and 2*abs(Force_{normal}(t))');  
% legend('\psi_{sim}(t)','\phi_{pres}(t)','\alpha_{sim}(t)','2*abs(Force_{normal}(t))');  
% title('被动转动角\psi_{sim}(t), 拍打角phi_{pres}(t) and 2*abs(Force_{normal}(t))随时间的变化')  
% grid on  % 被动转动角psi_sim(t)和扭转角phi_pres(t)随时间的变化
% % axis([1.6,2.7,-110,110])
% axis([0.9,2.1,-110,110])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) 将截取的单个稳定的周期——以攻角变化周期为区间的运动学数据输出到Excel表——弧度制
% 以备kenimatics_wing_and_AoA_fruitfly函数的调用用于计算气动力系数等
A=[t1,psi1,dpsi1,ddpsi1,phi1,dphi1,ddphi1];     % size(A)  % (2000*7)     % ——弧度制
xlswrite('wingmotion_oneT_simpres.xlsx',A,'sheet1','A1:G2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) 绘制单片翅弦微元的三维空间位置——冲程中点时刻
R_wingeff =3.0040;        % 单位是 mm   
xr=0.3289;                     % x-root offset  \mm
R_eff=R_wingeff+xr;
C_avereff =0.8854;          % 单位是 mm
C_025=0.25*C_avereff;
C_075=0.75*C_avereff;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % psi_sim=0;
% % alpha_sim=0;
% psi_sim=pi/4;
% alpha_sim=pi/4;
% % % 前缘——实心圆点绘制
% x1=R_eff;
% y1=C_025*sin(psi_sim);
% z1=C_025*cos(psi_sim);
% % 扭转轴心——实心圆点绘制
% x2=R_eff;
% y2=0;
% z2=0;
% % 弦向压心点————实心圆点绘制 
% d_cp=(0.82*abs(alpha_sim)/pi+0.05)*C_avereff;
% if d_cp>C_025
%     x3=R_eff;
%     y3=-(d_cp-C_025)*cos(psi_sim);
%     z3=-(d_cp-C_025)*sin(psi_sim);
% elseif d_cp<=C_025;
%     x3=R_eff;
%     y3=(C_025-d_cp)*sin(psi_sim);
%     z3=(C_025-d_cp)*cos(psi_sim);
% end
% % 中弦点—力失的起点—实心圆绘制
% x4=R_eff;
% y4=-C_025*sin(psi_sim);
% z4=-C_025*cos(psi_sim);
% % 后缘——实心圆点绘制
% x5=R_eff;
% y5=-C_075*sin(psi_sim);
% z5=-C_075*cos(psi_sim);
% % 片元弦向的特殊点阵
% x=[x1,x2,x3,x4,x5];
% y=[y1,y2,y3,y4,y5];
% z=[z1,z2,z3,z4,z5];
% wing_chord=[x;y;z];  % size(wing_chord)   % (3*5)
% % 绘制片元弦向的特殊点阵的分布
% figure(3)
% plot3(x1,y1,z1,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',8);   % 前缘——实心圆点绘制
% hold on
% plot3(x2,y2,z2,'ro','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);   % 扭转轴心——实心圆点绘制
% hold on
% plot3(x3,y3,z3,'bo','LineWidth',1.5,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',4);  % 弦向压心点————实心圆点绘制
% hold on
% plot3(x4,y4,z4,'go','LineWidth',1.5,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);  % 中弦点—力失的起点—实心圆绘制
% hold on
% plot3(x5,y5,z5,'ko','LineWidth',1.5,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3)   % 后缘——实心圆点绘制
% hold on
% plot3(x,y,z,'k-','LineWidth',2)
% xlabel('侧向(x)')
% ylabel('纵向(y)')
% zlabel('垂直方向(z)')
% grid on
% % axis square              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) 绘制单片翅弦微元的三维空间位置
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% animation——动画
% 动画片段程序开始段
Nstep=num- locat_max+1;    % Nstep=334; 
np=2;                     %  调节总的帧数和时间步长
% Nframe = np*10.4375;   % number of movie frames -------------------------- for change
Nframe = np*16.7;   % number of movie frames -------------------------- for change
Nskip = floor(Nstep/Nframe);
Nspeed = 2;  % number of frames per second ------------------------ for change
animation = 1;       %  1=save animation file --------------------------- for change
duration = Nframe/Nspeed;  % duration of movie  % 80 seconds @Nframe=4*40
display(['duratin of movie will be ', num2str(duration), ' seconds ']); %输出显示动画持续的时间(s)
%输出如下：duratin of movie will be 10 seconds 
if animation == 1
    aviobj = avifile('ElmentMotion_3D.avi','fps',Nspeed); % filename --- for change
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nskip=50;
% Nstep=2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nskip:Nstep   % L=[1:Nskip:Nstep]; % size(L)=1*40;
    %% 下面一段是坐标系之间的变换
    phi=phi1(i,1);            % 拍打角——单位是rad
    psi=psi1(i,1);             % 扭转角——单位是rad
    Yaw = [ cos(phi)   -sin(phi)   0; ...
                 sin(phi)   cos(phi)   0; ...
                 0              0            1];         % yaw matrix——符号可能需要方向
    Pitch = [1    0             0; ...
                 0    cos(psi)   -sin(psi); ...
                 0    sin(psi)    cos(psi)];        % pitch matrix——符号可能需要方向
    R_w2b = Yaw*Pitch;     % from wing frame to body frame——从翅坐标系到体坐标系
    %     R_b2w=R_w2b';                    % from body frame to wing frame——从体坐标系到翅坐标系
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 前缘——实心圆点绘制
    x1=R_eff;
    y1=0;
    z1=C_025;
    p_lead=[x1;y1;z1];
    p_lead=R_w2b*p_lead;
    % 扭转轴心——实心圆点绘制
    x2=R_eff;
    y2=0;
    z2=0;
    p_rotaxis=[x2;y2;z2];
    p_rotaxis=R_w2b*p_rotaxis;
    % 弦向压心点————实心圆点绘制 
    alpha=alpha1(i,1);             % 扭转角——单位是rad
    d_cp=(0.82*abs(alpha)/pi+0.05)*C_avereff;
    y_dcp=0;
    z_dcp=-d_cp;
%     if d_cp>C_025
%         x3=R_eff;
%         y3=y_dcp-y1;
%         z3=z_dcp-z1;
%     elseif d_cp<=C_025;
%         x3=R_eff;
%         y3=y_dcp+y1;
%         z3=z_dcp+z1;
%     end
    if d_cp>C_025
        x3=R_eff;
        y3=0;
        z3=-(d_cp-C_025);
    elseif d_cp<=C_025;
        x3=R_eff;
        y3=0;
        z3=(C_025-d_cp);
    end
   p_cop=[x3;y3;z3];
   p_cop=R_w2b*p_cop;
    % 中弦点—力失的起点—实心圆绘制
    x4=R_eff;
    y4=0;
    z4=-C_025;
    p_mid=[x4;y4;z4];
    p_mid=R_w2b*p_mid;
    % 后缘——实心圆点绘制
    x5=R_eff;
    y5=0;
    z5=-C_075;
    p_tail=[x5;y5;z5];
    p_tail=R_w2b*p_tail;
    % 片元弦向的特殊点阵
    x=[p_lead(1,1),p_rotaxis(1,1),p_cop(1,1),p_mid(1,1),p_tail(1,1)];
    y=[p_lead(2,1),p_rotaxis(2,1),p_cop(2,1),p_mid(2,1),p_tail(2,1)];
    z=[p_lead(3,1),p_rotaxis(3,1),p_cop(3,1),p_mid(3,1),p_tail(3,1)];
    wing_C=[x;y;z];  % size(wing_chord)   % (3*5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    theta=(-90:10:90)'*pi/180;   n=length(theta);     %  size(theta)=(19*1)
%    p1=[1:n;zeros(1,n)]';   % size(p1)=(19*2)            % y=0时x轴上分布的19个起始点的x坐标
%    r=2*ones(n,1);            % size(r)=(19*1)
%    % 极坐标转换为笛卡尔坐标
%    [u,v]=pol2cart(theta,r)  % 调用函数pol2cart:Transform polar or cylindrical coordinates to Cartesian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for j = 1:4
%         wing_C(:,j) = xl_B(i,:)'+bodyt(:,j);   % body shape in lab frame——可用于添加额外位矢
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
xlabel('侧向(x)')
ylabel('纵向(y)')
zlabel('垂直方向(z)')
grid on
axis([-0.5 3.5 -4 4 -1 1.5]);
% 绘制前缘实心圆点
plot3(p_lead(1,1),p_lead(2,1),p_lead(3,1),'ro-','LineWidth',1.5,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)
hold on
plot3([p_lead(1,1),p_tail(1,1)],[p_lead(2,1),p_tail(2,1)],[p_lead(3,1),p_tail(3,1)],'r-','LineWidth',1.5) % 绘制前缘后缘实线
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3(Fnormal (1,1:2),Fnormal (2,1:2),Fnormal (3,1:2),'b-','LineWidth',1.5) % 绘制作用于中弦点的法向力失
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) 使用arrow3绘制三维箭头
% % mid_chord=[(x_lead+x_tail)/2,(y_lead+y_tail)/2,(z_lead+z_tail)/2];  %作用于中弦点的法向力失——起点
% % F_normal_termini=[(x_lead+x_tail)/2,Forcenormal(i,1),(z_lead+z_tail)/2];  % 作用于中弦点的法向力失——终点
% x1=[(x_lead+x_tail)/2,(x_lead+x_tail)/2];
% y1=[(y_lead+y_tail)/2,-Forcenormal(i,1)];
% z1=[(z_lead+z_tail)/2,(z_lead+z_tail)/2];
% F_normal=[x1;y1;z1];  % size(wing_chord)   % (3*2)
% Fnormal = R_w2b*F_normal;        % wing_chord in body frame  % (3*2)
% F_mid=[Fnormal(1,1),Fnormal(2,1),Fnormal(3,1)];
% F_termini=[Fnormal(1,2),Fnormal(2,2),Fnormal(3,2)];
% % % pbaspect([4 4 1])    % 调用函数pbaspect: Set or query plot box aspect ratio(图箱纵横比)
% % % pbaspect('auto'); % daspect('auto')
% arrow3(F_mid,F_termini,'b',0.9)     % 调用函数arrow3绘制次作用于中弦点的法向力失
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 注意quiver的使用: 终点对应的是轴分量
% A three-dimensional quiver plot displays vectors with components (u,v,w) at the points (x,y,z), 
% where u,v,w,x,y, and z all have real (non-complex) values.
%% (a) 弦向压心点—绘制力失
x_cop=[x3,0];
y_cop=[y3,F_ytran(i,1)];  % 弦向压心点力失的起点——>终点
z_cop=[z3,0];
F_ytran_cop=[x_cop;y_cop;z_cop];  % size(wing_chord)   % (3*2)
Fytran_cop= R_w2b*F_ytran_cop;        % wing_chord in body frame  % (3*2)
quiver3(Fytran_cop(1,1),Fytran_cop(2,1),Fytran_cop(3,1),...                                                               % 起点
             Fytran_cop(1,2),Fytran_cop(2,2),Fytran_cop(3,2),0.5,'b-','LineWidth',2.0); % 比例因子0.5      % 终点
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (b) 中弦点—绘制力失
x_mid=[x4,0];
y_mid=[y4,F_yadd(i,1)];  % 中弦点力失的起点——>终点
z_mid=[z4,0];
F_yadd_mid=[x_mid;y_mid;z_mid];  % size(wing_chord)   % (3*2)
Fyadd_mid= R_w2b*F_yadd_mid;        % wing_chord in body frame  % (3*2)
quiver3(Fyadd_mid(1,1),Fyadd_mid(2,1),Fyadd_mid(3,1),...                                                               % 起点
             Fyadd_mid(1,2),Fyadd_mid(2,2),Fyadd_mid(3,2),0.5,'g-','LineWidth',2.0); % 比例因子0.5      % 终点
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 动画片段程序结束段
    if animation == 1
        frame = getframe;
        aviobj = addframe(aviobj,frame);
    else
        pause(0.5)
    end
end
if animation == 1
    aviobj = close(aviobj);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%